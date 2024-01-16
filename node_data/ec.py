from __future__ import annotations
from pypath.share import curl, settings
from contextlib import ExitStack

from bioregistry import normalize_curie
from tqdm import tqdm
from time import time

import pandas as pd
import os

from enum import Enum, auto
from typing import Union

from pypath.inputs import expasy, uniprot

from biocypher._logger import logger

from pydantic import BaseModel, DirectoryPath, validate_call

logger.debug(f"Loading module {__name__}.")

class ECNodeField(Enum):
    NAME = "name"

class ECEdgeType(Enum):
    EC_HIERARCHY = auto()
    PROTEIN_TO_EC = auto()

class ECModel(BaseModel):
    ec_node_fields: Union[list[ECNodeField], None] = None
    edge_types: Union[list[ECEdgeType], None] = None
    test_mode: bool = False
    export_csv: bool = False
    output_dir: DirectoryPath | None = None
    add_prefix: bool = True

class EC:
    def __init__(self,
                 ec_node_fields: Union[list[ECNodeField], None] = None,
                 edge_types: Union[list[ECEdgeType], None] = None,
                 test_mode: bool = False,
                 export_csv: bool = False,
                 output_dir: DirectoryPath | None = None,
                 add_prefix: bool = True):
        
        model = ECModel(ec_node_fields=ec_node_fields,
                        edge_types=edge_types,
                        test_mode=test_mode,
                        export_csv=export_csv,
                        output_dir=output_dir,
                        add_prefix=add_prefix).model_dump()
        
        self.export_csv = model["export_csv"]
        self.output_dir = model["output_dir"]
        self.add_prefix = model["add_prefix"]
        self.swissprots = set(uniprot._all_uniprots("*", True))

        # set node fields
        self.set_node_fields(ec_node_fields=model["ec_node_fields"])

        # set edge types
        self.set_edge_types(edge_types=model["edge_types"])

        # set early_stopping, if test_mode true
        self.early_stopping = None
        if model["test_mode"]:
            self.early_stopping = 100
    
    @validate_call
    def download_ec_data(
        self,
        cache: bool = False,
        debug: bool = False,
        retries: int = 3,
    ) -> None:
        """
        Wrapper function to download ec data from various databases using pypath.
        Args
            cache: if True, it uses the cached version of the data, otherwise
            forces download.
            debug: if True, turns on debug mode in pypath.
            retries: number of retries in case of download error.
        """
        
        with ExitStack() as stack:
            stack.enter_context(settings.context(retries=retries))

            if debug:
                stack.enter_context(curl.debug_on())

            if not cache:
                stack.enter_context(curl.cache_off())
                
            logger.debug("Started downloading Expasy EC number data")
            t0 = time()
            
            self.enzymes = expasy.expasy_enzymes()
            
            self.enzyme_classes = expasy.expasy_enzyme_classes()

            self.prepare_ec_hierarchy_dict()
            
            t1 = time()
            logger.info(f"Expasy EC number data is downloaded in {round((t1-t0) / 60, 2)} mins")

    @validate_call    
    def get_nodes(self, label: str = "ec_number") -> list[tuple]:
        if not hasattr(self, "enzymes") or not hasattr(self, "enzyme_classes"):
            self.download_ec_data()
        if not hasattr(self, "ec_dict"):
            self.prepare_ec_hierarchy_dict()
        
        logger.info("Started writing ec number nodes")
        
        node_list = []
        
        for index, (level_1_entry, level_1_dict) in tqdm(enumerate(self.ec_dict.items())):
            level_1_id = self.add_prefix_to_id(prefix="eccode", identifier=level_1_entry)
            props = {}
            if ECNodeField.NAME in self.ec_node_fields:
                props[ECNodeField.NAME.value] = level_1_dict["name"].replace("|",",").replace("'","^")
            
            node_list.append((level_1_id, label, props))

            for level_2_entry, level_2_dict in level_1_dict.items():
                if level_2_entry != "name":
                    level_2_id = self.add_prefix_to_id(prefix="eccode", identifier=level_2_entry)
                    props = {}
                    if ECNodeField.NAME in self.ec_node_fields:
                        props[ECNodeField.NAME.value] = level_2_dict["name"].replace("|",",").replace("'","^")

                    node_list.append((level_2_id, label, props))

                    for level_3_entry, level_3_dict in level_2_dict.items():
                        if level_3_entry != "name":
                            level_3_id = self.add_prefix_to_id(prefix="eccode", identifier=level_3_entry)
                            props = {}
                            if ECNodeField.NAME in self.ec_node_fields:
                                props[ECNodeField.NAME.value] = level_3_dict["name"].replace("|",",").replace("'","^")
                            
                            node_list.append((level_3_id, label, props))

                            if level_3_dict["entries"]:                        
                                for level_4_entry in level_3_dict["entries"]:
                                    level_4_id = self.add_prefix_to_id(prefix="eccode", identifier=level_4_entry)
                                    props = {}
                                    if ECNodeField.NAME in self.ec_node_fields:
                                        props[ECNodeField.NAME.value] = self.enzymes[level_4_entry]['de'].replace(".","").replace("|",",").replace("'","^")
                                    
                                    node_list.append((level_4_id, label, props))
            
            if self.early_stopping and index+1 == self.early_stopping:
                break

        # write ec node data to cv
        if self.export_csv:
            if self.output_dir:
                full_path = os.path.join(self.output_dir, "Ec.csv")
            else:
                full_path = os.path.join(os.getcwd(), "Ec.csv")

            df_list = []
            for _id, _, props in node_list:
                row = {"ec_number":_id} | props
                df_list.append(row)

            df = pd.DataFrame.from_records(df_list)
            df.to_csv(full_path, index=False)
            logger.info(f"EC node data is written: {full_path}")

        return node_list
    
    def get_edges(self) -> list[tuple]:
        
        edge_list = []
        
        if ECEdgeType.PROTEIN_TO_EC in self.edge_types:
            edge_list.extend(self.get_protein_ec_edges())
        
        if ECEdgeType.EC_HIERARCHY in self.edge_types:
            edge_list.extend(self.get_ec_hierarchy_edges())
        
        return edge_list
    
    @validate_call
    def get_ec_hierarchy_edges(self, label: str = "ec_number_is_a_ec_number") -> list[tuple]:
        if not hasattr(self, "enzymes") or not hasattr(self, "enzyme_classes"):
            self.download_ec_data()
        if not hasattr(self, "ec_dict"):
            self.prepare_ec_hierarchy_dict()
            
        logger.info("Started writing ec number hierarchical edges")
        
        edge_list = []
        for index, (level_1_entry, level_1_dict) in tqdm(enumerate(self.ec_dict.items())):
            level_1_id = self.add_prefix_to_id(prefix="eccode", identifier=level_1_entry)

            for level_2_entry, level_2_dict in level_1_dict.items():
                if level_2_entry != "name":
                    level_2_id = self.add_prefix_to_id(prefix="eccode", identifier=level_2_entry)
                    
                    edge_list.append((None, level_2_id, level_1_id, label, {}))

                    for level_3_entry, level_3_dict in level_2_dict.items():
                        if level_3_entry != "name":
                            level_3_id = self.add_prefix_to_id(prefix="eccode", identifier=level_3_entry)
                            
                            edge_list.append((None, level_2_id, level_3_id, label, {}))

                            if level_3_dict["entries"]:                        
                                for level_4_entry in level_3_dict["entries"]:
                                    level_4_id = self.add_prefix_to_id(prefix="eccode", identifier=level_4_entry)
                                    
                                    edge_list.append((None, level_4_id, level_3_id, label, {}))

            if self.early_stopping and index+1 == self.early_stopping:
                break

        # write ec hierarchy data to csv
        if self.export_csv:
            if self.output_dir:
                full_path = os.path.join(self.output_dir, "Ec_hierarchy.csv")
            else:
                full_path = os.path.join(os.getcwd(), "Ec_hierarchy.csv")

            df_list = []
            for _, child, parent, label, _ in edge_list:
                df_list.append({"child_id":child, "parent_id":parent, "label":label})

            df = pd.DataFrame.from_records(df_list)
            df.to_csv(full_path, index=False)
            logger.info(f"EC hierarchy edge data is written: {full_path}")

        return edge_list
    
    @validate_call
    def get_protein_ec_edges(self, label: str = "protein_catalyzes_ec_number") -> list[tuple]:
        if not hasattr(self, "enzymes"):
            self.download_ec_data()
        
        logger.info("Started writing protein-ec number edges")
        
        edge_list = []
        for index, (ec_number, ec_number_items) in tqdm(enumerate(self.enzymes.items())):
            if ec_number_items.get("uniprots"):
                for protein in ec_number_items["uniprots"]:
                    if protein in self.swissprots:
                        protein_id = self.add_prefix_to_id(prefix="uniprot", identifier=protein)
                        ec_id = self.add_prefix_to_id(prefix="eccode", identifier=ec_number)
                        edge_list.append((None, protein_id, ec_id, label, {}))

            if self.early_stopping and index+1 == self.early_stopping:
                break


        # write protein-ec data to csv
        if self.export_csv:
            if self.output_dir:
                full_path = os.path.join(self.output_dir, "Protein_to_ec.csv")
            else:
                full_path = os.path.join(os.getcwd(), "Protein_to_ec.csv")

            df_list = []
            for _, protein_id, ec_id, _, _ in edge_list:
                df_list.append({"protein_id":protein_id, "ec_id":ec_id})

            df = pd.DataFrame.from_records(df_list)
            df.to_csv(full_path, index=False)
            logger.info(f"Protein-ec edge data is written: {full_path}")

        return edge_list                        
            
    def prepare_ec_hierarchy_dict(self) -> None:

        if not hasattr(self, "enzyme_classes"):
            self.enzyme_classes = expasy.expasy_enzyme_classes()

        if not hasattr(self, "enzymes"):
            self.enzymes = expasy.expasy_enzymes()
        
        logger.debug("Started preparing ec hierarchy dictionary")

        self.ec_dict = {}

        for entry, name in self.enzyme_classes:
            entry = entry.replace(" ", "")
            
            # if there is 3 - in the entry, it is a level 1 entry
            if entry.count("-") == 3:
                self.ec_dict[entry] = {"name":name}
            # if there is 2 - in the entry, it is a level 2 entry
            elif entry.count("-") == 2:
                level_1_entry = entry.split(".")[0] + ".-.-.-"
                if level_1_entry not in self.ec_dict:
                    self.ec_dict[level_1_entry] = {}
                self.ec_dict[level_1_entry][entry] = {"name":name}
            # if there is 1 - in the entry, it is a level 3 entry
            elif entry.count("-") == 1:
                level_1_entry = entry.split(".")[0] + ".-.-.-"
                level_2_entry = level_1_entry.split(".")[0] + "." + entry.split(".")[1] + ".-.-"
                if level_1_entry not in self.ec_dict:
                    self.ec_dict[level_1_entry] = {}
                if level_2_entry not in self.ec_dict[level_1_entry]:
                    self.ec_dict[level_1_entry][level_2_entry] = {}
                self.ec_dict[level_1_entry][level_2_entry][entry] = {"name":name, 'entries':[]}
        
        for level_4 in self.enzymes.keys():
            level_1_entry = level_4.split(".")[0] + ".-.-.-"
            level_2_entry = level_1_entry.split(".")[0] + "." + level_4.split(".")[1] + ".-.-"
            level_3_entry = level_2_entry.split(".")[0] + "." + level_4.split(".")[1] + "." + level_4.split(".")[2] + ".-"
            if not self.enzymes[level_4]["de"].startswith("Transferred entry") and not self.enzymes[level_4]["de"].startswith("Deleted"):
                self.ec_dict[level_1_entry][level_2_entry][level_3_entry]["entries"].append(level_4)

    @validate_call           
    def add_prefix_to_id(self, prefix: str =None, identifier: str = None, sep: str = ":") -> str:
        """
        Adds prefix to database id
        """
        if self.add_prefix and identifier:
            return normalize_curie(prefix + sep + str(identifier))
        
        return identifier
    
    def set_node_fields(self, ec_node_fields):
        if ec_node_fields:
            self.ec_node_fields = ec_node_fields
        else:
            self.ec_node_fields = [field.value for field in ECNodeField]
    
    def set_edge_types(self, edge_types):
        if edge_types:
            self.edge_types = edge_types
        else:
            self.edge_types = [field for field in ECEdgeType]
