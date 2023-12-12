from __future__ import annotations
from pypath.share import curl, settings
from contextlib import ExitStack

from bioregistry import normalize_curie
from tqdm.notebook import tqdm
from time import time

from pypath.inputs import expasy, uniprot

class EC:
    def __init__(self, add_prefix = True):
        self.add_prefix = add_prefix
        self.swissprots = list(uniprot._all_uniprots("*", True))
        
    def download_ec_data(
        self,
        cache=False,
        debug=False,
        retries=3,
    ):
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
                
            print("Started downloading Expasy EC number data")
            t0 = time()
            
            self.enzymes = expasy.enzymes()
            
            self.enzyme_classes = expasy.enzyme_classes()
            
            t1 = time()
            print(f"Expasy EC number data is downloaded in {round((t1-t0) / 60, 2)} mins")
            
    def get_nodes(self, label= "ec_number"):
        if not hasattr(self, "enzymes") or not hasattr(self, "enzyme_classes"):
            self.download_ec_data()
        if not hasattr(self, "ec_dict"):
            self.prepare_ec_hierarchy_dict()
        
        print("Started writing ec number nodes")
        
        node_list = []
        for level_1_entry, level_1_dict in tqdm(self.ec_dict.items()):
            level_1_id = self.add_prefix_to_id(prefix="eccode", identifier=level_1_entry)
            node_list.append((level_1_id, label, {"name":level_1_dict["name"].replace("|",",").replace("'","^")},))

            for level_2_entry, level_2_dict in level_1_dict.items():
                if level_2_entry != "name":
                    level_2_id = self.add_prefix_to_id(prefix="eccode", identifier=level_2_entry)
                    node_list.append((level_2_id, label, {"name":level_2_dict["name"].replace("|",",").replace("'","^")},))

                    for level_3_entry, level_3_dict in level_2_dict.items():
                        if level_3_entry != "name":
                            level_3_id = self.add_prefix_to_id(prefix="eccode", identifier=level_3_entry)
                            node_list.append((level_3_id, label, {"name":level_3_dict["name"].replace("|",",").replace("'","^")},))

                            if level_3_dict["entries"]:                        
                                for level_4_entry in level_3_dict["entries"]:
                                    level_4_id = self.add_prefix_to_id(prefix="eccode", identifier=level_4_entry)
                                    node_list.append((level_4_id, label,
                                                      {"name":self.enzymes[level_4_entry]["name"][0].replace(".","").replace("|",",").replace("'","^")},))
                                    
        return node_list
    
    def get_edges(self):
        
        edge_list = []
        
        edge_list.extend(self.get_protein_ec_edges())
        
        edge_list.extend(self.get_ec_hierarchy_edges())
        
        return edge_list
    
    def get_ec_hierarchy_edges(self, label="ec_number_is_a_ec_number"):
        if not hasattr(self, "enzymes") or not hasattr(self, "enzyme_classes"):
            self.download_ec_data()
        if not hasattr(self, "ec_dict"):
            self.prepare_ec_hierarchy_dict()
            
        print("Started writing ec number hierarchical edges")
        
        edge_list = []
        for level_1_entry, level_1_dict in tqdm(self.ec_dict.items()):
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
                                    
        return edge_list
    
    def get_protein_ec_edges(self, label="protein_catalyzes_ec_number"):
        if not hasattr(self, "enzymes"):
            self.download_ec_data()
        
        print("Started writing protein-ec number edges")
        
        edge_list = []
        for ec_number, ec_number_items in tqdm(self.enzymes.items()):
            if ec_number_items.get("annotations"):
                for protein in ec_number_items["annotations"]:
                    if protein in self.swissprots:
                        protein_id = self.add_prefix_to_id(prefix="uniprot", identifier=protein)
                        ec_id = self.add_prefix_to_id(prefix="eccode", identifier=ec_number)
                        edge_list.append((None, protein_id, ec_id, label, {}))
                        
        return edge_list                        
            
    def prepare_ec_hierarchy_dict(self):
        splitted = str(self.enzyme_classes).split("\n")
        
        self.ec_dict = {}

        for entry in splitted:
            if entry:        
                if "\t" in entry and entry.count("\t") == 1:
                    level_2_entry = level_1_entry.split(".")[0] + "." + entry.replace("\t", "").split(":")[0].strip() + ".-.-"
                    level_2_name = entry.replace("\t", "").split(":")[1].strip()
                    self.ec_dict[level_1_entry][level_2_entry] = {"name":level_2_name}
                elif "\t" in entry and entry.count("\t") == 2:
                    level_3_entry = ".".join(level_2_entry.split(".")[:2]) + "." + entry.replace("\t", "").split(":")[0].strip() + ".-"
                    level_3_name = entry.replace("\t", "").split(":")[1].strip()
                    search = ".".join(level_2_entry.split(".")[:2]) + "." + entry.replace("\t", "").split(":")[0].strip()
                    self.ec_dict[level_1_entry][level_2_entry][level_3_entry] = {"name":level_3_name, "entries":[]}

                    for level_4 in self.enzymes.keys():
                        if (not self.enzymes[level_4]["name"][0].startswith("Transferred entry") and search == ".".join(level_4.split(".")[:3]))\
                        and (not self.enzymes[level_4]["name"][0].startswith("Deleted") and search == ".".join(level_4.split(".")[:3])):
                            self.ec_dict[level_1_entry][level_2_entry][level_3_entry]["entries"].append(level_4)
                else:
                    level_1_entry = entry.split(":")[0].strip() + ".-.-.-"
                    level_1_name = entry.split(":")[1].strip()
                    self.ec_dict[level_1_entry] = {"name":level_1_name}
                    
    def add_prefix_to_id(self, prefix=None, identifier=None, sep=":") -> str:
        """
        Adds prefix to database id
        """
        if self.add_prefix and identifier:
            return normalize_curie(prefix + sep + str(identifier))
        
        return identifier
