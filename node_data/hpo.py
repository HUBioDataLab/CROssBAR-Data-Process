from __future__ import annotations

from pypath.share import curl, settings

from pypath.inputs import hpo, ontology

from contextlib import ExitStack

from bioregistry import normalize_curie
from tqdm.notebook import tqdm
from time import time

from enum import Enum, auto

import pandas as pd
import numpy as np

class HPO:
    def __init__(self, add_prefix = True):
        self.add_prefix = add_prefix
    
    def download_hpo_data(
        self,
        cache=False,
        debug=False,
        retries=3,
    ):
        """
        Wrapper function to download hpo data from various databases using pypath.
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
                
            print("Started downloading HPO data")
            t0 = time()
        
            self.protein_hpo_annotations = hpo.hpo_annotations()
            
            self.hpo_terms = hpo.hpo_terms()
            
            self.hpo_phenotype_disease = hpo.hpo_diseases()
            
            self.hpo_ontology = hpo.hpo_ontology()
            
            t1 = time()
            print(f"HPO data is downloaded in {round((t1-t0) / 60, 2)} mins")            
    
    def process_phenotype_disease(self):
        if not hasattr(self, "hpo_phenotype_disease"):
            self.download_hpo_data()
        if not hasattr(self, "mondo_mappings"):
            self.prepare_mondo_mappings()
            
        print("Started processing HPO phenotype-disease interaction data")
        t0 = time()            
        
        df_list = []
        for hpo_id, diseases in self.hpo_phenotype_disease.items():
            if hpo_id == "hpo_id":
                continue
                
            for disease in diseases:
                 if disease.omim.split(":")[0] == "OMIM" and self.mondo_mappings.get(disease.omim.split(":")[1]):
                        if disease.pmid:
                            if ";" in disease.pmid:
                                pmid = "|".join([i.replace("PMID:", "") for i in disease.pmid.split(";")])
                            else:
                                pmid = disease.pmid.replace("PMID:", "")
                        else:
                            pmid = None
                            
                        df_list.append((hpo_id, self.mondo_mappings.get(disease.omim.split(":")[1]),
                                       pmid, disease.evidence))
        
        
        df = pd.DataFrame(df_list, columns=["hpo_id", "disease_id", "pubmed_ids", "evidence"])
        df.fillna(value=np.nan, inplace=True)
            
        df = df.groupby(["hpo_id", "disease_id"], sort=False, as_index=False).aggregate({"hpo_id":"first",
                                                                                          "disease_id":"first",
                                                                                          "pubmed_ids":self.merge_source_column,
                                                                                          "evidence":"first"})
        df.replace("", np.nan, inplace=True)
            
        t1 = time()
        print(f"HPO phenotype-disease interaction data is processed in {round((t1-t0) / 60, 2)} mins")
            
        return df            
        
    def get_nodes(self, label="phenotype"):
        if not hasattr(self, "hpo_terms"):
            self.download_hpo_data()
        
        node_list = []
        
        for term, name in tqdm(self.hpo_terms.items()):
            hpo_id = self.add_prefix_to_id(prefix="hp", identifier=term)
            
            props = {}
            props["name"] = name
            if self.hpo_ontology["synonyms"].get(term):
                if len(self.hpo_ontology["synonyms"].get(term)) == 1:
                    props["synonyms"] = list(self.hpo_ontology["synonyms"].get(term))[0]
                else:                    
                    props["synonyms"] = list(self.hpo_ontology["synonyms"].get(term))
                    
            node_list.append((hpo_id, label, props))
            
        return node_list
    
    def get_edges(self):
        
        print("Preparing all edge types")
        
        edge_list = []
        
        edge_list.extend(self.get_protein_phenotype_edges())
        
        edge_list.extend(self.get_phenotype_hierarchical_edges())
        
        edge_list.extend(self.get_phenotype_disease_edges())
        
        return edge_list
    
    def get_protein_phenotype_edges(self, label="protein_is_associated_with_phenotype"):
        if not hasattr(self, "protein_hpo_annotations"):
            self.download_hpo_data()
            
        print("Preparing protein-phenotype edges")
        
        edge_list = set()
        
        for uniprot_id, annotations in tqdm(self.protein_hpo_annotations.items()):
            protein_id = self.add_prefix_to_id(prefix="uniprot", identifier=uniprot_id)
            for annot in annotations:
                hpo_id = self.add_prefix_to_id(prefix="hp", identifier=annot.hpo_id)
                edge_list.add((None, protein_id, hpo_id, label))
            
        return [i + tuple([{}]) for i in edge_list]
    
    def get_phenotype_hierarchical_edges(self, label="phenotype_is_a_phenotype"):
        if not hasattr(self, "hpo_ontology"):
            self.download_hpo_data()
        
        print("Preparing phenotype hierarchical edges")
        
        edge_list = []
        for child, parents in tqdm(self.hpo_ontology["parents"].items()):
            child_id = self.add_prefix_to_id(prefix="hp", identifier=child)
            for parent in parents:
                parent_id = self.add_prefix_to_id(prefix="hp", identifier=parent)
                edge_list.append((None, child_id, parent_id, label, {}))
                
        return edge_list
    
    def get_phenotype_disease_edges(self, label="phenotype_is_associated_with_disease"):
        phenotype_disease_df = self.process_phenotype_disease()
        
        print("Preparing phenotype-disease edges")
        
        edge_list = []
        
        for index, row in tqdm(phenotype_disease_df.iterrows(), total=phenotype_disease_df.shape[0]):
            _dict = row.to_dict()
            
            hpo_id = self.add_prefix_to_id(prefix="hp", identifier=_dict["hpo_id"])
            disease_id = self.add_prefix_to_id(prefix="MONDO", identifier=_dict["disease_id"])
            
            del _dict["disease_id"], _dict["hpo_id"]
            
            props = {}
            for k, v in _dict.items():
                if str(v) != "nan":
                    if isinstance(v, str) and "|" in v:
                        props[k] = v.split("|")
                    else:
                        props[k] = v
                        
            edge_list.append((None, hpo_id, disease_id, label, props))
            
        return edge_list
    
    def prepare_mondo_mappings(self):
        mondo = ontology.ontology(ontology="mondo", fields=["is_obsolete", "obo_xref"])

        self.mondo_mappings = {}

        mapping_db_list = ["OMIM"]

        for term in mondo:
            if not term.is_obsolete and term.obo_id and "MONDO" in term.obo_id and term.obo_xref:
                for xref in term.obo_xref:
                    if xref.get("database") in mapping_db_list:
                        self.mondo_mappings[xref["id"]] = term.obo_id
            
    def add_prefix_to_id(self, prefix=None, identifier=None, sep=":") -> str:
        """
        Adds prefix to database id
        """
        if self.add_prefix and identifier:
            return normalize_curie(prefix + sep + str(identifier))
        
        return identifier
    
    def merge_source_column(self, element, joiner="|"):        
        _list = []
        for e in list(element.dropna().values):
            if joiner in e:
                for i in e.split(joiner):
                    _list.append(i)
            else:
                _list.append(e)

        return joiner.join(list(dict.fromkeys(_list).keys()))
