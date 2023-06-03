from __future__ import annotations
from pypath.share import curl, settings

from pypath.inputs import collectri, dorothea, trrust, uniprot
import pypath.utils.mapping as mapping

from contextlib import ExitStack

from bioregistry import normalize_curie
from tqdm.notebook import tqdm
from time import time

import pandas as pd
import numpy as np

class TFGen:
    def __init__(self, organism=9606, add_prefix = True):
        self.organism = organism
        self.add_prefix = add_prefix
        
        
        self.effect_mapping = {0:"Unknown", 1:"Activation", -1:"Repression"}
    
    def download_tfgen_data(
        self,
        cache=False,
        debug=False,
        retries=3,
    ):
        """
        Wrapper function to download tf-gen relation data from various databases using pypath.
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
                
    
    def download_dorothea_data(self):
        
        print("Started downloading DoRothEA data")
        t0 = time()
        
        # ONLY FOR HUMAN. LATER ADD MOUSE AS WELL
        self.dorothea_interactions = list(dorothea.dorothea_interactions(organism = self.organism, levels={"A", "B", "C"}))
        
        t1 = time()
        print(f"DoRothEA data is downloaded in {round((t1-t0) / 60, 2)} mins")
        
    def download_collectri_data(self):
        print("Started downloading CollecTRI data")
        t0 = time()
        
        self.collectri_interactions = list(collectri.collectri_interactions())
        
        self.uniprot_to_entrez = {k:v.strip(";").split(";")[0] for k,v in uniprot.uniprot_data("database(GeneID)", self.organism, True).items()}
        
        t1 = time()
        print(f"CollecTRI data is downloaded in {round((t1-t0) / 60, 2)} mins")
        
    def download_trrust_data(self):
        print("Started downloading TRRUST data")
        t0 = time()
        
        # ONLY FOR HUMAN. LATER ADD MOUSE AS WELL
        self.trrust_gene_symbol_to_entrez_id = {entry["gene_symbol"]:entry["entrez_id"] for entry in trrust.scrape_human()}
        
        # ONLY FOR HUMAN. LATER ADD MOUSE AS WELL
        self.trrust_interactions = trrust.trrust_human()        
        
        t1 = time()
        print(f"CollecTRI data is downloaded in {round((t1-t0) / 60, 2)} mins")
        
    def process_dorothea_tf_gene(self):
        
        if not hasattr(self, "dorothea_interactions"):
            self.download_dorothea_data()
            
        print("Started processing DoRothEA tf-gen data")
        t0 = time()
        
        df_list = []
        for interaction in self.dorothea_interactions:
            tf = self.map_gene_symbol_to_entrez_id(interaction.tf)
            target = self.map_gene_symbol_to_entrez_id(interaction.target)
            if tf and target:
                tf = list(tf)[0]
                target = list(target)[0]
                effect = self.effect_mapping[interaction.effect]
                
                pubmed_id = None
                if interaction.pubmed:
                    pubmed_id = interaction.pubmed
                    
                df_list.append((tf, target, pubmed_id, effect, interaction.level))
                
        df = pd.DataFrame(df_list, columns=["tf", "target", "pubmed_id", "tf_effect", "dorothea_confidence_level"])
        
        df.drop_duplicates(subset=["tf", "target"], ignore_index=True, inplace=True)
        
        df.fillna(value=np.nan, inplace=True)
        
        df["source"] = "DoRothEA"
        
        t1 = time()
        print(f"DoRothEA tf-gen data is processed in {round((t1-t0) / 60, 2)} mins")
        
        return df
    
    def process_collectri_tf_gene(self):
        if not hasattr(self, "collectri_interactions"):
            self.download_collectri_data()
            
        print("Started processing CollecTRI tf-gen data")
        t0 = time()
        
        df_list = []
        for interaction in self.collectri_interactions:
            if str(interaction.tf).startswith("COMPLEX"):
                for tf in list(interaction.tf):
                    if self.uniprot_to_entrez.get(interaction.tf) and self.uniprot_to_entrez.get(interaction.target):
                        tf_entrez = self.uniprot_to_entrez[interaction.tf]
                        target_entrez = self.uniprot_to_entrez[interaction.target]
                        effect = self.effect_mapping[interaction.effect]
                        
                        pubmed_id = None
                        if interaction.pubmed:
                            pubmed_id = interaction.pubmed.replace(";", "|")
                            
                        df_list.append((tf_entrez, target_entrez, pubmed_id, effect))
            else:
                if self.uniprot_to_entrez.get(interaction.tf) and self.uniprot_to_entrez.get(interaction.target):
                    tf_entrez = self.uniprot_to_entrez[interaction.tf]
                    target_entrez = self.uniprot_to_entrez[interaction.target]
                    effect = self.effect_mapping[interaction.effect]
                    
                    pubmed_id = None
                    if interaction.pubmed:                        
                        pubmed_id = interaction.pubmed.replace(";", "|")
                        
                    df_list.append((tf_entrez, target_entrez, pubmed_id, effect))
                    
                    
        df = pd.DataFrame(df_list, columns=["tf", "target", "pubmed_id", "tf_effect"])
        
        df.fillna(value=np.nan, inplace=True)
        
        df = df.groupby(["tf", "target"], sort=False, as_index=False).aggregate({"tf":"first",
                                                                    "target":"first",
                                                                    "pubmed_id":self.merge_source_column,
                                                                    "tf_effect":self.find_conflicting_tf_effects})
        
        
        df.dropna(subset="tf_effect", inplace=True)
                               
        df["source"] = "CollecTRI"
        
        t1 = time()
        print(f"CollecTRI tf-gen data is processed in {round((t1-t0) / 60, 2)} mins")
        
        return df
    
    def process_trrust_tf_gene(self):
        if not hasattr(self, "trrust_interactions"):
            self.download_trrust_data()
            
        print("Started processing TRRUST tf-gen data")
        t0 = time()
        
        df_list = []
        for tf, genes in self.trrust_interactions.items():
            
            if self.trrust_gene_symbol_to_entrez_id.get(tf):
                tf_entrez = self.trrust_gene_symbol_to_entrez_id[tf]
                
                for target in genes:
                    if self.trrust_gene_symbol_to_entrez_id.get(target["gene_symbol"]):
                        target_entrez = self.trrust_gene_symbol_to_entrez_id[target["gene_symbol"]]
                        df_list.append((tf_entrez, target_entrez, target["effect"],))
                        
        df = pd.DataFrame(df_list, columns=["tf", "target", "tf_effect"])
        
        df = df.groupby(["tf", "target"], sort=False, as_index=False).aggregate({"tf":"first",
                                                                    "target":"first",
                                                                    "tf_effect":self.find_conflicting_tf_effects})
        df.dropna(subset="tf_effect", inplace=True)
                               
        df["source"] = "TRRUST"
        
        t1 = time()
        print(f"TRRUST tf-gen data is processed in {round((t1-t0) / 60, 2)} mins")
        
        return df
    
    def merge_tf_gen_data(self):
        trrust_df = self.process_trrust_tf_gene()
        
        dorothea_df = self.process_dorothea_tf_gene()
        
        collectri_df = self.process_collectri_tf_gene()
        
        print("Started merging tf-gen edge data")
        t0 = time()
        
        # merge dorothea and collectri
        merged_df = pd.merge(dorothea_df, collectri_df, how="outer", on=["tf", "target"])
        merged_df.replace("", np.nan, inplace=True)
        
        # merge source column
        merged_df["source"] = merged_df[["source_x", "source_y"]].apply(self.merge_source_column, axis=1)
        merged_df.drop(columns=["source_x", "source_y"], inplace=True)
        
        # merge pubmed_id column
        merged_df["pubmed_id"] = merged_df[["pubmed_id_x", "pubmed_id_y"]].apply(self.merge_source_column, axis=1)
        merged_df.drop(columns=["pubmed_id_x", "pubmed_id_y"], inplace=True)
        
        # merge tf_effect column
        merged_df["tf_effect"] = merged_df[["tf_effect_x", "tf_effect_y"]].apply(self.find_conflicting_tf_effects, axis=1)
        merged_df.drop(columns=["tf_effect_x", "tf_effect_y"], inplace=True)
        merged_df.dropna(subset="tf_effect", inplace=True)
        
        # merge dorothea+collectri and trrust
        merged_df = merged_df.merge(trrust_df, how="outer", on=["tf", "target"])
        merged_df.replace("", np.nan, inplace=True)
        
        # merge source column
        merged_df["source"] = merged_df[["source_x", "source_y"]].apply(self.merge_source_column, axis=1)
        merged_df.drop(columns=["source_x", "source_y"], inplace=True)
        
        # merge tf_effect column
        merged_df["tf_effect"] = merged_df[["tf_effect_x", "tf_effect_y"]].apply(self.find_conflicting_tf_effects, axis=1)
        merged_df.drop(columns=["tf_effect_x", "tf_effect_y"], inplace=True)
        merged_df.dropna(subset="tf_effect", inplace=True)
        
        t1 = time()
        print(f"Tf-gene edge data is merged in {round((t1-t0) / 60, 2)} mins")
        
        return merged_df
    
    def get_edges(self, label="gene_regulates_gene"):
        
        tf_gen_edges_df = self.merge_tf_gen_data()
        
        print("Started writing tf-gen edges")
        
        edge_list = []
        for index, row in tqdm(tf_gen_edges_df.iterrows(), total=tf_gen_edges_df.shape[0]):
            _dict = row.to_dict()
            
            tf_id = self.add_prefix_to_id(prefix="ncbigene", identifier=_dict["tf"])
            target_id = self.add_prefix_to_id(prefix="ncbigene", identifier=_dict["target"])
            
            del _dict["tf"], _dict["target"]
            
            props = {}
            for k, v in _dict.items():
                if str(v) != "nan":
                    if isinstance(v, str) and "|" in v:
                        props[k] = v.split("|")
                    else:
                        props[k] = v
                        
            edge_list.append((None, tf_id, target_id, label, props))
            
        return edge_list            
            
    def add_prefix_to_id(self, prefix=None, identifier=None, sep=":") -> str:
        """
        Adds prefix to database id
        """
        if self.add_prefix and identifier:
            return normalize_curie(prefix + sep + str(identifier))
        
        return identifier
        
        
    def map_gene_symbol_to_entrez_id(self, gene_symbol):
        return mapping.map_name(gene_symbol, "genesymbol", "entrez")
    
    def merge_source_column(self, element, joiner="|"):
        
        _list = []
        for e in list(element.dropna().values):
            if joiner in e:
                for i in e.split(joiner):
                    _list.append(i)
            else:
                _list.append(e)
                
        return joiner.join(list(dict.fromkeys(_list).keys()))
    
    def find_conflicting_tf_effects(self, element):
        tf_effects = set()
        for e in list(element.dropna().values):
            tf_effects.add(e)

        if len(tf_effects) > 2:
            return np.nan
        elif len(tf_effects) == 1:
            return list(tf_effects)[0]
        elif tf_effects == {'Unknown', 'Activation'} or tf_effects == {'Unknown', 'Repression'}:
            _list = list(tf_effects)
            _list.remove("Unknown")
            return _list[0]
        else:
            np.nan
