from __future__ import annotations
from pypath.share import curl, settings

from pypath.inputs import sider, drugbank, offsides, adrecs

from contextlib import ExitStack

from bioregistry import normalize_curie
from tqdm.notebook import tqdm
from time import time

import pandas as pd
import numpy as np


class SideEffect:
    def __init__(self, drugbank_user, drugbank_passwd, add_prefix=True):
        self.drugbank_user = drugbank_user
        self.drugbank_passwd = drugbank_passwd
        
        self.add_prefix = add_prefix
        
        
    def download_side_effect_data(
        self,
        cache=False,
        debug=False,
        retries=3,
    ):
        """
        Wrapper function to download side effect data from various databases using pypath.
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
                
    
    def dowload_sider_data(self):
        
        print("Started downloading Sider data")
        t0 = time()
        
        self.drugbank_data = drugbank.DrugbankFull(user = self.drugbank_user, passwd = self.drugbank_passwd)
        drugbank_drug_names = self.drugbank_data.drugbank_drugs_full(fields = ["name"])
        self.drugbank_name_to_drugbank_id_dict = {i.name.lower():i.drugbank_id for i in drugbank_drug_names}
        
        sider_drug = sider.sider_drug()
        self.cid_to_sider_drug_name = {i.cid:i.drug_name for i in sider_drug}
        
        sider_meddra_tsv = sider.sider_meddra_tsv()
        self.umls_to_meddra_id = {i.cid:{"meddra_id":i.meddra_id, "name":i.side_effect_name} for i in sider_meddra_tsv}
        
        self.sider_meddra_with_freq = sider.sider_meddra_with_freq()
        
        t1 = time()
        print(f"Sider data is downloaded in {round((t1-t0) / 60, 2)} mins")
        
    def download_offsides_data(self):
        print("Started downloading OffSides data")
        t0 = time()
        
        if not hasattr(self, "drugbank_data"):
            self.drugbank_data = drugbank.DrugbankFull(user = self.drugbank_user, passwd = self.drugbank_passwd)
            
        drugbank_external_ids = self.drugbank_data.drugbank_external_ids_full()
        self.rxcui_to_drugbank = {v.get("RxCUI"):k for k, v in drugbank_external_ids.items() if v.get("RxCUI")}
        
        self.offsides_data = offsides.offsides()
        
        t1 = time()
        print(f"OffSides data is downloaded in {round((t1-t0) / 60, 2)} mins")
        
    def download_adrecs_data(self):
        print("Started downloading ADReCS data")
        t0 = time()
        
        adrecs_drug_information = adrecs.adrecs_drug_information()
        self.adrecs_drug_id_to_drugbank_id = {dr.drug_id: dr.drugbank_id for dr in adrecs_drug_information if dr.drugbank_id}
        
        self.adrecs_terminology = adrecs.adrecs_terminology()
        
        self.adrecs_adr_id_to_adrecs_drug_id = {dr.adr_id: dr.drug_id for dr in adrecs.adrecs_drugs()}
        
        t1 = time()
        print(f"ADReCS data is downloaded in {round((t1-t0) / 60, 2)} mins")
        
    def process_sider_drug_side_effect(self):
        if not hasattr(self, "sider_meddra_with_freq"):
            self.dowload_sider_data()
            
        print("Started processing Sider drug-side effect data")
        t0 = time()
        
        df_list = []
        for interaction in self.sider_meddra_with_freq:
            if self.drugbank_name_to_drugbank_id_dict.get(self.cid_to_sider_drug_name.get(interaction.cid)) and self.umls_to_meddra_id.get(interaction.umls_concept_id_on_MedDRA):
                drugbank_id = self.drugbank_name_to_drugbank_id_dict.get(self.cid_to_sider_drug_name.get(interaction.cid))
                meddra_id = self.umls_to_meddra_id[interaction.umls_concept_id_on_MedDRA]["meddra_id"]
                df_list.append((drugbank_id, meddra_id, interaction.frequency))
                
        df = pd.DataFrame(df_list, columns=["drugbank_id", "meddra_id", "frequency"])
        
        df.drop_duplicates(subset=["drugbank_id", "meddra_id"], ignore_index=True, inplace=True)
        
        df["source"] = "Sider"
        
        t1 = time()
        print(f"Sider drug-side effect data is processed in {round((t1-t0) / 60, 2)} mins")
        
        return df
    
    def process_offsides_drug_side_effect(self):
        if not hasattr(self, "offsides_data"):
            self.download_offsides_data()
            
        print("Started processing OffSides drug-side effect data")
        t0 = time()
        
        df_list = []
        for interaction in self.offsides_data:
            if self.rxcui_to_drugbank.get(interaction.drug_rxnorn_id):
                df_list.append((self.rxcui_to_drugbank[interaction.drug_rxnorn_id], interaction.condition_meddra_id,
                               round(float(interaction.PRR), 3),))
                
        df = pd.DataFrame(df_list, columns=["drugbank_id", "meddra_id", "proportional_reporting_ratio"])
        
        df.drop_duplicates(subset=["drugbank_id", "meddra_id"], ignore_index=True, inplace=True)
        
        df["source"] = "OffSides"
        
        t1 = time()
        print(f"OffSides drug-side effect data is processed in {round((t1-t0) / 60, 2)} mins")
        
        return df
    
    def process_adrecs_drug_side_effect(self):
        if not hasattr(self, "adrecs_terminology"):
            self.download_adrecs_data()
            
        print("Started processing ADReCS drug-side effect data")
        t0 = time()
        
        df_list = []
        for interaction in self.adrecs_terminology:
            if self.adrecs_drug_id_to_drugbank_id.get(self.adrecs_adr_id_to_adrecs_drug_id.get(interaction.adr_id)):
                drugbank_id = self.adrecs_drug_id_to_drugbank_id[self.adrecs_adr_id_to_adrecs_drug_id[interaction.adr_id]]
                df_list.append((drugbank_id, str(interaction.meddra_code), ))
                
        df = pd.DataFrame(df_list, columns=["drugbank_id", "meddra_id"])
        
        df.drop_duplicates(subset=["drugbank_id", "meddra_id"], ignore_index=True, inplace=True)
        
        df["source"] = "ADReCS"           
        
        t1 = time()
        print(f"ADReCS drug-side effect data is processed in {round((t1-t0) / 60, 2)} mins")
        
        return df
    
    def merge_drug_side_effect_data(self):
        adrecs_df = self.process_adrecs_drug_side_effect()
        
        sider_df = self.process_sider_drug_side_effect()
        
        offsides_df = self.process_offsides_drug_side_effect()
        
        print("Started merging drug-side effect edge data")
        t0 = time()
        
        merged_df = pd.merge(adrecs_df, sider_df, how="outer", on=["drugbank_id", "meddra_id"])
        
        merged_df["source"] = merged_df[["source_x", "source_y"]].apply(self.merge_source_column, axis=1)
        
        merged_df.drop(columns=["source_x", "source_y"], inplace=True)
        
        merged_df = merged_df.merge(offsides_df, how="outer", on=["drugbank_id", "meddra_id"])
        
        merged_df["source"] = merged_df[["source_x", "source_y"]].apply(self.merge_source_column, axis=1)
        
        merged_df.drop(columns=["source_x", "source_y"], inplace=True)
        
        t1 = time()
        print(f"Drug-side effect edge data is merged in {round((t1-t0) / 60, 2)} mins")
        
        return merged_df
    
    def drug_side_effect_edges(self, label="drug_has_side_effect"):
        drug_side_effect_edges_df = self.merge_drug_side_effect_data()
        
        print("Started writing drug-side effect edges")
        
        edge_list = []
        for index, row in tqdm(drug_side_effect_edges_df.iterrows(), total=drug_side_effect_edges_df.shape[0]):
            _dict = row.to_dict()
            
            meddra_id = self.add_prefix_to_id(prefix="meddra", identifier=_dict["meddra_id"])
            drug_id = self.add_prefix_to_id(prefix="drugbank", identifier=_dict["drugbank_id"])
            
            del _dict["meddra_id"], _dict["drugbank_id"]
            
            props = {}
            for k, v in _dict.items():
                if str(v) != "nan":
                    if isinstance(v, str) and "|" in v:
                        props[k] = v.split("|")
                    else:
                        props[k] = v
            
            edge_list.append((None, drug_id, meddra_id, label, props))
            
        return edge_list
        
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
