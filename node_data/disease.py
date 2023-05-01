from __future__ import annotations

from pypath.share import curl, settings
from pypath.inputs import (
    pathophenodb,
    ctdbase,
    clinvar,
    disgenet,
    pharos,
    chembl,
    diseases,
    opentargets,
    drugbank,
    uniprot,
    unichem,
)
import kegg_local

from pypath.inputs import ontology
from pypath.formats import obo
from pypath.utils import mapping

from typing import Union
from contextlib import ExitStack

from bioregistry import normalize_curie
from tqdm.notebook import tqdm
from time import time
import collections

from enum import Enum, auto

import numpy as np
import pandas as pd

class DiseaseNodeField(Enum):
    NAME = "name"
    SYNONYMS = "obo_synonym"
    
    # xrefs
    UMLS = "UMLS"
    DOID = "DOID"
    MESH = "MESH"
    OMIM = "OMIM"
    EFO = "EFO"
    ORPHANET = "Orphanet"
    HP = "HP" # Human Phenotype Ontology
    ICD10CM = "ICD10CM"
    NCIT = "NCIT"
    ICD9 = "ICD9"
    MEDDRA = "MedDRA"
    
    @classmethod
    def get_xrefs(cls):
        return [cls.UMLS.value, cls.DOID.value, cls.MESH.value, cls.OMIM.value, cls.EFO.value,
               cls.ORPHANET.value, cls.HP.value, cls.ICD10CM.value, cls.NCIT.value, cls.ICD9.value,
               cls.MEDDRA.value]

class DiseaseNodeType(Enum):
    pass

class DiseaseEdgeType(Enum):
    MONDO_HIERARCHICAL_RELATIONS = auto()
    ORGANISM_TO_DISEASE = auto()
    GENE_TO_DISEASE = auto()
    DISEASE_TO_DRUG = auto()
    

class GENE_TO_DISEASE_INTERACTION_FIELD(Enum):
    PUBMED_ID = ""
    

class Disease:
    """
    Class that downloads disease data using pypath and reformats it to be ready
    for import into a BioCypher database.
    """
    
    def __init__(self, drugbank_user, drugbank_passwd, node_types: Union[list[DiseaseNodeType], None] = None, 
                 disease_node_fields: Union[list[DiseaseNodeField], None] = None,
                edge_types: Union[list[DiseaseEdgeType], None] = None,
                 add_prefix = True,):
        
        self.drugbank_user = drugbank_user
        self.drugbank_passwd = drugbank_passwd
        self.add_prefix = add_prefix
        
        
        self.set_node_and_edge_types(edge_types=edge_types)
        self.set_node_and_edge_fields(disease_node_fields=disease_node_fields)
    
    def download_disease_data(
        self,
        cache=False,
        debug=False,
        retries=3,
    ):
        """
        Wrapper function to download disease data from various databases using pypath.
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
                
            
            self.download_mondo_data()
            
            self.prepare_mappings()
            
            self.download_pathophenodb_data()
            
            self.download_ctd_data()
            
            
    def download_mondo_data(self):
        fields = ["is_obsolete"]
        if DiseaseNodeField.SYNONYMS.value in self.disease_node_fields:
            fields.append(DiseaseNodeField.SYNONYMS.value)
            
        if set(DiseaseNodeField.get_xrefs()).intersection(set(self.disease_node_fields)):
            fields.append("obo_xref")
         
        t0 = time()
        
        self.mondo = ontology.ontology(ontology="mondo", fields=fields)
        
        t1 = time()
        print(f"Mondo data is downloaded in {round((t1-t0) / 60, 2)} mins")
        
        if DiseaseEdgeType.MONDO_HIERARCHICAL_RELATIONS in self.edge_types:
            t0 = time()
            
            mondo_obo_reader = obo.Obo("http://purl.obolibrary.org/obo/mondo.obo")

            mondo_obo_reader.parent_terms()

            self.mondo_hierarchical_relations = {k:v for k,v in mondo_obo_reader.parents.items() if v}
            
            t1 = time()
            print(f"Mondo hierarchical relations data is downloaded in {round((t1-t0) / 60, 2)} mins")
            
            
    def download_pathophenodb_data(self):
        if DiseaseEdgeType.ORGANISM_TO_DISEASE in self.edge_types:
            t0 = time()

            self.pathopheno_organism_disease_int = pathophenodb.disease_pathogen_interactions()

            t1 = time()
            print(f"PathophenoDB organism-disease interaction data is downloaded in {round((t1-t0) / 60, 2)} mins")
            
    def download_ctd_data(self):
        if DiseaseEdgeType.GENE_TO_DISEASE in self.edge_types:
            t0 = time()
            
            self.ctdbase_gda = ctdbase.ctdbase_relations(relation_type='gene_disease')
            
            t1 = time()
            print(f"CTD gene-disease interaction data is downloaded in {round((t1-t0) / 60, 2)} mins")
            
            
        if DiseaseEdgeType.DISEASE_TO_DRUG in self.edge_types:
            t0 = time()
            
            self.ctdbase_cd = ctdbase.ctdbase_relations(relation_type='chemical_disease')
            
            self.drugbank_data = drugbank.DrugbankFull(user = self.drugbank_user, passwd = self.drugbank_passwd)            
            drugbank_drugs_detailed = self.drugbank_data.drugbank_drugs_full(fields = ["cas_number"])
            self.cas_to_drugbank = {drug.cas_number:drug.drugbank_id for drug in drugbank_drugs_detailed if drug.cas_number}
            
            t1 = time()
            print(f"CTD chemical-disease interaction data is downloaded in {round((t1-t0) / 60, 2)} mins")
            
    def download_chembl_data(self):
        if DiseaseEdgeType.DISEASE_TO_DRUG in self.edge_types:
            t0 = time()
            
            self.chembl_disease_drug = chembl.chembl_drug_indications()
            
            self.chembl_to_drugbank = {k:list(v)[0] for k,v in unichem.unichem_mapping("chembl", "drugbank").items()}
            
            t1 = time()
            print(f"CHEMBL drug indication data is downloaded in {round((t1-t0) / 60, 2)} mins")
            
            
    def download_kegg_data(self):
        if DiseaseEdgeType.DISEASE_TO_DRUG in self.edge_types:
            t0 = time()
            
            self.kegg_drug_disease = kegg_local.drug_to_disease()
            
            if not hasattr(self, "drugbank_data"):
                self.drugbank_data = drugbank.DrugbankFull(user = self.drugbank_user, passwd = self.drugbank_passwd)
            
            drugbank_drugs_external_ids = self.drugbank_data.drugbank_external_ids_full()
            self.kegg_drug_to_drugbank = {v.get("KEGG Drug"):k for k, v in drugbank_drugs_external_ids.items() if v.get("KEGG Drug")}
            
            kegg_disease_ids = kegg_local._Disease()._data.keys()
            
            self.kegg_diseases_mappings = {}
            for dis in kegg_disease_ids:
                result = kegg_local.get_diseases(dis)
                self.kegg_diseases_mappings[dis] = result[0].db_links
                
            t1 = time()
            print(f"KEGG drug indication data is downloaded in {round((t1-t0) / 60, 2)} mins")
            
    def download_diseases_data(self):
        if DiseaseEdgeType.GENE_TO_DISEASE in self.edge_types:
            t0 = time()
            
            self.diseases_knowledge = diseases.knowledge_filtered()
            self.diseases_experimental = diseases.experiments_filtered()

            
            if not hasattr(self, "ensembl_gene_to_uniprot"):
                uniprot_to_ensembl = uniprot.uniprot_data("database(Ensembl)", "9606", True)

                self.ensembl_gene_to_uniprot = {self.ensembl_transcript_to_ensembl_gene(ensts):uniprot_id for uniprot_id, ensts in uniprot_to_ensembl.items() if self.ensembl_transcript_to_ensembl_gene(ensts)}

            
            self.ensembl_protein_to_uniprot = {}
            for k,v in self.ensembl_gene_to_uniprot.items():
                ensps = self.ensembl_gene_to_ensembl_protein(k)
                for p in ensps:
                    self.ensembl_protein_to_uniprot[p] = v
            
            t1 = time()
            print(f"DISEASES gene-disease interaction data is downloaded in {round((t1-t0) / 60, 2)} mins")        
        
    def prepare_mappings(self):
        self.mondo_mappings = collections.defaultdict(dict)
        mapping_db_list = ["UMLS", "DOID", "MESH", "OMIM", "EFO", "Orphanet", "HP", "ICD10CM"]
        
        if not hasattr(self, "mondo"):
            self.download_mondo_data()
            
        for term in self.mondo:
            if not term.is_obsolete and term.obo_id and "MONDO" in term.obo_id and term.obo_xref:
                for xref in term.obo_xref:
                    if xref.get("database") in mapping_db_list:
                        db = mapping_db_list[mapping_db_list.index(xref.get("database"))]
                        self.mondo_mappings[db][xref["id"]] = term.obo_id
                        
    def process_ctd_chemical_disease(self): # DUPLICATELERE BAK     
        
        if not hasattr(self, "mondo_mappings"):
            self.prepare_mappings()
            
        if not hasattr(self, "ctdbase_cd"):
            self.download_ctd_data()
        
        print("Processing CTD chemical-disease data")
        t0 = time()
        
        df_list = []
        for interaction in tqdm(self.ctdbase_cd):
            if interaction.CasRN and interaction.DirectEvidence and interaction.DirectEvidence == "therapeutic"\
            and interaction.PubMedIDs and self.cas_to_drugbank.get(interaction.CasRN):                
                db = interaction.DiseaseID.split(":")[0]
                disease_id = interaction.DiseaseID.split(":")[1]
                if self.mondo_mappings[db].get(disease_id, None):
                    disease_id = self.mondo_mappings[db].get(disease_id)
                    drug_id = self.cas_to_drugbank.get(interaction.CasRN)                    
                    if isinstance(interaction.PubMedIDs, list):
                        pubmed_ids = "|".join(interaction.PubMedIDs)
                    else:
                        pubmed_ids = interaction.PubMedIDs
                        
                    df_list.append((disease_id, drug_id, pubmed_ids))
        
        df = pd.DataFrame(df_list, columns=["disease_id", "drug_id", "pubmed_ids"])
        df["source"] = "CTD"
        
        t1 = time()
        print(f"CTD chemical-disease interaction data is processed in {round((t1-t0) / 60, 2)} mins")
        
        return df
    
    def process_chembl_drug_indication(self):        
        
        if not hasattr(self, "mondo_mappings"):
            self.prepare_mappings()
        if not hasattr(self, "chembl_disease_drug"):
            self.download_chembl_data()
        
        print("Processing CHEMBL drug indication data")
        t0 = time()
        
        df_list = []
        for dd in tqdm(self.chembl_disease_drug):
            if dd.efo_id and self.chembl_to_drugbank.get(dd.molecule_chembl) and dd.efo_id.split(":")[0] in list(self.mondo_mappings.keys())+["MONDO"]\
            and dd.max_phase > 0.0:
                db = dd.efo_id.split(":")[0]
                disease_id = dd.efo_id.split(":")[1]
                drug_id = self.chembl_to_drugbank.get(dd.molecule_chembl)
                
                if db == "MONDO":
                    df_list.append((disease_id, drug_id, dd.max_phase))
                else:
                    if self.mondo_mappings[db].get(disease_id):
                        disease_id = self.mondo_mappings[db].get(disease_id)
                        df_list.append((disease_id, drug_id, dd.max_phase))
        
        df = pd.DataFrame(df_list, columns=["disease_id", "drug_id", "max_phase"])
        df["source"] = "ChEMBL"
        
        df.drop_duplicates(subset=["disease_id", "drug_id"], ignore_index=True, inplace=True)
        
        t1 = time()
        print(f"CHEMBL drug indication data is processed in {round((t1-t0) / 60, 2)} mins")
        
        return df
    
    def process_kegg_drug_indication(self):
        if not hasattr(self, "mondo_mappings"):
            self.prepare_mappings()
        if not hasattr(self, "kegg_drug_disease"):
            self.download_kegg_data()
        
        print("Processing CHEMBL drug indication data")
        t0 = time()
        
        kegg_dbs_to_mondo_dbs = {"MeSH":"MESH", "OMIM":"OMIM", "ICD-10":"ICD10CM",}
        
        df_list = []
        for drug, kegg_diseases in self.kegg_drug_disease.items():
            if self.kegg_drug_to_drugbank.get(drug):
                for interaction in kegg_diseases.DiseaseEntries:
                    disease_id = None
                    if self.kegg_diseases_mappings.get(interaction.disease_id):
                        found = False
                        
                        for db in kegg_dbs_to_mondo_dbs.keys():
                            if found:
                                break
                            
                            if self.kegg_diseases_mappings[interaction.disease_id].get(db):
                                for ref in self.ensure_iterable(self.kegg_diseases_mappings[interaction.disease_id][db]):
                                    if found:
                                        break
                                        
                                    if self.mondo_mappings[kegg_dbs_to_mondo_dbs[db]].get(ref):
                                        disease_id = self.mondo_mappings[kegg_dbs_to_mondo_dbs[db]].get(ref)
                                        found = True
                    
                    if disease_id:
                        drug_id = self.kegg_drug_to_drugbank.get(drug)
                        df_list.append((disease_id, drug_id))
                        
                        
        df = pd.DataFrame(df_list, columns=["disease_id", "drug_id",])
        df["source"] = "KEGG"
        
        df.drop_duplicates(subset=["disease_id", "drug_id"], ignore_index=True, inplace=True)
        
        return df
    
    def process_diseases_gene_disease(self):

        if not hasattr(self, "mondo_mappings"):
            self.prepare_mappings()
        if not hasattr(self, "diseases_knowledge") or not hasattr(self, "diseases_experimental"):
            self.download_diseases_data()
        
        print("Processing DISEASES gene-disease data")
        t0 = time()
        
        df_list = []
        for dg in tqdm(self.diseases_knowledge):
            if self.ensembl_protein_to_uniprot.get(dg.gene_identifier) and self.mondo_mappings["DOID"].get(dg.disease_identifier.split(":")[1]):
                
                gene_id = self.ensembl_protein_to_uniprot.get(dg.gene_identifier)
                disease_id = self.mondo_mappings["DOID"].get(dg.disease_identifier.split(":")[1])
                df_list.append((gene_id, disease_id))

        diseases_knowledge_df = pd.DataFrame(df_list, columns=["gene_id", "disease_id",])
        diseases_knowledge_df["source"] = "DISEASES Knowledge"

        df_list = []
        for dg in tqdm(self.diseases_experimental):
            if self.ensembl_protein_to_uniprot.get(dg.gene_identifier) and self.mondo_mappings["DOID"].get(dg.disease_identifier.split(":")[1]):
                gene_id = self.ensembl_protein_to_uniprot.get(dg.gene_identifier)
                disease_id = self.mondo_mappings["DOID"].get(dg.disease_identifier.split(":")[1])
                df_list.append((gene_id, disease_id))

        diseases_experimental_df = pd.DataFrame(df_list, columns=["gene_id", "disease_id",])
        diseases_experimental_df["source"] = "DISEASES Experimental"
        
        t1 = time()
        print(f"DISEASES gene-disease data is processed in {round((t1-t0) / 60, 2)} mins")

        return diseases_knowledge_df, diseases_experimental_df 
        
    def merge_disease_drug_edge_data(self):
        # Prepare dataframes for merging
        ctd_df = self.process_ctd_chemical_disease()
        
        chembl_df = self.process_chembl_drug_indication()
        
        kegg_df = self.process_kegg_drug_indication()
        
        print("Started merging disease-drug data")
        t0 = time()
        
        # MERGE CHEMBL AND CTD
        merged_df = chembl_df.merge(ctd_df, how="outer", on=["disease_id", "drug_id"])
        
        merged_df["source"] = merged_df[["source_x", "source_y"]].apply(self.merge_source_column, axis=1)
        
        merged_df.drop(columns=["source_x", "source_y"], inplace=True)
        
        # MERGE CHEMBL+CTD AND KEGG
        merged_df = merged_df.merge(kegg_df, how="outer", on=["disease_id", "drug_id"])
        
        merged_df["source"] = merged_df[["source_x", "source_y"]].apply(self.merge_source_column, axis=1)
        
        merged_df.drop(columns=["source_x", "source_y"], inplace=True)
        
        t1 = time()
        print(f"Disease-drug edge data is merged in {round((t1-t0) / 60, 2)} mins")
        
        return merged_df
            
    def get_nodes(self, label="disease"):
        print("Preparing Disease nodes")
        
        node_list = []
        
        xref_dbs = list(set(DiseaseNodeField.get_xrefs()).intersection(self.disease_node_fields))
        
        for term in tqdm(self.mondo):
            if not term.is_obsolete and term.obo_id and "MONDO" in term.obo_id:
                disease_id = self.add_prefix_to_id(prefix="MONDO", identifier=term.obo_id)
                props = {}
                
                if DiseaseNodeField.NAME.value in self.disease_node_fields and term.label:
                    props[DiseaseNodeField.NAME.value] = term.label
                    
                if DiseaseNodeField.SYNONYMS.value in self.disease_node_fields and term.obo_synonym:
                    synonym_set = set()
                    for syn in term.obo_synonym:
                        synonym_set.add(syn["name"])

                    props[DiseaseNodeField.SYNONYMS.value] = list(synonym_set)
                    
                if xref_dbs and term.obo_xref:
                    for xref in term.obo_xref:
                        if xref["database"] in xref_dbs:
                            props[xref["database"]] = xref["id"]
                            
                            
                node_list.append((disease_id, label, props))
                
        return node_list
        
    def get_mondo_hiererchical_edges(self, label="disease_is_a_disease") -> list:
        print("Preparing Mondo hiererchical edges")
        
        edge_list = []
        
        for source, target_list in tqdm(self.mondo_hierarchical_relations.items()):
            source_id = self.add_prefix_to_id(prefix="MONDO", identifier=source)
            for target in target_list:
                target_id = self.add_prefix_to_id(prefix="MONDO", identifier=target)
                edge_list.append((None, source_id, target_id, label, {}))

        return edge_list
        
    def get_organism_disease_edges(self, label="causes") -> list:
        print("Preparing organism-disease edges")
        
        if not hasattr(self, "pathopheno_organism_disease_int"):
            self.download_pathophenodb_data()
        
        if not hasattr(self, "mondo_mappings"):
            self.prepare_mappings()
            
        edge_list = []
        
        for interaction in tqdm(self.pathopheno_organism_disease_int):
            if interaction.evidence == "manual assertion" and self.mondo_mappings["DOID"].get(interaction.disease_id.split(":")[1], None):

                disease_id = self.add_prefix_to_id(prefix="MONDO", identifier=self.mondo_mappings["DOID"].get(interaction.disease_id.split(":")[1]))
                organism_id = self.add_prefix_to_id(prefix="ncbitaxon", identifier=interaction.pathogen_taxid)

                edge_list.append((None, organism_id, disease_id, label, {}))
                    
        return edge_list
    
    def get_disease_drug_edges(self, label="disease_is_treated_by_drug"):
        print("Started writing disease-drug edges")
        
        disease_drug_edges_df = self.merge_disease_drug_edge_data()
        
        edge_list = []
        for index, row in tqdm(disease_drug_edges_df.iterrows(), total=disease_drug_edges_df.shape[0]):
            _dict = row.to_dict()
            
            disease_id = self.add_prefix_to_id(prefix="MONDO", identifier=_dict["disease_id"])
            drug_id = self.add_prefix_to_id(prefix="drugbank", identifier=_dict["drug_id"])
            
            del _dict["disease_id"], _dict["drug_id"]
            
            props = {}
            for k, v in _dict.items():
                if str(v) != "nan":
                    if isinstance(v, str) and "|" in v:
                        props[k] = v.split("|")
                    else:
                        props[k] = v
                        
            edge_list.append((None, disease_id, drug_id, label, props))            
        
        return edge_list
        
    def get_gene_disease_edges(self): # GENE-DISEASE'E PROCESS GEREKİYOR
        print("Preparing gene-disease edges")
        
        if not hasattr(self, "mondo_mappings"):
            self.prepare_mappings()
            
        for interaction in tqdm(self.ctdbase_gda):
            if interaction.PubMedIDs:
                db = interaction.DiseaseID.split(":")[0]
                if self.mondo_mappings[db].get(interaction.DiseaseID.split(":")[1], None):
                    disease_id = self.mondo_mappings[db].get(interaction.DiseaseID.split(":")[1])
                    gene_id = interaction.GeneID
                    
                    _dict = interaction._asdict()
        
    
    def add_prefix_to_id(self, prefix=None, identifier=None, sep=":") -> str:
        """
        Adds prefix to uniprot id
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
    
    def ensure_iterable(self, element):
        if isinstance(element, (list, tuple, set)):
            return element
        else:
            return [element]
        
    def ensembl_transcript_to_ensembl_gene(self, enst_ids):
        enst_id_list = [e.split(" ")[0].split(".")[0] for e in enst_ids.split(";") if e]
        if enst_id_list:
            ensg_id_list = set([list(mapping.map_name(_id, "enst_biomart", "ensg_biomart"))[0] for _id in enst_id_list if mapping.map_name(_id, "enst_biomart", "ensg_biomart")])
            if ensg_id_list:
                return list(ensg_id_list)[0]
            else:
                return None
        else:
            return None
        
    def ensembl_gene_to_ensembl_protein(self, ensg_id):
        return mapping.map_name(ensg_id, "ensg_biomart", "ensp_biomart")
    
    def set_node_and_edge_types(self, edge_types):
        if edge_types:
            self.edge_types = edge_types
        else:
            self.edge_types = [field for field in DiseaseEdgeType]
        
    def set_node_and_edge_fields(self, disease_node_fields):
        if disease_node_fields:
            self.disease_node_fields = [field.value for field in disease_node_fields]
        else:
            self.disease_node_fields = [field.value for field in DiseaseNodeField]
