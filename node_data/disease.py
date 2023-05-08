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
    humsavar,
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
    DISEASE_TO_DISEASE = auto()
    

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
            
    def download_opentargets_data(self):
        if DiseaseEdgeType.GENE_TO_DISEASE in self.edge_types:
            t0 = time()
            
            self.opentargets_direct = opentargets.overall_direct_score()
                           
            uniprot_to_entrez = uniprot.uniprot_data("database(GeneID)", "9606", True)
            self.uniprot_to_entrez = {k:v.strip(";").split(";")[0] for k,v in uniprot_to_entrez.items()}
            
            if not hasattr(self, "ensembl_gene_to_uniprot"):
                uniprot_to_ensembl = uniprot.uniprot_data("database(Ensembl)", "9606", True)

                self.ensembl_gene_to_uniprot = {self.ensembl_transcript_to_ensembl_gene(ensts):uniprot_id for uniprot_id, ensts in uniprot_to_ensembl.items() if self.ensembl_transcript_to_ensembl_gene(ensts)}
            
            t1 = time()
            print(f"Open Targets direct gene-disease interaction data is downloaded in {round((t1-t0) / 60, 2)} mins")
            
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
            
    def download_clinvar_data(self):
        if DiseaseEdgeType.GENE_TO_DISEASE in self.edge_types:
            t0 = time()
            
            self.clinvar_variant_disease = clinvar.clinvar_raw()
            
            self.clinvar_citation = clinvar.clinvar_citations()
            
            t1 = time()
            print(f"Clinvar variant-disease interaction data is downloaded in {round((t1-t0) / 60, 2)} mins")
            
    def download_humsavar_data(self):
        if DiseaseEdgeType.GENE_TO_DISEASE in self.edge_types:
            t0 = time()
            
            self.humsavar_data = humsavar.humsavar()
            
            if not hasattr(self, "uniprot_to_entrez"):
                uniprot_to_entrez = uniprot.uniprot_data("database(GeneID)", "9606", True)
                self.uniprot_to_entrez = {k:v.strip(";").split(";")[0] for k,v in uniprot_to_entrez.items()}
            
            t1 = time()
            print(f"Humsavar variant-disease interaction data is downloaded in {round((t1-t0) / 60, 2)} mins")
            
    def download_disgenet_data(self, from_csv = True):
        if from_csv:
            print("Skipping downloading part. Will directly process from csv file")
            
        if DiseaseEdgeType.DISEASE_TO_DISEASE in self.edge_types and not from_csv:
            t0 = time()
            
            if not hasattr(self, "disgenet_api"):
                self.disgenet_api = disgenet.DisgenetApi()
                
            if not hasattr(self, "disgenet_disease_ids"):                
                self.disgenet_disease_ids = disgenet.disease_id_mappings().keys()
                
            if not hasattr(self, "disgenet_id_mappings_dict"):
                self.prepare_disgenet_id_mappings()
                
            
            self.disgenet_dda_gene = []
            self.disgenet_dda_variant = []
            for disease_id in tqdm(self.disgenet_disease_ids):
                try:
                    self.disgenet_dda_gene.extend(
                        self.disgenet_api.get_ddas_that_share_genes(disease_id)
                    )
                    self.disgenet_dda_variant.extend(
                        self.disgenet_api.get_ddas_that_share_variants(disease_id)
                    )
                except TypeError:
                    print(f'{disease_id} not available')
                    
            t1 = time()
            print(f"Disgenet disease-disease interaction data is downloaded in {round((t1-t0) / 60, 2)} mins")
            
        if DiseaseEdgeType.GENE_TO_DISEASE in self.edge_types and not from_csv:
            t0 = time()
            
            if not hasattr(self, "disgenet_api"):
                self.disgenet_api = disgenet.DisgenetApi()
                
            if not hasattr(self, "disgenet_disease_ids") or not hasattr(self, "disgenet_id_mappings_dict"):             
                self.disgenet_disease_ids = disgenet.disease_id_mappings().keys()
                
            if not hasattr(self, "disgenet_id_mappings_dict"):
                self.prepare_disgenet_id_mappings()
                
            if not hasattr(self, "uniprot_to_entrez"):
                uniprot_to_entrez = uniprot.uniprot_data("database(GeneID)", "9606", True)
                self.uniprot_to_entrez = {k:v.strip(";").split(";")[0] for k,v in uniprot_to_entrez.items()}
                
            self.gene_symbol_to_uniprot = {}
            for k, v in uniprot.uniprot_data("genes", "9606", True).items():
                for symbol in v.split(" "):
                    self.gene_symbol_to_uniprot[symbol] = k
                    
            
            self.disgenet_gda = []
            self.disgenet_vda = []
            for disease_id in tqdm(self.disgenet_disease_ids):
                try:
                    self.disgenet_gda.extend(
                        self.disgenet_api.get_gdas_by_diseases(disease_id)
                    )
                    self.disgenet_vda.extend(
                        self.disgenet_api.get_vdas_by_diseases(disease_id)
                    )

                except TypeError:
                    print(f'{disease_id} not available')
    
    def prepare_mappings(self):
        self.mondo_mappings = collections.defaultdict(dict)
        mapping_db_list = ["UMLS", "DOID", "MESH", "OMIM", "EFO", "Orphanet", "HP", "ICD10CM", "NCIT"]
        
        
        if not hasattr(self, "mondo"):
            self.download_mondo_data()
            
        for term in self.mondo:
            if not term.is_obsolete and term.obo_id and "MONDO" in term.obo_id and term.obo_xref:
                for xref in term.obo_xref:
                    if xref.get("database") in mapping_db_list:
                        db = mapping_db_list[mapping_db_list.index(xref.get("database"))]
                        self.mondo_mappings[db][xref["id"]] = term.obo_id
                        
                        
    def prepare_disgenet_id_mappings(self):
        
        disgenet_id_mappings = disgenet.disease_id_mappings()
                
        selected_dbs = ["DO", "EFO", "HPO", "MONDO", "MSH", "NCI", "ICD10CM", "OMIM"]
        self.disgenet_id_mappings_dict = collections.defaultdict(dict)

        for disg_id, mappings in disgenet_id_mappings.items():
            map_dict = {}
            for m in mappings.vocabularies:
                if m.vocabulary in selected_dbs and m.vocabulary not in map_dict.keys():
                    map_dict[m.vocabulary] = m.code.split(":")[1] if ":" in m.code else m.code

            self.disgenet_id_mappings_dict[disg_id] = map_dict
                        
    def process_ctd_chemical_disease(self):     
        
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
        
        df = df.groupby(["disease_id", "drug_id"], sort=False, as_index=False).aggregate({"disease_id":"first",
                                                                                          "drug_id":"first",
                                                                                          "pubmed_ids":self.merge_source_column,
                                                                                          "source":"first"})        
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
                    df_list.append((dd.efo_id, drug_id, dd.max_phase))
                else:
                    if self.mondo_mappings[db].get(disease_id):
                        disease_id = self.mondo_mappings[db].get(disease_id)
                        df_list.append((disease_id, drug_id, dd.max_phase))
        
        df = pd.DataFrame(df_list, columns=["disease_id", "drug_id", "max_phase"])
        df["source"] = "ChEMBL"
        
        df.sort_values(by="max_phase", ascending=False, ignore_index=True, inplace=True)
        
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
        
        # Doesnt look necessary
        df.drop_duplicates(subset=["disease_id", "drug_id"], ignore_index=True, inplace=True)
        
        return df
    
    def process_opentargets_gene_disease(self):
        if not hasattr(self, "mondo_mappings"):
            self.prepare_mappings()
        if not hasattr(self, "opentargets_direct"):
            self.download_opentargets_data()
        
        df_list = []
        for disease, targets in self.opentargets_direct.items():
            db = disease.split("_")[0]
            for target in targets:
                if self.uniprot_to_entrez.get(self.ensembl_gene_to_uniprot.get(target["targetId"]))\
                and self.mondo_mappings[db].get(disease.split("_")[1]):
                    gene_id = self.uniprot_to_entrez.get(self.ensembl_gene_to_uniprot.get(target["targetId"]))
                    disease_id = self.mondo_mappings[db].get(disease.split("_")[1])
                    score = round(target["score"], 3) # NE YAPILACAK BU?
                    df_list.append((gene_id, disease_id, score))
                    
        
        df = pd.DataFrame(df_list, columns=["gene_id", "disease_id", "opentargets_score"])
        df["source"] = "Open Targets"
        
        df.sort_values(by="opentargets_score", ascending=False, ignore_index=True, inplace=True)
        
        df.drop_duplicates(subset=["gene_id", "disease_id"], ignore_index=True, inplace=True)
        
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
        
        # drop duplicates
        diseases_knowledge_df.drop_duplicates(subset=["gene_id", "disease_id"], ignore_index=True, inplace=True)

        df_list = []
        for dg in tqdm(self.diseases_experimental):
            if self.ensembl_protein_to_uniprot.get(dg.gene_identifier) and self.mondo_mappings["DOID"].get(dg.disease_identifier.split(":")[1]):
                gene_id = self.ensembl_protein_to_uniprot.get(dg.gene_identifier)
                disease_id = self.mondo_mappings["DOID"].get(dg.disease_identifier.split(":")[1])
                score = float(dg.confidence_score)
                df_list.append((gene_id, disease_id, score))

        diseases_experimental_df = pd.DataFrame(df_list, columns=["gene_id", "disease_id", "diseases_confidence_score"])
        diseases_experimental_df["source"] = "DISEASES Experimental"
        
        df.sort_values(by="diseases_confidence_score", ascending=False, ignore_index=True, inplace=True)
        
        # drop duplicates
        diseases_experimental_df.drop_duplicates(subset=["gene_id", "disease_id"], ignore_index=True, inplace=True)
        
        t1 = time()
        print(f"DISEASES gene-disease data is processed in {round((t1-t0) / 60, 2)} mins")

        return diseases_knowledge_df, diseases_experimental_df
    
    def process_clinvar_gene_disease(self):
        if not hasattr(self, "mondo_mappings"):
            self.prepare_mappings()
        if not hasattr(self, "clinvar_variant_disease"):
            self.download_clinvar_data()
            
        print("Processing Clinvar variant-disease data")
        t0 = time()
        
        selected_clinical_significances = ["Pathogenic", "Likely pathogenic", "Pathogenic/Likely pathogenic"]
        selected_review_status = ["criteria provided, multiple submitters, no conflicts", "reviewed by expert panel",
                         "practice guideline"]
        review_status_dict = {"criteria provided, multiple submitters, no conflicts":2,
                             "reviewed by expert panel":3,
                             "practice guideline":4}
        clinvar_dbs_to_mondo_dbs = {"MONDO":"MONDO", "OMIM":"OMIM", "Orphanet":"Orphanet", "HP":"HP", "MeSH":"MESH"}
        
        df_list = []
        for var in tqdm(self.clinvar_variant_disease):
            if var.entrez and var.clinical_significance in selected_clinical_significances and var.review_status in selected_review_status:
                diseases_set = set()
                for phe in var.phenotype_ids:
                    dbs_and_ids = list(dict.fromkeys(phe.split(":")).keys())
                    
                    if len(dbs_and_ids) > 2:
                        dbs_and_ids = phe.split(":")[1:]
                        
                    if dbs_and_ids[0] in clinvar_dbs_to_mondo_dbs.keys():
                        if dbs_and_ids[0] == "MONDO":
                            diseases_set.add("MONDO:"+dbs_and_ids[1])
                        else:
                            if self.mondo_mappings[clinvar_dbs_to_mondo_dbs[dbs_and_ids[0]]].get(dbs_and_ids[1]):
                                diseases_set.add(self.mondo_mappings[clinvar_dbs_to_mondo_dbs[dbs_and_ids[0]]].get(dbs_and_ids[1]))
                
                if diseases_set:
                    review_status = review_status_dict[var.review_status]
                    for d in diseases_set:
                        df_list.append((var.entrez, d, var.allele, var.clinical_significance, review_status,
                                       "rs"+str(var.rs), var.variation_id)) # NELER EKLENECEK?
                        
        df = pd.DataFrame(df_list, columns=["gene_id", "disease_id", "allele_id", "clinical_significance",
                                           "review_status", "dbsnp_id", "variation_id"])
        df["source"] = "Clinvar"
        df["variant_source"] = "Clinvar"
        
        
        df_list = []
        for cit in self.clinvar_citation:
            if cit.citation_source in ["PubMed", "PubMedCentral"]:
                df_list.append((cit.allele, cit.variation_id, cit.citation_id))


        clinvar_citation_df = pd.DataFrame(df_list, columns=["allele_id", "variation_id", "pubmed_ids"])
        clinvar_citation_df = clinvar_citation_df.groupby(["allele_id", "variation_id"], sort=False, as_index=False).aggregate(
                                                                                         {"allele_id":"first",
                                                                                          "variation_id":"first",
                                                                                          "pubmed_ids":self.merge_source_column,
                                                                                          })
        
        
        df = df.merge(clinvar_citation_df, how="left", on=["allele_id", "variation_id"])
        
        df.sort_values(by="review_status", ascending=False, ignore_index=True, inplace=True)
        
        df.drop_duplicates(subset=["gene_id", "disease_id"], ignore_index=True, inplace=True)
        
        t1 = time()
        print(f"Clinvar variant-disease data is processed in {round((t1-t0) / 60, 2)} mins")
        
        return df
    
    def process_humsavar_gene_disease(self):
        if not hasattr(self, "mondo_mappings"):
            self.prepare_mappings()
        if not hasattr(self, "humsavar_data"):
            self.download_humsavar_data()
            
        print("Processing Humsavar variant-disease data")
        t0 = time()

        df_list = []
        for v in self.humsavar_data:
            if v.variant_category == "LP/P" and v.disease_omim_id and v.dbSNP and self.uniprot_to_entrez.get(v.swiss_prot_ac)\
            and self.mondo_mappings["OMIM"].get(v.disease_omim_id.split(":")[1]):
                gene_id = self.uniprot_to_entrez.get(v.swiss_prot_ac)
                disease_id = self.mondo_mappings["OMIM"].get(v.disease_omim_id.split(":")[1])
                df_list.append((gene_id, disease_id, v.dbSNP))
                
        df = pd.DataFrame(df_list, columns=["gene_id", "disease_id", "dbsnp_id"])
        df["source"] = "Humsavar"
        df["variant_source"] = "Humsavar"
        
        df.drop_duplicates(subset=["gene_id", "disease_id"], ignore_index=True, inplace=True)
                
        t1 = time()
        print(f"Humsavar variant-disease data is processed in {round((t1-t0) / 60, 2)} mins")
        
        return df
    
    def process_disgenet_gene_disease(self, from_csv = True):
        
        if not hasattr(self, "mondo_mappings"):
            self.prepare_mappings()
        if not hasattr(self, "disgenet_id_mappings_dict"):
            self.prepare_disgenet_id_mappings()
        if not hasattr(self, "disgenet_gda") or not hasattr(self, "disgenet_vda"):
            self.download_disgenet_data()
            
        disgenet_dbs_to_mondo_dbs = {"DO":"DOID", "EFO":"EFO", "HPO":"HP", "MSH":"MESH", "NCI":"NCIT", "ICD10CM":"ICD10CM",
                            "OMIM":"OMIM"}
        
        if from_csv:
            print("Processing Disgenet gene-disease data from csv")
            t0 = time()
            
            disgenet_gda_df = pd.read_csv("disgenet_dga.csv")
            
            disgenet_vda_df = pd.read_csv("disgenet_vda.csv")
            
            t1 = time()
            print(f"Disgenet gene-disease data is processed in {round((t1-t0) / 60, 2)} mins")
            
            return something1, something2
        else:
            print("Processing Disgenet gene-disease data")
            t0 = time()
            
            df_list = []
            for gda in self.disgenet_gda:
                diseaseid = self.mondo_mappings["UMLS"].get(dga.diseaseid)
                
                if not diseaseid:
                    if self.disgenet_id_mappings_dict.get(dga.diseaseid):
                        if self.disgenet_id_mappings_dict.get(dga.diseaseid).get("MONDO"):
                            diseaseid = "MONDO:" + self.disgenet_id_mappings_dict.get(dga.diseaseid).get("MONDO")
                        else:
                            map_dict = self.disgenet_id_mappings_dict.get(dga.diseaseid)
                            for db, map_v in map_dict.items():                  
                                if self.mondo_mappings[disgenet_dbs_to_mondo_dbs[db]].get(map_v):                        
                                    diseaseid = self.mondo_mappings[disgenet_dbs_to_mondo_dbs[db]].get(map_v)
                                    break
                
                if diseaseid:
                    df_list.append((gda.geneid, diseaseid, gda.score))
                    
            disgenet_gda_df = pd.DataFrame(df_list, columns=["gene_id", "disease_id", "disgenet_gene_disease_score"])
            disgenet_gda_df["source"] = "Disgenet Gene-Disease"
            
            # DOES NOT LOOK NECESSARY
            disgenet_gda_df.sort_values(by="disgenet_gene_disease_score", ascending=False, ignore_index=True, inplace=True)
            disgenet_gda_df.drop_duplicates(subset=["gene_id", "disease_id"], ignore_index=True, inplace=True)
            
            df_list = []
            for vda in self.disgenet_vda:
                if vda.gene_symbol and self.uniprot_to_entrez.get(self.gene_symbol_to_uniprot.get(vda.gene_symbol)):
                    diseaseid = self.mondo_mappings["UMLS"].get(vda.diseaseid)

                    if not diseaseid:
                        if self.disgenet_id_mappings_dict.get(vda.diseaseid):
                            if self.disgenet_id_mappings_dict.get(vda.diseaseid).get("MONDO"):
                                diseaseid = "MONDO:" + self.disgenet_id_mappings_dict.get(vda.diseaseid).get("MONDO")
                            else:
                                map_dict = self.disgenet_id_mappings_dict.get(vda.diseaseid)
                                for db, map_v in map_dict.items():                  
                                    if self.mondo_mappings[disgenet_dbs_to_mondo_dbs[db]].get(map_v):                        
                                        diseaseid = self.mondo_mappings[disgenet_dbs_to_mondo_dbs[db]].get(map_v)
                                        break

                    if diseaseid:
                        df_list.append((vda.geneid, diseaseid, vda.score, vda.variantid))
                        
            
            disgenet_vda_df = pd.DataFrame(df_list, columns=["gene_id", "disease_id", "disgenet_variant_disease_score", "dbsnp_id"])
            disgenet_vda_df["source"] = "Disgenet Variant-Disease"
            disgenet_vda_df["variant_source"] = "Disgenet Variant-Disease"
            
            disgenet_vda_df.sort_values(by="disgenet_variant_disease_score", ascending=False, ignore_index=True, inplace=True)
            disgenet_vda_df.drop_duplicates(subset=["gene_id", "disease_id"], ignore_index=True, inplace=True)
                
            t1 = time()
            print(f"Disgenet gene-disease data is processed in {round((t1-t0) / 60, 2)} mins")
            
            return disgenet_gda_df, disgenet_vda_df
            
    def process_disgenet_disease_disease(self, from_csv = True):
        
        if not hasattr(self, "mondo_mappings"):
            self.prepare_mappings()
        if not hasattr(self, "disgenet_id_mappings_dict"):
                self.prepare_disgenet_id_mappings()
        if not hasattr(self, "disgenet_dda_gene") or not hasattr(self, "disgenet_dda_variant"):
            self.download_disgenet_data()
        
        disgenet_dbs_to_mondo_dbs = {"DO":"DOID", "EFO":"EFO", "HPO":"HP", "MSH":"MESH", "NCI":"NCIT", "ICD10CM":"ICD10CM",
                            "OMIM":"OMIM"}
        if from_csv:
            print("Processing Disgenet disease-disease data from csv")
            t0 = time()
            
            # DISEASE-DISEASE BY GENE
            disgenet_dda_gene_df = pd.read_csv("disgenet_dda_gene.csv")
            disgenet_dda_gene_df = disgenet_dda_gene_df[["diseaseid1", "diseaseid2", "jaccard_genes"]]
                        
            df_list = []
            for _, dda in disgenet_dda_gene_df.iterrows():
                if round(dda["jaccard_genes"], 3) != 0.0:
                    diseaseid1 = self.mondo_mappings["UMLS"].get(dda["diseaseid1"])
                    diseaseid2 = self.mondo_mappings["UMLS"].get(dda["diseaseid2"])
                    
                    if not diseaseid1:
                        if self.disgenet_id_mappings_dict.get(dda["diseaseid1"]):
                            if self.disgenet_id_mappings_dict.get(dda["diseaseid1"]).get("MONDO"):                    
                                diseaseid1 = "MONDO:" + self.disgenet_id_mappings_dict.get(dda["diseaseid1"]).get("MONDO")
                            else:
                                map_dict = self.disgenet_id_mappings_dict.get(dda["diseaseid1"])
                                for db, map_v in map_dict.items():
                                    if self.mondo_mappings[disgenet_dbs_to_mondo_dbs[db]].get(map_v):
                                        diseaseid1 = self.mondo_mappings[disgenet_dbs_to_mondo_dbs[db]].get(map_v)
                                        break
                                        
                    if not diseaseid2:
                        if self.disgenet_id_mappings_dict.get(dda["diseaseid2"]):
                            if self.disgenet_id_mappings_dict.get(dda["diseaseid2"]).get("MONDO"):
                                diseaseid2 = "MONDO:" + self.disgenet_id_mappings_dict.get(dda["diseaseid2"]).get("MONDO")
                            else:
                                map_dict = self.disgenet_id_mappings_dict.get(dda["diseaseid2"])
                                for db, map_v in map_dict.items():
                                    if self.mondo_mappings[disgenet_dbs_to_mondo_dbs[db]].get(map_v):
                                        diseaseid2 = self.mondo_mappings[disgenet_dbs_to_mondo_dbs[db]].get(map_v)
                                        break
                                            
                    if diseaseid1 and diseaseid2:
                        df_list.append((diseaseid1, diseaseid2, round(dda["jaccard_genes"], 3),))
                
            disgenet_dda_gene_df = pd.DataFrame(df_list, columns=["disease_id1", "disease_id2", "disgenet_jaccard_genes_score"])
            disgenet_dda_gene_df["source"] = "Disgenet Disease-Disease Gene"
            
            disgenet_dda_gene_df.sort_values(by="disgenet_jaccard_genes_score", ascending=False, ignore_index=True, inplace=True)
            disgenet_dda_gene_df.drop_duplicates(subset=["disease_id1", "disease_id2"], ignore_index=True, inplace=True)
            
            # DISEASE-DISEASE BY VARIANT
            disgenet_dda_variant_df = pd.read_csv("disgenet_dda_variant.csv")
            disgenet_dda_variant_df = disgenet_dda_variant_df[["diseaseid1", "diseaseid2", "jaccard_variants"]]
            
            df_list = []
            for _, dda in disgenet_dda_variant_df.iterrows():
                if round(dda["jaccard_variants"], 3) != 0.0:
                    diseaseid1 = self.mondo_mappings["UMLS"].get(dda["diseaseid1"])
                    diseaseid2 = self.mondo_mappings["UMLS"].get(dda["diseaseid2"])
                    
                    if not diseaseid1:
                        if self.disgenet_id_mappings_dict.get(dda["diseaseid1"]):
                            if self.disgenet_id_mappings_dict.get(dda["diseaseid1"]).get("MONDO"):                    
                                diseaseid1 = "MONDO:" + self.disgenet_id_mappings_dict.get(dda["diseaseid1"]).get("MONDO")
                            else:
                                map_dict = self.disgenet_id_mappings_dict.get(dda["diseaseid1"])
                                for db, map_v in map_dict.items():
                                    if self.mondo_mappings[disgenet_dbs_to_mondo_dbs[db]].get(map_v):
                                        diseaseid1 = self.mondo_mappings[disgenet_dbs_to_mondo_dbs[db]].get(map_v)
                                        break
                                        
                    if not diseaseid2:
                        if self.disgenet_id_mappings_dict.get(dda["diseaseid2"]):
                            if self.disgenet_id_mappings_dict.get(dda["diseaseid2"]).get("MONDO"):
                                diseaseid2 = "MONDO:" + self.disgenet_id_mappings_dict.get(dda["diseaseid2"]).get("MONDO")
                            else:
                                map_dict = self.disgenet_id_mappings_dict.get(dda["diseaseid2"])
                                for db, map_v in map_dict.items():
                                    if self.mondo_mappings[disgenet_dbs_to_mondo_dbs[db]].get(map_v):
                                        diseaseid2 = self.mondo_mappings[disgenet_dbs_to_mondo_dbs[db]].get(map_v)
                                        break

                    if diseaseid1 and diseaseid2:
                        df_list.append((diseaseid1, diseaseid2, round(dda["jaccard_variants"], 3),))

            disgenet_dda_variant_df = pd.DataFrame(df_list, columns=["disease_id1", "disease_id2", "disgenet_jaccard_variants_score"])
            disgenet_dda_variant_df["source"] = "Disgenet Disease-Disease Variant"
            
            disgenet_dda_variant_df.sort_values(by="disgenet_jaccard_variants_score", ascending=False, ignore_index=True, inplace=True)
            disgenet_dda_variant_df.drop_duplicates(subset=["disease_id1", "disease_id2"], ignore_index=True, inplace=True)
            
            t1 = time()
            print(f"Disgenet disease-disease data is processed in {round((t1-t0) / 60, 2)} mins")
            
            return disgenet_dda_gene_df, disgenet_dda_variant_df
        else:
            print("Processing Disgenet disease-disease data")
            t0 = time()
            
            df_list = []
            # DISEASE-DISEASE BY GENE
            for dda in self.disgenet_dda_gene:
                if round(dda.jaccard_genes, 3) != 0.0:
                    diseaseid1 =  self.mondo_mappings["UMLS"].get(dda.diseaseid1)
                    diseaseid2 = self.mondo_mappings["UMLS"].get(dda.diseaseid2)
                    
                    if not diseaseid1:
                        if self.disgenet_id_mappings_dict.get(dda.diseaseid1):
                            if self.disgenet_id_mappings_dict.get(dda.diseaseid1).get("MONDO"):                    
                                diseaseid1 = "MONDO:" + self.disgenet_id_mappings_dict.get(dda.diseaseid1).get("MONDO")
                            else:
                                map_dict = self.disgenet_id_mappings_dict.get(dda.diseaseid1)
                                for db, map_v in map_dict.items():
                                    if self.mondo_mappings[disgenet_dbs_to_mondo_dbs[db]].get(map_v):
                                        diseaseid1 = self.mondo_mappings[disgenet_dbs_to_mondo_dbs[db]].get(map_v)
                                        break
                                        
                    if not diseaseid2:
                        if self.disgenet_id_mappings_dict.get(dda.diseaseid2):
                            if self.disgenet_id_mappings_dict.get(dda.diseaseid2).get("MONDO"):
                                diseaseid2 = "MONDO:" + self.disgenet_id_mappings_dict.get(dda.diseaseid2).get("MONDO")
                            else:
                                map_dict = self.disgenet_id_mappings_dict.get(dda.diseaseid2)
                                for db, map_v in map_dict.items():
                                    if self.mondo_mappings[disgenet_dbs_to_mondo_dbs[db]].get(map_v):
                                        diseaseid2 = self.mondo_mappings[disgenet_dbs_to_mondo_dbs[db]].get(map_v)
                                        break
                                        
                    if diseaseid1 and diseaseid2:
                        df_list.append((diseaseid1, diseaseid2, round(dda.jaccard_genes, 3),))
            
            disgenet_dda_gene_df = pd.DataFrame(df_list, columns=["disease_id1", "disease_id2", "disgenet_jaccard_genes_score"])
            disgenet_dda_gene_df["source"] = "Disgenet Disease-Disease Gene"
            
            disgenet_dda_gene_df.sort_values(by="disgenet_jaccard_genes_score", ascending=False, ignore_index=True, inplace=True)
            disgenet_dda_gene_df.drop_duplicates(subset=["disease_id1", "disease_id2"], ignore_index=True, inplace=True)
            
            df_list = []
            # DISEASE-DISEASE BY VARIANT
            for dda in self.disgenet_dda_variant:
                if round(dda.jaccard_variants, 3) != 0.0:
                    diseaseid1 =  self.mondo_mappings["UMLS"].get(dda.diseaseid1)
                    diseaseid2 = self.mondo_mappings["UMLS"].get(dda.diseaseid2)
                    
                    if not diseaseid1:
                        if self.disgenet_id_mappings_dict.get(dda.diseaseid1):
                            if self.disgenet_id_mappings_dict.get(dda.diseaseid1).get("MONDO"):                    
                                diseaseid1 = "MONDO:" + self.disgenet_id_mappings_dict.get(dda.diseaseid1).get("MONDO")
                            else:
                                map_dict = self.disgenet_id_mappings_dict.get(dda.diseaseid1)
                                for db, map_v in map_dict.items():
                                    if self.mondo_mappings[disgenet_dbs_to_mondo_dbs[db]].get(map_v):
                                        diseaseid1 = self.mondo_mappings[disgenet_dbs_to_mondo_dbs[db]].get(map_v)
                                        break
                                        
                    if not diseaseid2:
                        if self.disgenet_id_mappings_dict.get(dda.diseaseid2):
                            if self.disgenet_id_mappings_dict.get(dda.diseaseid2).get("MONDO"):
                                diseaseid2 = "MONDO:" + self.disgenet_id_mappings_dict.get(dda.diseaseid2).get("MONDO")
                            else:
                                map_dict = self.disgenet_id_mappings_dict.get(dda.diseaseid2)
                                for db, map_v in map_dict.items():
                                    if self.mondo_mappings[disgenet_dbs_to_mondo_dbs[db]].get(map_v):
                                        diseaseid2 = self.mondo_mappings[disgenet_dbs_to_mondo_dbs[db]].get(map_v)
                                        break
                                        
                    if diseaseid1 and diseaseid2:
                        df_list.append((diseaseid1, diseaseid2, round(dda.jaccard_variants, 3),))
                        
            disgenet_dda_variant_df = pd.DataFrame(df_list, columns=["disease_id1", "disease_id2", "disgenet_jaccard_variants_score"])
            disgenet_dda_variant_df["source"] = "Disgenet Disease-Disease Variant"
            
            disgenet_dda_variant_df.sort_values(by="disgenet_jaccard_variants_score", ascending=False, ignore_index=True, inplace=True)
            disgenet_dda_variant_df.drop_duplicates(subset=["disease_id1", "disease_id2"], ignore_index=True, inplace=True)
            
            t1 = time()
            print(f"Disgenet disease-disease data is processed in {round((t1-t0) / 60, 2)} mins")
            
            return disgenet_dda_gene_df, disgenet_dda_variant_df
            
    def merge_disease_drug_edge_data(self) -> pd.DataFrame:
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
            
    def get_nodes(self, label="disease") -> list:
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
