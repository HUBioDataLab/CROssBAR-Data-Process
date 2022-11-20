import os
import sys
import pandas as pd
import numpy as np
import time
import collections
import argparse
from pathlib import Path
from time import time

from pypath.inputs import intact
from pypath.inputs import string
from pypath.inputs import biogrid
from pypath.share import curl
from pypath.inputs import uniprot
from pypath.utils import mapping

parser = argparse.ArgumentParser(description='CROssBAR v2 PPI data retrieval arguments')

parser.add_argument(
    '--output_dir',
    type=str,
    required=True)

parser.add_argument(
    '--n_rows_in_file',
    type=int,
    default=100000,
    help='approximate number of rows in each splitted output csvs (except string)')

class PPI_data:
    def __init__(self, output_dir, n_rows_in_file):
        Path(output_dir).mkdir(parents=True, exist_ok=True)
        self.output_dir = output_dir
        self.n_rows_in_file = n_rows_in_file
        self.swissprots = list(uniprot._all_uniprots("*", True))

    def export_dataframe(self, dataframe, data_label):
        output_base =  os.path.join(self.output_dir, data_label)
        Path(output_base).mkdir(parents=True, exist_ok=True)
        n_chunks = round(len(dataframe) / self.n_rows_in_file)
        for id, chunk in  enumerate(np.array_split(dataframe, n_chunks)):
            chunk.to_csv(os.path.join(output_base, f"crossbar_ppi_data_{data_label}_{id+1}.csv"), index=False)

        return output_base

    def intact_process(self):
        print("Started downloading IntAct data")

        t0 = time()

        with curl.cache_off():
            intact_ints = intact.intact_interactions(miscore=0, organism=None, complex_expansion=True, only_proteins=True)

        t1 = time()

        print(f'IntAct data is downloaded in {round((t1-t0) / 60, 2)} mins, now started processing')

        intact_df = pd.DataFrame.from_records(intact_ints, columns=intact_ints[0]._fields)
        
        # turn list columns to string
        for list_column in ["pubmeds", "methods", "interaction_types"]:
            intact_df[list_column] = [';'.join(map(str, l)) for l in intact_df[list_column]]
        
        intact_df["source"] = "IntAct"
        intact_df = intact_df[['source', 'id_a', 'id_b', 'pubmeds', 'mi_score', 'methods',  'interaction_types']]
        intact_df.columns = ['source', 'uniprot_a', 'uniprot_b', 'intact_pubmed_id', 'intact_score', 'intact_methods', 'intact_interaction_types']
        
        # drop rows if uniprot_a or uniprot_b is not a swiss-prot protein
        intact_df = intact_df[(intact_df["uniprot_a"].isin(self.swissprots)) & (intact_df["uniprot_b"].isin(self.swissprots))]
        intact_df.reset_index(drop=True, inplace=True)
        
        # assing pubmed ids that contain unassigned to NaN value 
        intact_df["intact_pubmed_id"].loc[intact_df["intact_pubmed_id"].astype(str).str.contains("unassigned", na=False)] = np.nan
        
        # drop duplicates if same a x b pair exists multiple times 
        # keep the pair with the highest score and collect pubmed ids of duplicated a x b pairs in that pair's pubmed id column
        # if a x b pair has same interaction type with b x a pair, drop b x a pair
        intact_df.sort_values(by=['intact_score'], ascending=False, inplace=True)
        intact_df_unique = intact_df.dropna(subset=["uniprot_a", "uniprot_b"]).reset_index(drop=True)
        intact_df_unique = intact_df_unique.groupby(["uniprot_a", "uniprot_b"], sort=False, as_index=False).aggregate({"source":"first", "uniprot_a":"first", "uniprot_b":"first", 
                                                    "intact_pubmed_id": lambda x: "|".join([str(e) for e in set(x.dropna())]),
                                                   "intact_score":"first", "intact_methods":"first", 
                                                    "intact_interaction_types":"first"})
        intact_df_unique["intact_pubmed_id"].replace("", np.nan, inplace=True) # replace empty string with NaN
        intact_df_unique = intact_df_unique[~intact_df_unique[["uniprot_a", "uniprot_b", "intact_interaction_types"]].apply(frozenset, axis=1).duplicated()].reset_index(drop=True)
        
        intact_output_base = self.export_dataframe(intact_df_unique, "intact")

        t2 = time()
        print(f'IntAct data is processed and written in {round((t2-t1) / 60, 2)} mins: {intact_output_base}')

        self.final_intact_ints = intact_df_unique

            

    def biogrid_process(self):
        t0 = time()
        print("Started downloading BioGRID data")

        with curl.cache_off():
            biogrid_ints = biogrid.biogrid_all_interactions(None, 9999999999, False)

        t1 = time()
        
        print(f'BioGRID data is downloaded in {round((t1-t0) / 60, 2)} mins, now started processing')
            
        biogrid_df = pd.DataFrame.from_records(biogrid_ints, columns=biogrid_ints[0]._fields)

        # biogrid id (gene symbols) to uniprot id mapping
        biogrid_df['partner_a'] = biogrid_df['partner_a'].str.upper()
        biogrid_df['partner_b'] = biogrid_df['partner_b'].str.upper()
        with curl.cache_off():
            uniprot_to_gene = uniprot.uniprot_data("genes", "*", True)
            uniprot_to_tax = uniprot.uniprot_data("organism-id", "*", True)

        gene_to_uniprot = collections.defaultdict(list)
        for k,v in uniprot_to_gene.items():
            for gene in v.split():
                gene_to_uniprot[gene.upper()].append(k)

        prot_a_uniprots = []
        for prot, tax in zip(biogrid_df['partner_a'], biogrid_df['tax_a']):
            uniprot_id_a = ";".join([_id for _id in gene_to_uniprot[prot] if tax == uniprot_to_tax[_id]])
            prot_a_uniprots.append(uniprot_id_a)

        prot_b_uniprots = []
        for prot, tax in zip(biogrid_df['partner_b'], biogrid_df['tax_b']):
            uniprot_id_b = ";".join([_id for _id in gene_to_uniprot[prot] if tax == uniprot_to_tax[_id]])
            prot_b_uniprots.append(uniprot_id_b)

        biogrid_df["uniprot_a"] = prot_a_uniprots
        biogrid_df["uniprot_b"] = prot_b_uniprots
        
        # drop rows that have semicolon (";")
        biogrid_df.drop(biogrid_df[(biogrid_df["uniprot_a"].str.contains(";")) | (biogrid_df["uniprot_b"].str.contains(";"))].index, axis=0, inplace=True)
        biogrid_df.reset_index(drop=True, inplace=True)
        
        biogrid_df["source"] = "BioGRID"
        biogrid_df = biogrid_df[['source', 'uniprot_a', 'uniprot_b', 'partner_a', 'partner_b', 'pmid', 'experimental_system',  'experimental_system_type', 'tax_a', 'tax_b']]
        biogrid_df.columns = ['source', 'uniprot_a', 'uniprot_b', 'biogrid_partner_a', 'biogrid_partner_b', 'biogrid_pubmed_id', 'biogrid_experimental_system', 'biogrid_experimental_system_type', 'biogrid_tax_a', 'biogrid_tax_b']
        
        # drop duplicates if same a x b pair exists multiple times 
        # keep the first pair and collect pubmed ids of duplicated a x b pairs in that pair's pubmed id column
        # if a x b pair has same experimental system type with b x a pair, drop b x a pair
        biogrid_df_unique = biogrid_df.dropna(subset=["uniprot_a", "uniprot_b"]).reset_index(drop=True)
        biogrid_df_unique.groupby(["uniprot_a", "uniprot_b"], sort=False, as_index=False).aggregate({"source":"first", "uniprot_a":"first",
                                                                                             "uniprot_b":"first", "biogrid_partner_a":"first",
                                                                                             "biogrid_partner_b":"first", 
                                                                                             "biogrid_pubmed_id":lambda x: "|".join([str(e) for e in set(x.dropna())]),
                                                                                             "biogrid_experimental_system":"first", "biogrid_experimental_system_type":"first",
                                                                                             "biogrid_tax_a":"first", "biogrid_tax_b":"first"})
        biogrid_df_unique = biogrid_df_unique[~biogrid_df_unique[["uniprot_a", "uniprot_b", "biogrid_experimental_system"]].apply(frozenset, axis=1).duplicated()].reset_index(drop=True)
        
        # drop rows if uniprot_a or uniprot_b is not a swiss-prot protein
        biogrid_df_unique = biogrid_df_unique[(biogrid_df_unique["uniprot_a"].isin(self.swissprots)) & (biogrid_df_unique["uniprot_b"].isin(self.swissprots))]
        biogrid_df_unique.reset_index(drop=True, inplace=True)
        
        biogrid_output_base = self.export_dataframe(biogrid_df_unique, "biogrid")

        t2 = time()
        print(f'BioGRID data is processed and written in {round((t2-t1) / 60, 2)} mins: {biogrid_output_base}')

        self.final_biogrid_ints = biogrid_df_unique


    def string_process(self):
        string_output_base =  os.path.join(self.output_dir, "string")
        Path(string_output_base).mkdir(parents=True, exist_ok=True)

        log_path = os.path.join(string_output_base, "skipped_tax_in_string.log")
        logfile = open(log_path, 'w')

        with curl.cache_off():
            uniprot_tax=uniprot.uniprot_data("organism-id", "*")
        tax_ids=set(uniprot_tax.values())

        #map string ids to swissprot ids
        with curl.cache_off():
            uniprot_to_string = uniprot.uniprot_data("database(STRING)", "*", True)
        string_to_uniprot = collections.defaultdict(list)
        for k,v in uniprot_to_string.items():
            for string_id in list(filter(None, v.split(";"))):
                string_to_uniprot[string_id.split(".")[1]].append(k)

        print("Started downloading and processing STRING data")
        print(f"Check the progress from {log_path} for the rest of the process")

        # download organism specific string data and collect into a single list
        ranges = range(0, len(tax_ids), 500)
        organisms_with_string_ints_total = 0

        self.final_string_ints = pd.DataFrame(columns = ['source', 'uniprot_a', 'uniprot_b', 'string_partner_a', 'string_partner_b', 'string_combined_score', 'string_physical_combined_score'])
        
        t0 = time()
        for idx, r in enumerate(ranges):
            t0_0 = time()
            tax_chunk = list(tax_ids)[r:r + 500]
            string_ints = []
            organisms_with_string_ints = 0

            for tax in tax_chunk:
                try: 
                    with curl.cache_off():
                        organism_string_ints = [i for i in string.string_links_interactions(ncbi_tax_id=int(tax), score_threshold="high_confidence")]
                    string_ints.extend(organism_string_ints)
                    logfile.write(f"{tax} done\n")
                    logfile.flush()
                    organisms_with_string_ints += 1
                    organisms_with_string_ints_total += 1
                except Exception as e: # there are no interactions for this tax
                    pass
            
            if string_ints:
                string_df = pd.DataFrame.from_records(string_ints, columns=string_ints[0]._fields)

                # map string ids to uniprot ids
                prot_a_uniprots = []
                for protein in string_df['protein_a']:
                    id_a= (
                        ";".join(string_to_uniprot[protein])
                            if protein in string_to_uniprot else
                        None
                    )
                    prot_a_uniprots.append(id_a)

                prot_b_uniprots = []
                for protein in string_df['protein_b']:
                    id_b= (
                        ";".join(string_to_uniprot[protein])
                            if protein in string_to_uniprot else
                        None
                    )
                    prot_b_uniprots.append(id_b)

                string_df["uniprot_a"] = prot_a_uniprots
                string_df["uniprot_b"] = prot_b_uniprots

                # drop duplicates if same a x b pair exists in b x a format
                # keep the one with the highest combined score
                string_df.sort_values(by=['combined_score'], ascending=False, inplace=True)
                string_df_unique = string_df.dropna(subset=["uniprot_a", "uniprot_b"]).drop_duplicates(subset=["uniprot_a", "uniprot_b"], keep="first").reset_index(drop=True)
                string_df_unique = string_df_unique[~string_df_unique[["uniprot_a", "uniprot_b", "combined_score"]].apply(frozenset, axis=1).duplicated()].reset_index(drop=True)

                string_df_unique["source"] = "STRING"
                string_df_unique = string_df_unique[['source', 'uniprot_a', 'uniprot_b', 'protein_a', 'protein_b', 'combined_score', 'physical_combined_score']]
                string_df_unique.columns = ['source', 'uniprot_a', 'uniprot_b', 'string_partner_a', 'string_partner_b', 'string_combined_score', 'string_physical_combined_score']

                # filter with swissprot ids
                string_df_unique = string_df_unique[(string_df_unique["uniprot_a"].isin(self.swissprots)) & (string_df_unique["uniprot_b"].isin(self.swissprots))]
                string_df_unique.reset_index(drop=True, inplace=True)
                
                t0_1 = time()
                string_output = os.path.join(string_output_base, f"crossbar_ppi_data_string_{idx+1}.csv")
                string_df_unique.to_csv(string_output, index=False)

                logfile.write(f'STRING data for {organisms_with_string_ints}/{len(tax_chunk)} organisms in the batch {idx+1} is processed and written in {round((t0_1-t0_0) / 60, 2)} mins: {string_output_base}\n')
                logfile.flush()

            else:
                logfile.write(f'There are no STRING interactions for {len(tax_chunk)} organisms in the batch {idx+1}\n')
                logfile.flush()

            self.final_string_ints = pd.concat([self.final_string_ints, string_df_unique])

        t1 = time()
        logfile.write(f'STRING data is processed and written in {round((t1-t0) / 60, 2)} mins for total {organisms_with_string_ints_total} organisms')
        logfile.flush()


    def merge_all(self):

        print("started merging interactions from all sources")
        temp = self.final_intact_ints.merge(self.final_biogrid_ints, on=['uniprot_a', 'uniprot_b'], how='outer' )

        print("merged intact and biogrid")
        temp = temp.merge(self.final_string_ints, on=['uniprot_a', 'uniprot_b'], how='outer')
        print("merged all")

        temp["source_all"] = temp[["source", "source_x", "source_y"]].apply(lambda x: ';'.join(x.dropna()), axis=1)
        temp.drop(['source_x', 'source_y', 'source'], axis=1, inplace=True)
        all_output = os.path.join(self.output_dir, "all_ppis.csv")

        temp.to_csv(all_output, index=False)
        
        all_output_base = self.export_dataframe(temp, "all_ppi_splitted")
      
    def merge_mall(self):
        # select and define fields of intact dataframe
        intact_refined_df_selected_features = self.final_intact_ints.rename(columns={"intact_methods":"method", "intact_interaction_types":"interaction_type", "intact_pubmed_id":"pubmed_id"})
        intact_refined_df_selected_features = intact_refined_df_selected_features.reindex(columns=["source", "uniprot_a", "uniprot_b", "pubmed_id", "method", "interaction_type", "intact_score"])
        intact_refined_df_selected_features["string_combined_score"] = np.nan
        intact_refined_df_selected_features["string_physical_combined_score"] = np.nan
        
        # select and define fields of biogrid dataframe
        biogrid_refined_df_selected_features = self.final_biogrid_ints.drop(columns=["biogrid_partner_a", "biogrid_partner_b", "biogrid_experimental_system_type", "biogrid_tax_a", "biogrid_tax_b"])
        biogrid_refined_df_selected_features = biogrid_refined_df_selected_features.rename(columns={"biogrid_experimental_system":"method", "biogrid_pubmed_id":"pubmed_id"})
        biogrid_refined_df_selected_features = biogrid_refined_df_selected_features.reindex(columns=["source","uniprot_a","uniprot_b", "pubmed_id", "method"])
        biogrid_refined_df_selected_features[["interaction_type", "intact_score", "string_combined_score", "string_physical_combined_score"]] = np.nan
        
        # select and define fields of string dataframe
        string_refined_df_selected_features = self.final_string_ints.drop(columns=["string_partner_a","string_partner_b"])
        string_refined_df_selected_features[["pubmed_id", "method", "interaction_type", "intact_score"]] = np.nan
        string_refined_df_selected_features = string_refined_df_selected_features.reindex(columns=["source", "uniprot_a", "uniprot_b",
        "pubmed_id", "method", "interaction_type", "intact_score", "string_combined_score", "string_physical_combined_score"])
        
        # merge intact and biogrid
        intact_plus_biogrid_selected_features_df = pd.merge(intact_refined_df_selected_features, biogrid_refined_df_selected_features,
                                                   on=["uniprot_a", "uniprot_b"], how="outer")
        
        # merge source_x and source_y columns
        intact_plus_biogrid_selected_features_df["source"] = intact_plus_biogrid_selected_features_df[["source_x", "source_y"]].apply(lambda x: '|'.join(x.dropna()), axis=1)
        
        # drop redundant columns
        intact_plus_biogrid_selected_features_df.drop(columns=["source_x", "source_y"], inplace=True)
        
        # merge pubmed_id_x and pubmed_id_y columns
        intact_plus_biogrid_selected_features_df["pubmed_id"] = intact_plus_biogrid_selected_features_df[["pubmed_id_x", "pubmed_id_y"]].apply(lambda x: int(x.dropna().tolist()[0]) if len(x.dropna().tolist())>0 else np.nan, axis=1)
        
        # drop redundant columns
        intact_plus_biogrid_selected_features_df.drop(columns=["pubmed_id_x", "pubmed_id_y"], inplace=True)
        
        # merge method_x and method_y columns
        intact_plus_biogrid_selected_features_df["method"] = intact_plus_biogrid_selected_features_df[["method_x", "method_y"]].apply(lambda x: x.dropna().tolist()[0] if len(x.dropna().tolist())>0 else np.nan, axis=1)
        
        # drop redundant columns
        intact_plus_biogrid_selected_features_df.drop(columns=["method_x", "method_y"], inplace=True)
        intact_plus_biogrid_selected_features_df.drop(columns=["interaction_type_y", "intact_score_y", 
                                                       "string_combined_score_y", "string_physical_combined_score_y"],
                                             inplace=True)
        
        # rename and reorder columns
        intact_plus_biogrid_selected_features_df.rename(columns={"interaction_type_x":"interaction_type", "intact_score_x":"intact_score",
                                                        "string_combined_score_x":"string_combined_score",
                                                        "string_physical_combined_score_x":"string_physical_combined_score"},
                                               inplace=True)        
        intact_plus_biogrid_selected_features_df = intact_plus_biogrid_selected_features_df.reindex(columns=['source', 'uniprot_a', 'uniprot_b', 'pubmed_id', 
                                                          'method', 'interaction_type', 'intact_score', 
                                                          'string_combined_score', 'string_physical_combined_score'])
        
        
        # merge intact+biogrid with string
        all_selected_features_df = pd.merge(intact_plus_biogrid_selected_features_df, string_refined_df_selected_features, on=["uniprot_a", "uniprot_b"], how="outer")
        
        # merge source_x and source_y columns
        all_selected_features_df["source"] = all_selected_features_df[["source_x", "source_y"]].apply(lambda x: '|'.join(x.dropna()), axis=1)
        
        # drop redundant columns
        all_selected_features_df.drop(columns=["source_x", "source_y"], inplace=True)
        all_selected_features_df.drop(columns=["string_combined_score_x", "string_physical_combined_score_x", "pubmed_id_y",
                                      "method_y", "interaction_type_y", "intact_score_y"], inplace=True)
        
        
        # rename and reorder columns
        all_selected_features_df.rename(columns={"interaction_type_x":"interaction_type", "intact_score_x":"intact_score",
                                        "pubmed_id_x":"pubmed_id", "method_x":"method", 
                                         "string_combined_score_y":"string_combined_score",
                                        "string_physical_combined_score_y":"string_physical_combined_score"},
                                inplace=True)
        all_selected_features_df = all_selected_features_df.reindex(columns=['source', 'uniprot_a', 'uniprot_b', 'pubmed_id', 
                                                          'method', 'interaction_type', 'intact_score', 
                                                          'string_combined_score', 'string_physical_combined_score'])
        
        
        # during the merging of 2 pubmed_id columns, it changes datatype from int to float. So it needs to be reverted
        def float_to_int(element):
            if "." in str(element):                
                dot_index = str(element).index(".")
                element = str(element)[:dot_index]
                return element
            else:
                return element
        
        # first make their datatype as string
        all_selected_features_df["pubmed_id"] = all_selected_features_df["pubmed_id"].astype(str, errors="ignore")
        all_selected_features_df["string_physical_combined_score"] = all_selected_features_df["string_physical_combined_score"].astype(str, errors="ignore")
        
        # then revert back them
        all_selected_features_df["pubmed_id"] = all_selected_features_df["pubmed_id"].apply(float_to_int)
        all_selected_features_df["string_physical_combined_score"] = all_selected_features_df["string_physical_combined_score"].apply(float_to_int)
        
        return all_selected_features_df      
        
        
if __name__ == "__main__":
    
    args = parser.parse_args()
    t0 = time() 
    ppi_downloader = PPI_data(output_dir=args.output_dir, n_rows_in_file=args.n_rows_in_file) 
    intact_ints = ppi_downloader.intact_process()
    biogrid_ints =ppi_downloader.biogrid_process()
    string_ints = ppi_downloader.string_process()
    ppi_downloader.merge_all()   
    
    #check_df = ppi_downloader.merge_mall()
    
    t1 = time()
    print(f'Done in total {round((t1-t0) / 60, 2 )} mins')
