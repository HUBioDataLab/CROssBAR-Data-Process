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

        # drop duplicates if same a x b pair exists in b x a format
        # keep the one with the highest score
        # keep both if their interaction types are different
        intact_df.sort_values(by=['mi_score'], ascending=False, inplace=True)
        intact_df_unique = intact_df.dropna(subset=["id_a", "id_b"]).drop_duplicates(subset=["id_a", "id_b"], keep="first").reset_index(drop=True)
        intact_df_unique = intact_df_unique[~intact_df_unique[["id_a", "id_b", "interaction_types"]].apply(frozenset, axis=1).duplicated()].reset_index(drop=True)

        intact_df_unique["source"] = "IntAct"
        intact_df_unique = intact_df_unique[['source', 'id_a', 'id_b', 'pubmeds', 'mi_score', 'methods',  'interaction_types']]
        intact_df_unique.columns = ['source', 'uniprot_a', 'uniprot_b', 'intact_pubmed_id', 'intact_score', 'intact_methods', 'intact_interaction_types']
        
        # drop rows if uniprot_a or uniprot_b is not a swiss-prot protein
        intact_df_unique = intact_df_unique[(intact_df_unique["uniprot_a"].isin(self.swissprots)) & (intact_df_unique["uniprot_b"].isin(self.swissprots))]

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
            uniprot_id_a = ";".join([id for id in gene_to_uniprot[prot] if tax == uniprot_to_tax[id]])
            prot_a_uniprots.append(uniprot_id_a)

        prot_b_uniprots = []
        for prot, tax in zip(biogrid_df['partner_b'], biogrid_df['tax_b']):
            uniprot_id_b = ";".join([id for id in gene_to_uniprot[prot] if tax == uniprot_to_tax[id]])
            prot_b_uniprots.append(uniprot_id_b)

        biogrid_df["uniprot_a"] = prot_a_uniprots
        biogrid_df["uniprot_b"] = prot_b_uniprots

        # drop duplicates if same a x b pair exists in b x a format
        # keep both if their pubmed ids are different
        biogrid_df_unique = biogrid_df.dropna(subset=["uniprot_a", "uniprot_b"]).drop_duplicates(subset=["uniprot_a", "uniprot_b"], keep="first").reset_index(drop=True)
        biogrid_df_unique = biogrid_df_unique[~biogrid_df_unique[["uniprot_a", "uniprot_b", "pmid"]].apply(frozenset, axis=1).duplicated()].reset_index(drop=True)

        biogrid_df_unique["source"] = "BioGRID"
        biogrid_df_unique = biogrid_df_unique[['source', 'uniprot_a', 'uniprot_b', 'partner_a', 'partner_b', 'pmid', 'experimental_system',  'experimental_system_type', 'tax_a', 'tax_b']]
        biogrid_df_unique.columns = ['source', 'uniprot_a', 'uniprot_b', 'biogrid_partner_a', 'biogrid_partner_b', 'biogrid_pubmed_id', 'biogrid_experimental_system', 'biogrid_experimental_system_type', 'biogrid_tax_a', 'biogrid_tax_b']

        # drop rows if uniprot_a or uniprot_b is not a swiss-prot protein
        biogrid_df_unique = biogrid_df_unique[(biogrid_df_unique["uniprot_a"].isin(self.swissprots)) & (biogrid_df_unique["uniprot_b"].isin(self.swissprots))]

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

if __name__ == "__main__":
    
    args = parser.parse_args()
    t0 = time() 
    ppi_downloader = PPI_data(output_dir=args.output_dir, n_rows_in_file=args.n_rows_in_file) 
    intact_ints = ppi_downloader.intact_process()
    biogrid_ints =ppi_downloader.biogrid_process()
    string_ints = ppi_downloader.string_process()
    ppi_downloader.merge_all()
    t1 = time()
    print(f'Done in total {round((t1-t0) / 60, 2 )} mins')