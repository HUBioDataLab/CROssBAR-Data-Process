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
from pypath.share import curl, settings
from pypath.inputs import uniprot

from tqdm import tqdm # progress bar

from biocypher._logger import logger
from contextlib import ExitStack

logger.debug(f"Loading module {__name__}.")

# TODO: delete argument parser and add output_dir to config file as a setting; maybe 'ppi_output_dir'
parser = argparse.ArgumentParser(description='CROssBAR v2 PPI data retrieval arguments')

parser.add_argument(
    '--output_dir',
    type=str,
    required=True)


class PPI_data:
    def __init__(self, output_dir, export_csvs = False, split_output = False, cache=False, debug=False, retries=6):
        """
            Args:
                export_csvs: Flag for whether or not create csvs of outputs
                split_csvs: whether or not to split output csv files to multiple parts
                TODO: add n_rows_in_file setting to config file
                cache: if True, it uses the cached version of the data, otherwise
                forces download.
                debug: if True, turns on debug mode in pypath.
                retries: number of retries in case of download error.
        """
        
        self.export_csvs = export_csvs
        self.split_output = split_output
        self.swissprots = list(uniprot._all_uniprots("*", True))
        self.cache = cache
        self.debug = debug
        self.retries = retries
        
        if export_csvs:
            Path(output_dir).mkdir(parents=True, exist_ok=True)
            self.output_dir = output_dir

    
    def export_dataframe(self, dataframe, data_label):
        Path(self.output_dir).mkdir(parents=True, exist_ok=True)

        # TODO: activate this block after adding n_rows_in_file setting to config file
        # if self.split_output:
        #     n_chunks = round(len(dataframe) / n_rows_in_file)
            # output_path =  os.path.join(self.output_dir, data_label)
            # for id, chunk in  enumerate(np.array_split(dataframe, n_chunks)):
            #     chunk.to_csv(os.path.join(output_path, f"crossbar_ppi_data_{data_label}_{id+1}.csv"), index=False)
        # else:
        output_path = os.path.join(self.output_dir, f"crossbar_ppi_data_{data_label}.csv")
        dataframe.to_csv(output_path, index=False)

        return output_path

    def download_intact_data(self):
        """
        Wrapper function to download IntAct data using pypath; used to access
        settings.
            
        To do: Make arguments of intact.intact_interactions selectable for user.
        """
                     
        logger.debug("Started downloading IntAct data")
        t0 = time()
        with ExitStack() as stack:                         
            stack.enter_context(settings.context(retries=self.retries))
            
            if self.debug:                
                stack.enter_context(curl.debug_on())
            if not self.cache:
                stack.enter_context(curl.cache_off())

            self.intact_ints = intact.intact_interactions(miscore=0, organism=None, complex_expansion=True, only_proteins=True)
        t1 = time()
        
        logger.info(f'IntAct data is downloaded in {round((t1-t0) / 60, 2)} mins')


    def intact_process(self, selected_fields=['source', 'id_a', 'id_b', 'pubmeds', 'mi_score', 'methods',  'interaction_types']):
        """
        Processor function for IntAct data. It drops duplicate and reciprocal duplicate protein pairs and collects pubmed ids of duplicated pairs. Also, it filters
        protein pairs found in swissprot.
        
        Args:
            selected_fields: fields to be used in the data.
            
        """
        logger.debug("Started processing IntAct data")
        t1 = time()
                         
        # create dataframe            
        intact_df = pd.DataFrame.from_records(self.intact_ints, columns=self.intact_ints[0]._fields)
        
        # turn list columns to string
        for list_column in ["pubmeds", "methods", "interaction_types"]:
            intact_df[list_column] = [';'.join(map(str, l)) for l in intact_df[list_column]]
        
        intact_df.fillna(value=np.nan, inplace=True)
        
        # add source database info
        intact_df["source"] = "IntAct"
        # filter selected fields
        intact_df = intact_df[selected_fields]
        # rename columns
        intact_df.columns = ['source', 'uniprot_a', 'uniprot_b', 'pubmed_id', 'intact_score', 'method', 'interaction_type']
        
        # drop rows if uniprot_a or uniprot_b is not a swiss-prot protein
        intact_df = intact_df[(intact_df["uniprot_a"].isin(self.swissprots)) & (intact_df["uniprot_b"].isin(self.swissprots))]
        intact_df.reset_index(drop=True, inplace=True)
        
        # assing pubmed ids that contain unassigned to NaN value 
        intact_df["pubmed_id"].loc[intact_df["pubmed_id"].astype(str).str.contains("unassigned", na=False)] = np.nan
        
        # drop duplicates if same a x b pair exists multiple times 
        # keep the pair with the highest score and collect pubmed ids of duplicated a x b pairs in that pair's pubmed id column
        # if a x b pair has same interaction type with b x a pair, drop b x a pair
        intact_df.sort_values(by=['intact_score'], ascending=False, inplace=True)
        intact_df_unique = intact_df.dropna(subset=["uniprot_a", "uniprot_b"]).reset_index(drop=True)
        intact_df_unique = intact_df_unique.groupby(["uniprot_a", "uniprot_b"], sort=False, as_index=False).aggregate({"source":"first", "uniprot_a":"first", "uniprot_b":"first", 
                                                    "pubmed_id": lambda x: "|".join([str(e) for e in set(x.dropna())]),
                                                   "intact_score":"first", "method":"first", 
                                                    "interaction_type":"first"})
        intact_df_unique["pubmed_id"].replace("", np.nan, inplace=True) # replace empty string with NaN
        intact_df_unique = intact_df_unique[~intact_df_unique[["uniprot_a", "uniprot_b", "interaction_type"]].apply(frozenset, axis=1).duplicated()].reset_index(drop=True)
        
        if self.export_csvs:
            intact_output_path = self.export_dataframe(intact_df_unique, "intact")
            logger.info(f'Final IntAct data is written: {intact_output_path}')

        self.final_intact_ints = intact_df_unique
        
        t2 = time()
        logger.info(f'IntAct data is processed in {round((t2-t1) / 60, 2)} mins')
                         
    def download_biogrid_data(self):
        """
        Wrapper function to download BioGRID data using pypath; used to access
        settings.
            
        To do: Make arguments of biogrid.biogrid_all_interactions selectable for user. 
        """
        
             
        logger.debug("Started downloading BioGRID data")
        t0 = time()

        with ExitStack() as stack:                         
            stack.enter_context(settings.context(retries=self.retries))
            
            if self.debug:                
                stack.enter_context(curl.debug_on())
            if not self.cache:
                stack.enter_context(curl.cache_off())

            # download biogrid data
            self.biogrid_ints = biogrid.biogrid_all_interactions(None, 9999999999, False)
                        
            # download these fields for mapping from gene symbol to uniprot id          
            self.uniprot_to_gene = uniprot.uniprot_data("genes", "*", True)
            self.uniprot_to_tax = uniprot.uniprot_data("organism-id", "*", True)
                    
        t1 = time()
        logger.info(f'BioGRID data is downloaded in {round((t1-t0) / 60, 2)} mins')
                         

    def biogrid_process(self, selected_fields=['source', 'uniprot_a', 'uniprot_b', 'pmid', 'experimental_system']):
        """
        Processor function for BioGRID data. It drops duplicate and reciprocal duplicate protein pairs and collects pubmed ids of duplicated pairs. In addition, it
        maps entries to uniprot ids using gene name and tax id information in the BioGRID data. Also, it filters protein pairs found in swissprot.
        
        Args:
            selected_fields: fields to be used in the data.            
        """
        
        logger.debug("Started processing BioGRID data")

        t1 = time()
                         
        # create dataframe          
        biogrid_df = pd.DataFrame.from_records(self.biogrid_ints, columns=self.biogrid_ints[0]._fields)

        # biogrid id (gene symbols) to uniprot id mapping
        biogrid_df['partner_a'] = biogrid_df['partner_a'].str.upper()
        biogrid_df['partner_b'] = biogrid_df['partner_b'].str.upper()
                         
        gene_to_uniprot = collections.defaultdict(list)
        for k,v in self.uniprot_to_gene.items():
            for gene in v.split():
                gene_to_uniprot[gene.upper()].append(k)

        prot_a_uniprots = []
        for prot, tax in zip(biogrid_df['partner_a'], biogrid_df['tax_a']):
            uniprot_id_a = (
                ";".join([_id for _id in gene_to_uniprot[prot] if tax == self.uniprot_to_tax[_id]])
                    if prot in gene_to_uniprot else
                None)
            prot_a_uniprots.append(uniprot_id_a)

        prot_b_uniprots = []
        for prot, tax in zip(biogrid_df['partner_b'], biogrid_df['tax_b']):
            uniprot_id_b = (
                ";".join([_id for _id in gene_to_uniprot[prot] if tax == self.uniprot_to_tax[_id]])
                    if prot in gene_to_uniprot else
                None)
            prot_b_uniprots.append(uniprot_id_b)

        biogrid_df["uniprot_a"] = prot_a_uniprots
        biogrid_df["uniprot_b"] = prot_b_uniprots
        
        biogrid_df.fillna(value=np.nan, inplace=True)
        
        # add source database info
        biogrid_df["source"] = "BioGRID"
        # filter selected fields
        biogrid_df = biogrid_df[selected_fields]
        # rename columns
        biogrid_df.columns = ['source', 'uniprot_a', 'uniprot_b', 'pubmed_id', 'method']
        
        # drop rows that have semicolon (";")
        biogrid_df.drop(biogrid_df[(biogrid_df["uniprot_a"].str.contains(";")) | (biogrid_df["uniprot_b"].str.contains(";"))].index, axis=0, inplace=True)
        biogrid_df.reset_index(drop=True, inplace=True)
        
        # drop rows if uniprot_a or uniprot_b is not a swiss-prot protein
        biogrid_df = biogrid_df[(biogrid_df["uniprot_a"].isin(self.swissprots)) & (biogrid_df["uniprot_b"].isin(self.swissprots))]
        biogrid_df.reset_index(drop=True, inplace=True)
        
        # drop duplicates if same a x b pair exists multiple times 
        # keep the first pair and collect pubmed ids of duplicated a x b pairs in that pair's pubmed id column
        # if a x b pair has same experimental system type with b x a pair, drop b x a pair
        biogrid_df_unique = biogrid_df.dropna(subset=["uniprot_a", "uniprot_b"]).reset_index(drop=True)
        biogrid_df_unique.groupby(["uniprot_a", "uniprot_b"], sort=False, as_index=False).aggregate({"source":"first", "uniprot_a":"first",
                                                                                             "uniprot_b":"first", 
                                                                                             "pubmed_id":lambda x: "|".join([str(e) for e in set(x.dropna())]),
                                                                                             "method":"first"})
        biogrid_df_unique["pubmed_id"].replace("", np.nan, inplace=True)
        biogrid_df_unique = biogrid_df_unique[~biogrid_df_unique[["uniprot_a", "uniprot_b", "method"]].apply(frozenset, axis=1).duplicated()].reset_index(drop=True)
        
        if self.export_csvs:
            biogrid_output_path = self.export_dataframe(biogrid_df_unique, "biogrid")
            logger.info(f'Final BioGRID data is written: {biogrid_output_path}')

        self.final_biogrid_ints = biogrid_df_unique
        
        t2 = time()
        logger.info(f'BioGRID data is processed in {round((t2-t1) / 60, 2)} mins')
    

    def download_string_data(self):
        """
        Wrapper function to download STRING data using pypath; used to access
        settings.
            
        To do: Make arguments of string.string_links_interactions selectable for user.
        """

        t0 = time()

        with ExitStack() as stack:
            
            stack.enter_context(settings.context(retries=self.retries))
            
            if self.debug:                
                stack.enter_context(curl.debug_on())
            if not self.cache:
                stack.enter_context(curl.cache_off())
             
            string_species = string.string_species()
            self.tax_ids = list(string_species.keys())
        
            # map string ids to swissprot ids
            uniprot_to_string = uniprot.uniprot_data("database(STRING)", "*", True)
            
            self.string_to_uniprot = collections.defaultdict(list)
            for k,v in uniprot_to_string.items():
                for string_id in list(filter(None, v.split(";"))):
                    self.string_to_uniprot[string_id.split(".")[1]].append(k)

            self.string_ints = []
            
            logger.debug("Started downloading STRING data")

            tax_ids_to_be_skipped = ['4565', ]
            for tax in tqdm(self.tax_ids):
                if tax in tax_ids_to_be_skipped:
                    continue
                else:
                    # remove proteins that does not have swissprot ids
                    organism_string_ints = [
                        i for i in string.string_links_interactions(ncbi_tax_id=int(tax), score_threshold="high_confidence")
                        if i.protein_a in self.string_to_uniprot and i.protein_b in self.string_to_uniprot
                        ]
                    
                    logger.debug(f"Downloaded STRING data with taxonomy id {str(tax)}")

                    if organism_string_ints:
                        self.string_ints.extend(organism_string_ints)
            
        t1 = time()
        logger.info(f'STRING data is downloaded in {round((t1-t0) / 60, 2)} mins')
                         

    def string_process(self, selected_fields=['source', 'uniprot_a', 'uniprot_b', 'combined_score', 'physical_combined_score']):
        """
        Processor function for STRING data. It drops duplicate and reciprocal duplicate protein pairs. In addition, it maps entries to uniprot ids 
        using crossreferences to STRING in the Uniprot data. Also, it filters protein pairs found in swissprot.
        
        Args:
            selected_fields: fields to be used in the data.            
        """
        
        logger.debug("Started processing STRING data")
        t1 = time()
                                                  
        # create dataframe
        string_df = pd.DataFrame.from_records(self.string_ints, columns=self.string_ints[0]._fields)
        
                         
        prot_a_uniprots = []
        for protein in string_df['protein_a']:            
            id_a= (";".join(self.string_to_uniprot[protein]))
                   # if protein in self.string_to_uniprot else None) 
                   # now that we filtered interactions in line 307, we should not get KeyError here
            prot_a_uniprots.append(id_a)
                         
        prot_b_uniprots = []
        for protein in string_df['protein_b']:            
            id_b= (";".join(self.string_to_uniprot[protein]))
                   # if protein in self.string_to_uniprot else None)
                   # now that we filtered interactions in line 307, we should not get KeyError here
            prot_b_uniprots.append(id_b)
                         
        string_df["uniprot_a"] = prot_a_uniprots
        string_df["uniprot_b"] = prot_b_uniprots
        
        string_df.fillna(value=np.nan, inplace=True)
        
        # add source database info
        string_df["source"] = "STRING"
        # filter selected fields
        string_df = string_df[selected_fields]
        # rename columns
        string_df.columns = ['source', 'uniprot_a', 'uniprot_b', 'string_combined_score', 'string_physical_combined_score']
        
        # filter with swissprot ids
        # we already filtered interactions in line 307, we can remove this part or keep it for a double check
        # string_df = string_df[(string_df["uniprot_a"].isin(self.swissprots)) & (string_df["uniprot_b"].isin(self.swissprots))]
        # string_df.reset_index(drop=True, inplace=True)
                         
        # drop duplicates if same a x b pair exists in b x a format
        # keep the one with the highest combined score
        string_df.sort_values(by=['string_combined_score'], ascending=False, inplace=True)
        string_df_unique = string_df.dropna(subset=["uniprot_a", "uniprot_b"]).drop_duplicates(subset=["uniprot_a", "uniprot_b"], keep="first").reset_index(drop=True)
        string_df_unique = string_df_unique[~string_df_unique[["uniprot_a", "uniprot_b", "string_combined_score"]].apply(frozenset, axis=1).duplicated()].reset_index(drop=True)
        
        t2 = time()
        logger.info(f'STRING data is processed in {round((t2-t1) / 60, 2)} mins')
        
        if self.export_csvs:
            string_output_path = self.export_dataframe(string_df_unique, "string")
            logger.info(f'Final STRING data is written: {string_output_path}')

        self.final_string_ints = string_df_unique

      
    def merge_mall(self):
        """
        Merge function for all 3 databases. Merges dataframe according to uniprot_a and uniprot_b (i.e., protein pairs) columns.
                  
        """
        
        t1 = time()
        logger.debug("started merging interactions from all 3 databases (IntAct, BioGRID, STRING)")
                         
        # reorder columns of intact dataframe
        intact_refined_df_selected_features = self.final_intact_ints
        intact_refined_df_selected_features = intact_refined_df_selected_features.reindex(columns=["source", "uniprot_a", "uniprot_b", "pubmed_id", "method", "interaction_type", "intact_score"])
        
        # reorder columns of biogrid dataframe
        biogrid_refined_df_selected_features = self.final_biogrid_ints
        biogrid_refined_df_selected_features = biogrid_refined_df_selected_features.reindex(columns=["source", "uniprot_a", "uniprot_b", "pubmed_id", "method"])
        
        # reorder columns of string dataframe
        string_refined_df_selected_features = self.final_string_ints
        string_refined_df_selected_features = string_refined_df_selected_features.reindex(columns=["source", "uniprot_a", "uniprot_b",
        "string_combined_score", "string_physical_combined_score"])
        
        # merge intact and biogrid
        intact_plus_biogrid_selected_features_df = pd.merge(intact_refined_df_selected_features, biogrid_refined_df_selected_features,
                                                   on=["uniprot_a", "uniprot_b"], how="outer")
        
        # merge source_x and source_y columns
        intact_plus_biogrid_selected_features_df["source"] = intact_plus_biogrid_selected_features_df[["source_x", "source_y"]].apply(lambda x: '|'.join(x.dropna()), axis=1)
        
        # drop redundant columns
        intact_plus_biogrid_selected_features_df.drop(columns=["source_x", "source_y"], inplace=True)
        
        def merge_pubmed_ids(elem):            
            if len(elem.dropna().tolist()) > 0:                
                new_list = []
                for e in elem.dropna().tolist():                    
                    if "|" in e:                        
                        new_list.extend(e.split("|"))
                    else:
                        new_list.append(e)
        
                return "|".join(list(set(new_list))) 
            else:
                return np.nan
        
        # merge pubmed_id_x and pubmed_id_y columns
        intact_plus_biogrid_selected_features_df["pubmed_id"] = intact_plus_biogrid_selected_features_df[["pubmed_id_x", "pubmed_id_y"]].apply(merge_pubmed_ids, axis=1)
        
        # drop redundant columns
        intact_plus_biogrid_selected_features_df.drop(columns=["pubmed_id_x", "pubmed_id_y"], inplace=True)
        
        # merge method_x and method_y columns
        intact_plus_biogrid_selected_features_df["method"] = intact_plus_biogrid_selected_features_df[["method_x", "method_y"]].apply(lambda x: x.dropna().tolist()[0] if len(x.dropna().tolist())>0 else np.nan, axis=1)
        
        # drop redundant columns
        intact_plus_biogrid_selected_features_df.drop(columns=["method_x", "method_y"], inplace=True)
        
        # reorder columns      
        intact_plus_biogrid_selected_features_df = intact_plus_biogrid_selected_features_df.reindex(columns=['source', 'uniprot_a', 'uniprot_b', 'pubmed_id', 
                                                          'method', 'interaction_type', 'intact_score'])
        
        logger.debug("merged intact and biogrid interactions")
                         
        # merge intact+biogrid with string
        self.all_selected_features_df = pd.merge(intact_plus_biogrid_selected_features_df, string_refined_df_selected_features, on=["uniprot_a", "uniprot_b"], how="outer")
        
        # merge source_x and source_y columns
        self.all_selected_features_df["source"] = self.all_selected_features_df[["source_x", "source_y"]].apply(lambda x: '|'.join(x.dropna()), axis=1)
        
        # drop redundant columns
        self.all_selected_features_df.drop(columns=["source_x", "source_y"], inplace=True)       
        
        # reorder columns
        self.all_selected_features_df = self.all_selected_features_df.reindex(columns=['source', 'uniprot_a', 'uniprot_b', 'pubmed_id', 
                                                          'method', 'interaction_type', 'intact_score', 
                                                          'string_combined_score', 'string_physical_combined_score'])
        
        
        # during the merging, it changes datatypes of some columns from int to float. So it needs to be reverted
        # However, for future, this may not be case so we can delete this part
        def float_to_int(element):
            if "." in str(element):                
                dot_index = str(element).index(".")
                element = str(element)[:dot_index]
                return element
            else:
                return element
        
        # first make their datatype as string
        self.all_selected_features_df["string_physical_combined_score"] = self.all_selected_features_df["string_physical_combined_score"].astype(str, errors="ignore")
        
        # then revert back them
        self.all_selected_features_df["string_physical_combined_score"] = self.all_selected_features_df["string_physical_combined_score"].apply(float_to_int)
        
        logger.debug("merged all interactions")
        t2 = time()
        logger.info(f'All data is merged and processed in {round((t2-t1) / 60, 2)} mins')
            
        if self.export_csvs:            
            all_df_path = self.export_dataframe(self.all_selected_features_df, "ppi_all")
            logger.info(f'Final STRING data is written: {all_df_path}')
                                 
        
if __name__ == "__main__":
    
    t0 = time() 
    
    args = parser.parse_args()
    
    ppi_downloader = PPI_data(output_dir=args.output_dir)
    # intact
    ppi_downloader.download_intact_data()
    ppi_downloader.intact_process()
    # biogrid
    ppi_downloader.download_biogrid_data()
    ppi_downloader.biogrid_process()
    # string
    ppi_downloader.download_string_data()
    string_ints = ppi_downloader.string_process()
    # merge all 3 databases
    check_df = ppi_downloader.merge_mall()    
   
    t1 = time()
    print(f'Done in total {round((t1-t0) / 60, 2 )} mins')
