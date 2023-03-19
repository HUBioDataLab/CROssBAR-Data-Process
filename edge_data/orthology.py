from __future__ import annotations

from pypath.share import curl, settings
from pypath.inputs import oma, uniprot, pharos
from pypath.utils import taxonomy
from contextlib import ExitStack
from bioregistry import normalize_curie
from time import time

import pandas as pd
import numpy as np 

from biocypher._logger import logger
from tqdm.notebook import tqdm


OMA_ORGANISMS = {
    4932, # s. cerevisiae
    10090, # mouse
    3702,
    10116, # rat
    559292,
    9913, # cow
    1264690,
    83333,
    6239, # c. elegans
    1423,
    39947,
    44689,
    7227, # drosophila
    8355, # Xenopus laevis
    7955, # zebrafish
    9031, # chicken
    1773,
    9601,
    9598, # chimp - APES TOGETHER STRONG
    9544, # Macaca - APES TOGETHER STRONG
    9595, # GORILLA GORILLA GORILLA - APES TOGETHER STRONG
    9601, # orangutan - APES TOGETHER STRONG
}.union(set(taxonomy.taxids.keys()))


class Orthology:
    """
    Class that downloads orthology data using pypath and reformats it to be ready
    for import into a BioCypher database.
    """
    
    def download_orthology_data(self, cache=False, debug=False, retries=3,):
        """
        Wrapper function to download orthology data from various databases using pypath.

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

            self.download_oma_data()
            self.download_pharos_data()
            
            
    def process_orthology_data(self):
        
        self.process_oma_data()
        self.process_pharos_data()

    def download_oma_data(self, tax=OMA_ORGANISMS):
        """
        Downloads orthology data from OMA against human.

        Args
            tax: list of taxids to download data for.
        """

        self.entry_name_to_uniprot = uniprot.uniprot_data(field= 'entry name', reviewed = True, organism= '*')
        self.entry_name_to_uniprot = {v:k for k,v in self.entry_name_to_uniprot.items()}
        
        uniprot_to_entrez = uniprot.uniprot_data(field= 'database(GeneID)', reviewed = True, organism= '*')        
        self.uniprot_to_entrez = dict()
        for k, v in uniprot_to_entrez.items():
            self.uniprot_to_entrez[k] = v.strip(";").split(";")[0]

        logger.debug("Started downloading OMA orthology data")
        t0 = time()
        
        self.oma_orthology = []

        for t in tqdm(tax):
            tax_orthology = oma.oma_orthologs(organism_a = 9606, organism_b = t)
            tax_orthology = [i for i in tax_orthology if i.id_a in self.entry_name_to_uniprot and i.id_b in self.entry_name_to_uniprot]
            tax_orthology = [i for i in tax_orthology if self.uniprot_to_entrez.get(self.entry_name_to_uniprot[i.id_a], None) and self.uniprot_to_entrez.get(self.entry_name_to_uniprot[i.id_b], None)]
            self.oma_orthology.extend(tax_orthology)

        t1 = time()
        logger.info(f'OMA orthology data is downloaded in {round((t1-t0) / 60, 2)} mins')

    def process_oma_data(self):
        """
        Processes orthology data from OMA.
        """
        
        logger.debug("Started processing OMA orthology data")
        t0 = time()
        
        df_list = []
        for ortholog in self.oma_orthology:
            df_list.append((self.uniprot_to_entrez[self.entry_name_to_uniprot[ortholog.id_a]],
                            self.uniprot_to_entrez[self.entry_name_to_uniprot[ortholog.id_b]],
                           ortholog.rel_type, round(ortholog.score)))

        oma_orthology_df = pd.DataFrame(df_list, columns=["entrez_a", "entrez_b", "relation_type", "oma_orthology_score"])

        oma_orthology_df["source"] = "OMA"
        
        oma_orthology_df.sort_values(by="oma_orthology_score", ascending=False, inplace=True)
        
        self.oma_orthology_duplicate_removed_df = oma_orthology_df[~oma_orthology_df[["entrez_a", "entrez_b"]].apply(frozenset, axis=1).duplicated()].reset_index(drop=True)
        
        t1 = time()
        logger.info(f'OMA orthology data is processed in {round((t1-t0) / 60, 2)} mins')
        
    def download_pharos_data(self):
        """
        Downloads orthology data from Pharos.
        """

        logger.debug("Started downloading Pharos orthology data")
        t0 = time()

        uniprot_to_entrez = uniprot.uniprot_data(field= 'database(GeneID)', reviewed = True, organism= '*')
        self.entrez_to_uniprot = {}
        for k,v in uniprot_to_entrez.items():
            for entrez in v.strip(';').split(';'):
                if entrez:
                    self.entrez_to_uniprot[entrez] = k
   
        self.pharos_orthology_init = pharos.pharos_targets(orthologs=True)
    
        t1 = time()
        logger.info(f'Pharos orthology data is downloaded in {round((t1-t0) / 60, 2)} mins')
    
    def process_pharos_data(self):
        """
        Processes orthology data from Pharos.
        """
        logger.debug("Started processing Pharos orthology data")
        t0 = time()
        
        selected_species = ['Mouse', 'Cow', 'Xenopus', 'Zebrafish', 'Rat', 'C. elegans', 'S.cerevisiae', 'Chicken', 'Chimp', 'Fruitfly', 'Dog', 'Macaque', 'Pig', 'Horse',]

        df_list = []
        for protein in self.pharos_orthology_init:
            if protein["orthologs"]:
                for ortholog in protein["orthologs"]:
                    if ortholog['geneid'] and str(ortholog['geneid']) in self.entrez_to_uniprot and str(ortholog['species']) in selected_species\
                    and protein["uniprot"] in self.uniprot_to_entrez:
                        
                        df_list.append((self.uniprot_to_entrez[protein["uniprot"]], str(ortholog['geneid']) ))
                
        
        pharos_orthology_df = pd.DataFrame(df_list, columns=["entrez_a", "entrez_b"])

        pharos_orthology_df["source"] = "Pharos"
        
        self.pharos_orthology_duplicate_removed_df = pharos_orthology_df[~pharos_orthology_df[["entrez_a", "entrez_b"]].apply(frozenset, axis=1).duplicated()].reset_index(drop=True)

        t1 = time()
        logger.info(f'Pharos orthology data is processed in {round((t1-t0) / 60, 2)} mins')

    def merge_source_column(self, element, joiner="|"):
        """
        Merges source columns' entries in selected dataframes
        """
        _list = []
        for e in list(element.dropna().values):
            if joiner in e:
                for i in e.split(joiner):
                    _list.append(i)
            else:
                _list.append(e)            

        return joiner.join(list(dict.fromkeys(_list).keys()))
    
    def merge_orthology_data(self):
        """
        Merges orthology data from OMA and Pharos
        """
        logger.debug("Started merged OMA and Pharos orthology data")
        t0 = time()
        
        oma_plus_pharos_orthology_df = self.oma_orthology_duplicate_removed_df.merge(self.pharos_orthology_duplicate_removed_df, how="outer",
                                                                       on=["entrez_a", "entrez_b"])
        
        oma_plus_pharos_orthology_df["source"] = oma_plus_pharos_orthology_df[["source_x", "source_y"]].apply(self.merge_source_column, axis=1)
        
        oma_plus_pharos_orthology_df.drop(columns=["source_x", "source_y"], inplace=True)
        
        self.all_orthology_df = oma_plus_pharos_orthology_df
        
        t1 = time()
        logger.info(f'OMA and Pharos orthology data is merged in {round((t1-t0) / 60, 2)} mins')        
    
    def get_orthology_edges(self, early_stopping=500):
        """
        Reformats orthology data to be ready for import into a BioCypher database.
        """

        self.edge_list = []
        
        logger.info("Preparing orthology edges.")
        
        for index, row in tqdm(self.all_orthology_df.iterrows(), total=self.all_orthology_df.shape[0]):
            _dict = row.to_dict()
            
            source = normalize_curie('ncbigene:' + _dict["entrez_a"])
            target = normalize_curie('ncbigene:' + _dict["entrez_b"])
            
            del _dict["entrez_a"], _dict["entrez_b"]
            
            props = dict()
            for k, v in _dict.items():
                if str(v) != "nan":
                    if isinstance(v, str) and "|" in v:
                        props[str(k).replace(" ","_").lower()] = v.replace("'", "^").split("|")
                    else:
                        props[str(k).replace(" ","_").lower()] = str(v).replace("'", "^")


            self.edge_list.append((None, source, target, "gene_is_orthologous_with_gene", props))
            
            if early_stopping and (index+1) == early_stopping:
                break
