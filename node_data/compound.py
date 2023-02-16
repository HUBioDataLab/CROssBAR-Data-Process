from __future__ import annotations

from pypath.share import curl, settings
from pypath.inputs import chembl, uniprot, unichem
from contextlib import ExitStack
from bioregistry import normalize_curie

import pandas as pd
import numpy as np

from time import time

from tqdm import tqdm

class Compound:
    """
    Class that downloads compound data using pypath and reformats it to be ready
    for import into a BioCypher database.
    """

    def download_compound_data(self, cache=False, debug=False, retries=3,):
        """
        Wrapper function to download compound data from Chembl database using pypath.
        
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
            
            t0 = time()
            
            print('Downloading compound nodes from Chembl database...')
            self.compounds = chembl.chembl_molecules()

            print('Downloading bioactivities from Chembl database...')
            self.chembl_acts = chembl.chembl_activities(standard_relation='=', )

            print('Downloading document mapping from Chembl database...')
            self.document_to_pubmed = chembl.chembl_documents()

            print('Downloading target mapping from Chembl database...')
            chembl_targets = chembl.chembl_targets()

            # filter out TrEMBL targets
            swissprots = list(uniprot._all_uniprots(organism = '*', swissprot =True))
            chembl_targets = [i for i in chembl_targets if i.accession in swissprots]
            self.target_dict = {i.target_chembl_id: i.accession for i in chembl_targets}

            print('Downloading assays from Chembl database...')
            chembl_assays = chembl.chembl_assays()
            self.assay_dict = {i.assay_chembl_id:i for i in chembl_assays if i.assay_type == 'B'}
            
            # chembl to drugbank mapping
            self.chembl_to_drugbank = unichem.unichem_mapping('chembl', 'drugbank')
            self.chembl_to_drugbank = {k:list(v)[0] for k, v in self.chembl_to_drugbank.items()}
            
            t1 = time()
            print(f'Chembl data is downloaded in {round((t1-t0) / 60, 2)} mins')            
    
    def get_compound_nodes(self):
        """
        Reformats compound node data to be ready for import into a BioCypher database.
        """
        print('Writing compound nodes...')
        
        # will filter out chembl ids whether they have edge or not
        activities_chembl = set()
        for act in self.chembl_acts:
            if act.assay_chembl in self.assay_dict and all([True if item else False for item in [act.standard_value,
                                                                                            act.standard_type,
                                                                                           self.target_dict.get(act.target_chembl, None)]]):

                activities_chembl.add(act.chembl)
                
        
        self.compound_nodes = []
        for compound in tqdm(self.compounds):
            
            if compound.chembl not in self.chembl_to_drugbank and compound.chembl in activities_chembl:
                
                compound_id = normalize_curie('chembl:' + compound.chembl)
                props = {
                    'type': compound.type,
                    'full_mwt': compound.full_mwt,
                    'species': compound.species,
                    'heavy_atoms': compound.heavy_atoms,
                    'alogp': compound.alogp,
                    'inchi': compound.std_inchi,
                    'inchikey': compound.std_inchi_key,
                }
                
                self.compound_nodes.append((compound_id, 'compound', props))
    
    
    def get_median(self, element):
        return round(float(element.dropna().median()), 3)

    def get_middle_row(self, element):
        if len(list(element.index)) == 1:
            return element.values[0]
        elif len(list(element.dropna().index)) == 0:
            return np.nan
        elif len(list(element.dropna().index)) % 2 == 1:
            middle = len(list(element.dropna().index)) // 2
            return element.dropna().values[middle]
        else:
            middle = round((len(list(element.dropna().index))/2 + 0.00001))
            return element.dropna().values[middle]

    def aggregate_column_level(self, element, joiner="|"):
        import numpy as np

        _set = set()
        for e in set(element.dropna().values):
            if joiner in e:
                for i in e.split(joiner):
                    _set.add(i)
            else:
                _set.add(e)

        if _set:
            return joiner.join(_set)
        else:
            return np.nan
        
        
    def process_chembl_cti_data(self):
        
        print("Started Chembl processing compound-target interaction data")
        t0 = time()
        
        df_list = []
        
        # filter activities
        for act in self.chembl_acts:
            if act.assay_chembl in self.assay_dict and act.chembl not in self.chembl_to_drugbank and all([True if item else False for item in [act.standard_value,
                                                                                            act.standard_type,
                                                                                           self.target_dict.get(act.target_chembl, None)]]):

                df_list.append((act.chembl, act.pchembl, act.standard_value, act.standard_type, act.assay_chembl, 
                                self.target_dict.get(act.target_chembl, None), str(self.document_to_pubmed.get(act.document, None)),
                               self.assay_dict[act.assay_chembl].confidence_score,))

        
        # create dataframe
        chembl_cti_df = pd.DataFrame(df_list, columns=["chembl", "pchembl", "activity_value", "activity_type", "assay_chembl", "uniprot_id",
                                                      "references", "confidence_score"])

        chembl_cti_df.fillna(value=np.nan, inplace=True)
        chembl_cti_df.replace("None", np.nan, inplace=True)

        # add source
        chembl_cti_df["source"] = "ChEMBL"
        
        # sort by activity value
        chembl_cti_df.sort_values(by="activity_value", ignore_index=True, inplace=True)
        
        # multiple processing
        self.chembl_cti_duplicate_removed_df = chembl_cti_df.groupby(["uniprot_id", "chembl"], sort=False, as_index=False).aggregate({
                                                                                           "chembl":"first",
                                                                                           "pchembl":self.get_median,
                                                                                           "activity_value":self.get_median,
                                                                                           "activity_type":self.get_middle_row,
                                                                                           "assay_chembl":self.aggregate_column_level,
                                                                                           "uniprot_id":"first",
                                                                                           "references":self.aggregate_column_level,
                                                                                           "confidence_score":self.get_middle_row,
                                                                                           "source":"first"}).replace("", np.nan)
        
        self.chembl_cti_duplicate_removed_df.fillna(value=np.nan, inplace=True)
        
        t1 = time()
        print(f'Chembl data is processed in {round((t1-t0) / 60, 2)} mins') 
        
        
    def get_cti_edges(self):
        """
        Reformats compound-target edge data to be ready for import into a BioCypher database.
        """

        print('Writing compound-target edges...')
        self.cti_edge_list = []
        
        for _, row in tqdm(self.chembl_cti_duplicate_removed_df.iterrows(), total=self.chembl_cti_duplicate_removed_df.shape[0]):
            
            _dict = row.to_dict()
            source = normalize_curie('chembl:' + _dict["chembl"])
            target = normalize_curie('uniprot:' +_dict["uniprot_id"])

            del _dict["chembl"] 
            del _dict["uniprot_id"]

            props = dict()
            for k, v in _dict.items():
                if str(v) != "nan":
                    if isinstance(v, str) and "|" in v:
                        props[str(k).replace(" ","_").lower()] = v.split("|")
                    else:
                        props[str(k).replace(" ","_").lower()] = v


            self.cti_edge_list.append((None, source, target, "targets", props))
