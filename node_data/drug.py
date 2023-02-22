
from __future__ import annotations

from pypath.share import curl, settings, common
from pypath.inputs import drugbank, drugcentral, stitch, string, uniprot, dgidb, pharos, ctdbase, unichem, chembl, uniprot, ddinter
import kegg_local
from contextlib import ExitStack
from typing import Literal
from bioregistry import normalize_curie
from tqdm import tqdm
from time import time
import collections

import pandas as pd
#import modin.pandas as pd
import numpy as np


class Drug:
    """
    WARNING: Drugbank, Drug Central and STITCH database download urls contain
        version number/date. Please update this urls (or request an update) in 
        resources.urls module of the pypath library before using this script
        in order to access the latest versions of data.

    Class that downloads drug data using pypath and reformats it to be ready
    for import into a BioCypher database.
    """

    def __init__(self, drugbank_user, drugbank_passwd):
        """
        Args
            drugbank_user: drugbank username
            drugbank_passwd: drugbank password
        """

        self.edge_list = []
        self.user = drugbank_user
        self.passwd = drugbank_passwd
        
        self.swissprots = list(uniprot._all_uniprots(organism = '*', swissprot=True))



    def download_drug_data(
        self, cache=False, debug=False, retries=6, 
        ):

        """
        Wrapper function to download drug data from various databases using pypath.

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
            
            # DTI
            self.download_chembl_dti_data()
            self.download_drugbank_node_data(self.user, self.passwd)
            self.download_drugbank_dti_data()
            self.download_pharos_dti_data()
            self.download_dgidb_dti_data()
            self.download_kegg_dti_data()
            self.download_stitch_dti_data()
            
            # DDI
            self.download_kegg_ddi_data()
            self.download_ddinter_ddi_data()
            
            # DGI
            self.download_ctd_data()
            
    
    def process_drug_data(self):
        
        # DTI
        self.process_drugbank_node_data()
        self.process_drugbank_dti_data()
        self.process_chembl_dti_data()
        self.process_pharos_dti_data()
        self.process_dgidb_dti_data()
        self.process_stitch_dti_data()
        self.process_kegg_dti_data()
        
        # DDI
        self.process_kegg_ddi_data()
        self.process_ddinter_ddi_data()
        
        # DGI
        self.process_ctd_data()
        

    def download_drugbank_node_data(self, user, passwd, drugbank_node_fields=['InChI', 'InChIKey', 'cas_number', 'name', 'groups', 
                                                                    'general_references', 'atc_codes', 'zinc', 'chembl', 
                                                                    'bindingdb', 'clinicaltrials', 'chebi', 'pubchem',
                                                                   'KEGG Drug', 'RxCUI', 'PharmGKB', 'PDB', 'Drugcentral']):

        """
        Wrapper function to download DrugBank drug entries, DTI and DDI data using pypath
        Args
            user: drugbank username
            passwd: drugbank password
        """
        
        fields_list = ['cas_number', 'name', 'groups', 'general_references', 'atc_codes',]
        unichem_external_fields_list = ['zinc', 'chembl', 'bindingdb', 'clinicaltrials', 'chebi', 'pubchem']
        drugbank_external_fields_list = ['KEGG Drug', 'RxCUI', 'PharmGKB', 'PDB', 'Drugcentral']
        
        # 
        unichem_external_fields = []
        drugbank_external_fields = []
        fields = []
        self.add_inchi = False
        self.add_inchikey = False
        for f in drugbank_node_fields:
            if f == 'InChI':
                self.add_inchi = True
            elif f == 'InChIKey':
                self.add_inchikey = True
            elif f in fields_list:
                fields.append(f)
            elif f in unichem_external_fields_list:
                unichem_external_fields.append(f)
            elif f in drugbank_external_fields_list:
                drugbank_external_fields.append(f)
            else:
                raise ValueError(f" {f} is an inappropriate field name. Please provide a sutiable field name")            
        
        self.unichem_external_fields = unichem_external_fields
        self.drugbank_external_fields = drugbank_external_fields
        
        print('Downloading Drugbank drug node data')
        t0 = time()
        # Drugbank Object
        self.drugbank_data = drugbank.DrugbankFull(user = user, passwd = passwd)

        # external ids
        self.drugbank_drugs_external_ids = self.drugbank_data.drugbank_external_ids_full()
        
        # inchi and inchikey
        if self.add_inchi or self.add_inchikey:
            self.drugbank_properties = self.drugbank_data.drugbank_properties_full()
        
        # core properties
        self.drugbank_drugs_detailed = self.drugbank_data.drugbank_drugs_full(fields = fields)
        
        # combine all external id mappings to self.drug_mappings_dict
        self.get_external_database_mappings(unichem_external_fields=unichem_external_fields, 
                                           drugbank_external_fields=drugbank_external_fields)
        
        t1 = time()
        print(f'Drugbank drug node data is downloaded in {round((t1-t0) / 60, 2)} mins')
        
    def process_drugbank_node_data(self):
        
        print('Processing Drugbank drug node data')
        t0 = time()
        
        self.drugbank_drugs = {}
        all_fields = list(self.drugbank_drugs_detailed[0]._fields) + self.drugbank_external_fields + self.unichem_external_fields
        
        for drug in self.drugbank_drugs_detailed:
    
            drugbank_id = drug.drugbank_id
            temp_dict = (
                drug._asdict() | self.drug_mappings_dict[drugbank_id]
                    if drugbank_id in self.drug_mappings_dict else
                    drug._asdict()
                    )

            self.drugbank_drugs[drugbank_id] = {f: temp_dict.get(f, None) for f in all_fields}
            
            if self.add_inchi:
                self.drugbank_drugs[drugbank_id]["InChI"] = self.drugbank_properties.get(drugbank_id, {}).get("InChI", None)
                
            if self.add_inchikey:
                self.drugbank_drugs[drugbank_id]["InChIKey"] = self.drugbank_properties.get(drugbank_id, {}).get("InChIKey", None)

            del self.drugbank_drugs[drugbank_id]['drugbank_id']
            
        
        t1 = time()
        print(f'Drugbank drug node data is processed in {round((t1-t0) / 60, 2)} mins')

    def download_drugbank_dti_data(self):
        # edge data
        print('Downloading Drugbank DTI data')
        t0 = time()
        self.drugbank_dti = self.drugbank_data.drugbank_targets_full(fields=['drugbank_id', 'actions', 'references', 'known_action', 'polypeptide',])
        t1 = time()
        print(f'Drugbank DTI data is downloaded in {round((t1-t0) / 60, 2)} mins')

    def get_external_database_mappings(self, 
                                       unichem_external_fields: list | None = None, 
                                       drugbank_external_fields: list | None = None,):
        
        if not unichem_external_fields:
            unichem_external_fields = ['zinc', 'chembl', 'bindingdb', 'clinicaltrials', 'chebi', 'pubchem']
        
        if not drugbank_external_fields:
            drugbank_external_fields = ['KEGG Drug', 'RxCUI', 'PharmGKB', 'PDB', 'Drugcentral']
        
        
        # create dictionaries for every unichem external fields
        unichem_drugbank_to_zinc_mapping = unichem.unichem_mapping('drugbank', 'zinc')
        unichem_drugbank_to_chembl_mapping = unichem.unichem_mapping('chembl', 'drugbank')
        unichem_drugbank_to_chembl_mapping = {list(v)[0]:k for k, v in unichem_drugbank_to_chembl_mapping.items()}
        unichem_drugbank_to_bindingdb_mapping = unichem.unichem_mapping('drugbank', 'bindingdb')
        unichem_drugbank_to_clinicaltrials_mapping = unichem.unichem_mapping('drugbank', 'clinicaltrials')
        unichem_drugbank_to_chebi_mapping = unichem.unichem_mapping('drugbank', 'chebi')
        unichem_drugbank_to_pubchem_mapping = unichem.unichem_mapping('drugbank', 'pubchem')

        # store above dicts in a unified dict
        self.unichem_external_fields_dict = {'zinc':unichem_drugbank_to_zinc_mapping, 'chembl':unichem_drugbank_to_chembl_mapping,
                                'bindingdb':unichem_drugbank_to_bindingdb_mapping, 'clinicaltrials':unichem_drugbank_to_clinicaltrials_mapping,
                                'chebi':unichem_drugbank_to_chebi_mapping, 'pubchem':unichem_drugbank_to_pubchem_mapping}
        
        # arrange unichem dict for selected unichem fields
        for field in self.unichem_external_fields_dict.keys():
            if field not in unichem_external_fields:
                del self.unichem_external_fields_dict[field]
                
                
        unichem_drugs_id_mappings = collections.defaultdict(dict)
        for k in self.drugbank_drugs_external_ids.keys():
            for field_name, field_dict in self.unichem_external_fields_dict.items():
                mapping = field_dict.get(k, None)
                if mapping and field_name != "chembl":
                    mapping = list(mapping)[0]

                unichem_drugs_id_mappings[k][field_name] = mapping

        # get drugcentral mappings
        self.cas_to_drugbank = {drug.cas_number:drug.drugbank_id for drug in self.drugbank_drugs_detailed if drug.cas_number}
        drugcentral_to_cas = drugcentral.drugcentral_mapping(id_type='drugcentral', target_id_type='cas')
        
        # create drugbank-drugcentral, drugcentral-drugbank, chembl-drugbank mappings that will be used for the future processes
        chembl_to_drugbank = unichem.unichem_mapping('chembl', 'drugbank')
        self.chembl_to_drugbank = {k:list(v)[0] for k, v in chembl_to_drugbank.items()}
        
        # create kegg-drugbank mapping
        self.kegg_to_drugbank = {v['KEGG Drug']:k for k, v in self.drugbank_drugs_external_ids.items() if v.get('KEGG Drug', None)}
        
        self.drugcentral_to_drugbank = collections.defaultdict(list)
        self.drugbank_to_drugcentral = collections.defaultdict(None)
        for k, v in drugcentral_to_cas.items():
            if list(v)[0] in self.cas_to_drugbank and k:
                self.drugbank_to_drugcentral[self.cas_to_drugbank[list(v)[0]]] = k
                self.drugcentral_to_drugbank[k] = self.cas_to_drugbank[list(v)[0]]
                
                # add drugcentral id to drugbank_drugs_external_ids
                if self.cas_to_drugbank[list(v)[0]] in self.drugbank_drugs_external_ids:
                    self.drugbank_drugs_external_ids[self.cas_to_drugbank[list(v)[0]]]['Drugcentral'] = k
            
        # create final external id mapping dict
        self.drug_mappings_dict = collections.defaultdict(dict)
        for k in self.drugbank_drugs_external_ids.keys():
            drugbank_mappings = {field:self.drugbank_drugs_external_ids[k].get(field, None) for field in drugbank_external_fields}
            unichem_mappings = unichem_drugs_id_mappings[k]
            self.drug_mappings_dict[k] = (drugbank_mappings | unichem_mappings)
        

    def process_drugbank_dti_data(self):
        
        
        def get_uniprot_ids(element):
            if element:
                if isinstance(element, tuple) and element[1] == 'Swiss-Prot':
                    return [element[0]]
                elif isinstance(element, list):
                    _list = []
                    for x in element:
                        if x[1] == 'Swiss-Prot':
                            _list.append(x[0])

                    return _list
            else:
                return None

        def aggregate_one_field(element, joiner="|"):
            if element:
                return joiner.join(set(element))
            else:
                return None
        
        print('Processing Drugbank DTI data')
        t0 = time()
        
        # create list instance for pandas dataframes 
        df_list = []

        selected_fields = ["actions", "references", "known_action"]
        
        # process drugbank dti
        for dti in self.drugbank_dti:
            uniprot_ids = get_uniprot_ids(dti.polypeptide)

            if not uniprot_ids or not dti.drugbank_id:
                continue

            for _id in uniprot_ids:
                attributes = dti._asdict()

                _list = []

                _list.append(dti.drugbank_id)
                _list.append(_id)

                for field in selected_fields:

                    if field == "references":
                        if isinstance(attributes[field], list):
                            aggregated = aggregate_one_field([i for i in attributes[field] if i is not None])
                            _list.append(aggregated)
                        else:
                            _list.append(attributes[field])

                    elif field == "actions":
                        if isinstance(attributes[field], list):
                            _list.append(attributes[field][0])
                        else:
                            _list.append(attributes[field])

                    else:
                        _list.append(attributes[field])


                df_list.append(_list)


        # create pandas dataframe
        drugbank_dti_df = pd.DataFrame(df_list, columns=["drugbank_id", "uniprot_id", "mechanism_of_action_type", "references", "known_action"])
        drugbank_dti_df.fillna(value=np.nan, inplace=True)
        drugbank_dti_df.replace("", np.nan, inplace=True)

        # add source
        drugbank_dti_df["source"] = "DrugBank"
        
        # Aggregate references field and remove duplicates
        self.drugbank_dti_duplicate_removed_df = drugbank_dti_df.groupby(["drugbank_id", "uniprot_id"], sort=False, as_index=False).aggregate({"drugbank_id":"first",
                                                                                              "uniprot_id":"first",
                                                                                              "mechanism_of_action_type":"first",
                                                                                              "references":self.aggregate_column_level,
                                                                                              "known_action":"first",
                                                                                              "source":"first"}).replace("", np.nan)

        self.drugbank_dti_duplicate_removed_df.fillna(value=np.nan, inplace=True)

        t1 = time()
        print(f'Drugbank DTI data is processed in {round((t1-t0) / 60, 2)} mins')
        
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
        
    def merge_source_column(self, element, joiner="|"):
        _list = []
        for e in list(element.dropna().values):
            if joiner in e:
                for i in e.split(joiner):
                    _list.append(i)
            else:
                _list.append(e)

        return joiner.join(list(dict.fromkeys(_list).keys()))


    def download_dgidb_dti_data(self):

        """
        Wrapper function to download DGIdb DTI data using pypath

        It returns ChEMBL IDs as drug IDs and Entrez Gene IDs as protein IDs.
        """

        # edge data
        print('Downloading DGIdb DTI data')
        t0 = time()
        
        self.dgidb_dti = dgidb.dgidb_interactions()
        
        # map entrez gene ids to swissprot ids
        uniprot_to_entrez = uniprot.uniprot_data("database(GeneID)", "*", True)
        self.entrez_to_uniprot = collections.defaultdict(list)
        for k,v in uniprot_to_entrez.items():
            for entrez_id in list(filter(None, v.split(";"))):
                self.entrez_to_uniprot[entrez_id].append(k)
                
        t1 = time()
        print(f'Dgidb DTI data is downloaded in {round((t1-t0) / 60, 2)} mins')
                
    def process_dgidb_dti_data(self):
        
        print('Processing DGIdb DTI data')
        t0 = time()

        df_list = []

        for dti in self.dgidb_dti:
            if dti.entrez and dti.drug_chembl and self.entrez_to_uniprot[dti.entrez]:
                if dti.pmid:
                    pmid = "|".join(dti.pmid.split(","))
                else:
                    pmid = None


                df_list.append((self.entrez_to_uniprot[dti.entrez][0], 
                                self.chembl_to_drugbank.get(dti.drug_chembl.split(":")[1], None),
                                dti.type, dti.score, pmid))


        dgidb_dti_df = pd.DataFrame(df_list, columns=["uniprot_id", "drugbank_id", "mechanism_of_action_type", "dgidb_score", "references"])

        dgidb_dti_df.fillna(value=np.nan, inplace=True)

        # add source
        dgidb_dti_df["source"] = "DGIdb"
        
        # sort by dgidb_score
        dgidb_dti_df.sort_values(by="dgidb_score", ignore_index=True, inplace=True, ascending=False)
        
        # remove pairs without drugbank ids
        dgidb_dti_duplicate_removed_df = dgidb_dti_df.dropna(subset="drugbank_id", axis=0).reset_index(drop=True)
        
        # remove duplicates
        self.dgidb_dti_duplicate_removed_df = dgidb_dti_duplicate_removed_df.groupby(["drugbank_id", "uniprot_id"], 
                                                                        sort=False, as_index=False).aggregate({
                                                                                                    "uniprot_id":"first",
                                                                                                    "drugbank_id":"first",
                                                                                                    "mechanism_of_action_type":self.get_middle_row,
                                                                                                    "dgidb_score":"first",
                                                                                                    "references":self.aggregate_column_level,
                                                                                                    "source":"first"}).replace("", np.nan)

        self.dgidb_dti_duplicate_removed_df.fillna(value=np.nan, inplace=True)

        t1 = time()
        print(f'Dgidb DTI data is processed in {round((t1-t0) / 60, 2)} mins')
        
    def download_kegg_dti_data(
        self,
        organism: str | list = "hsa"
        ):

        """
        Wrapper function to download KEGG DTI and DDI data using pypath

        Args:
            organism: KEGG organism code. Default is "hsa" for human.
                    If None, it downloads DTI data for all organisms in KEGG.

        """
        
        # DTI
        if organism is None:
            organism = kegg_local._kegg_list('organism')
        
        organism = common.to_list(organism)

        print(f'Downloading KEGG DTI data for {len(organism)} organism(s)')
        t0 = time()

        self.kegg_dti = set()

        for org in tqdm(organism):
            organism_dti = kegg_local.drug_to_gene(org = org)
            for k, v in organism_dti.items():
                if self.kegg_to_drugbank.get(k, None) and not isinstance(v, str):
                    drugbank_id = self.kegg_to_drugbank[k]
                    for gene_entry in v.GeneEntries:
                        if isinstance(gene_entry.uniprot_ids, (tuple, list, set)):
                            for uniprot_id in list(gene_entry.uniprot_ids):
                                if uniprot_id in self.swissprots:
                                    self.kegg_dti.add((drugbank_id, uniprot_id))
                        else:
                            if str(gene_entry.uniprot_ids) in self.swissprots:
                                self.kegg_dti.add((drugbank_id, gene_entry.uniprot_ids))
                                
        
        self.kegg_dti = list(self.kegg_dti)
        
        t1 = time()
        print(f'KEGG DTI data is downloaded in {round((t1-t0) / 60, 2)} mins')
        
        
    def process_kegg_dti_data(self):
        
        print(f'Processing KEGG DTI data')
        t0 = time()

        self.kegg_dti_df = pd.DataFrame(self.kegg_dti, columns=["drugbank_id", "uniprot_id"])

        self.kegg_dti_df["source"] = "Kegg"
        
        t1 = time()
        print(f'KEGG DTI data is processed in {round((t1-t0) / 60, 2)} mins')
        
    def download_kegg_ddi_data(self, from_csv=True):
        # DDI
        print('Downloading KEGG DDI data, this may take around 12 hours')
        t0 = time()
        if from_csv:
            print("Skipping to processing part")
        else:
            self.kegg_ddi_data = kegg_local.drug_to_drug()
        t1 = time()
        print(f'KEGG DDI data is downloaded in {round((t1-t0) / 60, 2)} mins')
        
    def process_kegg_ddi_data(self, from_csv=True):
        
        print('Processing KEGG DDI data')
        t0 = time()
        
        if from_csv:
            # even with cache it takes a lot of time to download kegg ddi data. To interfere that, processed
            # version of it saved in a csv file
            self.kegg_ddi_duplicate_removed_df = pd.read_csv('kegg_ddi_duplicate_removed_df.csv')
        else:
            
            kegg_ddi = set()
            for k, v in self.kegg_ddi_data.items():
                if self.kegg_to_drugbank.get(k, None) and not isinstance(v, str):
                    drug1_drugbank_id = self.kegg_to_drugbank[k]
                    if v.interactions:
                        for interaction in list(v.interactions):
                            if interaction.type == "drug" and interaction.id and self.kegg_to_drugbank.get(interaction.id, None):                    
                                drug2_drugbank_id = self.kegg_to_drugbank[interaction.id]
                                                                
                                if interaction.contraindication:                        
                                    contraindication = "contraindication"
                                else:
                                    contraindication = ""

                                if interaction.precaution:                        
                                    precaution = "precaution"
                                else:
                                    precaution = ""

                                if contraindication and precaution:
                                    interaction_type = "|".join([contraindication, precaution])
                                else:                        
                                    interaction_type = contraindication or precaution


                                kegg_ddi.add((drug1_drugbank_id, drug2_drugbank_id, interaction_type))


            kegg_ddi_df = pd.DataFrame(list(kegg_ddi), columns=["drug1", "drug2", "interaction_type"])
            
            kegg_ddi_df.replace("", np.nan, inplace=True)

            kegg_ddi_df["source"] = "Kegg"

            self.kegg_ddi_duplicate_removed_df = kegg_ddi_df[~kegg_ddi_df[["drug1", "drug2"]].apply(frozenset, axis=1).duplicated()].reset_index(drop=True)

        t1 = time()
        print(f'KEGG DDI data is processed in {round((t1-t0) / 60, 2)} mins')
            
            
    def download_ddinter_ddi_data(self):
        
        print('Downloading DDInter DDI data')
        t0 = time()
        
        ddinter_mappings = ddinter.get_all_mappings()
        
        self.ddinter_to_drugbank = {mapping.ddinter_id:mapping.drugbank for mapping in ddinter_mappings if mapping.drugbank}
        
        self.ddinter_interactions = ddinter.get_all_interactions()
        
        t1 = time()
        print(f'DDInter DDI data is downloaded in {round((t1-t0) / 60, 2)} mins')
        
    def process_ddinter_ddi_data(self):
        
        print('Processing DDInter DDI data')
        t0 = time()
        
        df_list = []
        for interaction in self.ddinter_interactions:
            if self.ddinter_to_drugbank.get(interaction.drug1_id, None) and self.ddinter_to_drugbank.get(interaction.drug2_id, None):
                if isinstance(interaction.level, tuple):
                    if len(interaction.level) > 1:
                        level = "|".join(list(interaction.level))
                    else:
                        level = list(interaction.level)[0]
                else:
                    level = interaction.level

                if isinstance(interaction.actions, tuple):
                    if len(interaction.actions) > 1:
                        actions = "|".join(list(interaction.actions))
                    else:
                        actions = list(interaction.actions)[0]
                else:
                    actions = interaction.actions

                df_list.append((self.ddinter_to_drugbank[interaction.drug1_id], self.ddinter_to_drugbank[interaction.drug2_id],
                               level, actions))
                
        ddinter_df = pd.DataFrame(df_list, columns=["drug1", "drug2", "interaction_level", "interaction_type"])

        ddinter_df["source"] = "DDInter"
        
        self.ddinter_duplicate_removed_df = ddinter_df[~ddinter_df[["drug1", "drug2"]].apply(frozenset, axis=1).duplicated()].reset_index(drop=True)
        
        t1 = time()
        print(f'DDInter DDI data is processed in {round((t1-t0) / 60, 2)} mins')
        
    def download_pharos_dti_data(self):

        print('Downloading Pharos DTI data')
        t0 = time()
        
        self.pharos_dti = pharos.pharos_targets(ligands=True)
        
        t1 = time()
        print(f'Pharos DTI data is downloaded in {round((t1-t0) / 60, 2)} mins')
    
    def process_pharos_dti_data(self):
        
        print('Processing Pharos DTI data')
        t0 = time()

        df_list = []
        for dti in self.pharos_dti:

            if dti["ligands"]:
                for lig in dti["ligands"]:
                    synonyms = dict([(syn["name"], syn["value"].split(",")[0]) for syn in lig["synonyms"] if syn["name"] in ["ChEMBL", "DrugCentral"]])

                    for act in lig["activities"]:
                        # add dtis that have activity value and type
                        if act["value"] and act["type"] and act["type"] != "-":

                            if act["pubs"]:
                                _list = []
                                for pub in act["pubs"]:                        
                                    if pub["__typename"] == "PubMed" and pub["pmid"]:
                                        _list.append(pub["pmid"])

                                pubmeds = "|".join(set(_list))
                            else:
                                pubmeds = None


                            df_list.append((dti["uniprot"], act["type"], act["moa"], act["value"], 
                                            pubmeds, self.drugcentral_to_drugbank.get(synonyms.get("DrugCentral", None), None)))

        pharos_dti_df = pd.DataFrame(df_list, columns=["uniprot_id", "activity_type", "mechanism_of_action_type", 
                                                       "activity_value", "references", "drugbank_id"])

        pharos_dti_df.fillna(value=np.nan, inplace=True)
        pharos_dti_df.replace("", np.nan, inplace=True)

        # add source
        pharos_dti_df["source"] = "Pharos"
        
        
        # Remove rows without drugbank id 
        pharos_dti_df.dropna(axis=0, subset="drugbank_id", inplace=True)
        
        # Sort by activity_value
        pharos_dti_df.sort_values(by="activity_value", ignore_index=True, inplace=True)
        
        # For every drug-target pair instance the preprocess as follows:
        # - get middle row for activity_type and mechanism_of_action_type
        # - get median of activity_value
        # - aggregate all the references
        self.pharos_dti_duplicate_removed_df = pharos_dti_df.groupby(["uniprot_id", "drugbank_id"], sort=False, as_index=False).aggregate({"uniprot_id":"first",
                                                                                              "activity_type":self.get_middle_row,
                                                                                              "mechanism_of_action_type":self.get_middle_row,
                                                                                              "activity_value":self.get_median,
                                                                                        "references":self.aggregate_column_level,
                                                                                        "drugbank_id":"first",
                                                                                        "source":"first"}).replace("", np.nan)

        self.pharos_dti_duplicate_removed_df.fillna(value=np.nan, inplace=True)
        
        t1 = time()
        print(f'Pharos DTI data is processed in {round((t1-t0) / 60, 2)} mins')
        
    def download_chembl_dti_data(self):
        
        print('Downloading Chembl DTI data')
        t0 = time()
                
        self.chembl_acts = chembl.chembl_activities(standard_relation='=')
        self.chembl_document_to_pubmed = chembl.chembl_documents()
        self.chembl_targets = chembl.chembl_targets()
        self.chembl_assays = chembl.chembl_assays()
        self.chembl_mechanisms = chembl.chembl_mechanisms()
        
        t1 = time()
        print(f'Chembl DTI data is downloaded in {round((t1-t0) / 60, 2)} mins')
        
    def process_chembl_dti_data(self):
        
        print('Processing Chembl DTI data')
        t0 = time()
        
        mechanism_dict = {i.chembl:i._asdict() for i in self.chembl_mechanisms}
        targets = [i for i in self.chembl_targets if i.accession in self.swissprots]
        target_dict = {i.target_chembl_id:i.accession for i in targets}
        assay_dict = {i.assay_chembl_id:i for i in self.chembl_assays if i.assay_type == 'B'}
        
        df_list = []

        for act in self.chembl_acts:
            # if activity is belong to a binding assay and if it have activity_type, activity_value and target uniprot id
            if act.assay_chembl in assay_dict and all([True if item else False for item in [act.standard_value,
                                                                                            act.standard_type,
                                                                                           target_dict.get(act.target_chembl, None)]]):

                df_list.append((act.pchembl, act.standard_value, act.standard_type,
                                target_dict.get(act.target_chembl, None), str(self.chembl_document_to_pubmed.get(act.document, None)),
                               assay_dict[act.assay_chembl].confidence_score, self.chembl_to_drugbank.get(act.chembl, None),
                               mechanism_dict.get(act.chembl, {}).get("action_type", None), 
                                mechanism_dict.get(act.chembl, {}).get("direct_interaction", None),
                               mechanism_dict.get(act.chembl, {}).get("disease_efficacy", None),
                               mechanism_dict.get(act.chembl, {}).get("mechanism_of_action", None),))

        
        # create pandas dataframe
        chembl_cti_df = pd.DataFrame(df_list, columns=["pchembl", "activity_value", "activity_type", "uniprot_id",
                                                      "references", "confidence_score", "drugbank_id", "mechanism_of_action_type", "direct_interaction",
                                                      "disease_efficacy", "mechanism_of_action"])

        chembl_cti_df.fillna(value=np.nan, inplace=True)
        chembl_cti_df.replace("None", np.nan, inplace=True)

        # add source
        chembl_cti_df["source"] = "ChEMBL"
        
        # SORT BY activity_value
        chembl_cti_df.sort_values(by="activity_value", ignore_index=True, inplace=True)
        
        chembl_dti_df = chembl_cti_df.dropna(subset=["drugbank_id"], axis=0).reset_index(drop=True)

        self.chembl_dti_duplicate_removed_df = chembl_dti_df.groupby(["uniprot_id", "drugbank_id"], sort=False, as_index=False).aggregate({
                                                                                                   "pchembl":self.get_median,
                                                                                                   "activity_value":self.get_median,
                                                                                                   "activity_type":self.get_middle_row,
                                                                                                   "uniprot_id":"first",
                                                                                                   "references":self.aggregate_column_level,
                                                                                                   "confidence_score":self.get_middle_row,
                                                                                                   "drugbank_id":"first",
                                                                                                   "mechanism_of_action_type":self.get_middle_row,
                                                                                                   "direct_interaction":self.get_middle_row,
                                                                                                   "disease_efficacy":self.get_middle_row,
                                                                                                   "mechanism_of_action":self.get_middle_row,
                                                                                                   "source":"first"}).replace("", np.nan)

        self.chembl_dti_duplicate_removed_df.fillna(value=np.nan, inplace=True)
        
        t1 = time()
        print(f'Chembl DTI data is processed in {round((t1-t0) / 60, 2)} mins')
        
        
    def download_ctd_data(self):
        
        print('Downloading CTD DGI data')
        t0 = time()
        
        self.ctd_dgi = ctdbase.ctdbase_relations(relation_type='chemical_gene')
        
        t1 = time()
        print(f'CTD DGI data is downloaded in {round((t1-t0) / 60, 2)} mins')
        
    def process_ctd_data(self):
        
        print('Processing CTD DGI data')
        t0 = time()
        
        df_list = []
        
        for cgi in self.ctd_dgi:
            if cgi.GeneID and cgi.CasRN and any([True if item in [['increases', 'expression'], ['decreases', 'expression']] else False for item in cgi.InteractionActions])\
            and self.cas_to_drugbank.get(cgi.CasRN, None):

                # if both (increases and decreases expression) of them occur in same InteractionActions field, don't add to list
                if len([
                        item for item in cgi.InteractionActions if item in
                    [['increases', 'expression'], ['decreases', 'expression']]
                ]) > 1:
                    continue

                if isinstance(cgi.PubMedIDs, list):
                    pmid = "|".join(cgi.PubMedIDs)
                else:
                    pmid = cgi.PubMedIDs

                interaction_actions = "_".join([
                    item for item in cgi.InteractionActions if item in
                    [['increases', 'expression'], ['decreases', 'expression']]
                ][0])
                
                
                df_list.append(
                    (cgi.GeneID,
                     self.cas_to_drugbank.get(cgi.CasRN), 
                     interaction_actions, 
                     pmid))
                
                
        ctd_cgi_df = pd.DataFrame(df_list, columns=["entrez_id", "drugbank_id", "action_type", "references"])
        
        
        def detect_conflicting_action_type(element):
            # if decreases_expression and increases_expression occur in same drug-gene pair, the pair is probably a bad entry
            if len(set(element.dropna().values)) > 1:
                return np.nan
            else:
                return list(set(element.dropna().values))[0]
            
            
        self.ctd_cgi_duplicate_removed_df = ctd_cgi_df.groupby(["drugbank_id", "entrez_id"], sort=False, as_index=False).aggregate({"entrez_id":"first",
                                                                                        "drugbank_id":"first",
                                                                                        "action_type":detect_conflicting_action_type,
                                                                                        "references":"first"}).replace("", np.nan)
        
        self.ctd_cgi_duplicate_removed_df.dropna(subset="action_type", inplace=True)
        
        t1 = time()
        print(f'CTD DGI data is processed in {round((t1-t0) / 60, 2)} mins')
        

    def download_stitch_dti_data(
            self, 
            organism: str | list = "9606", 
            score_threshold: int | Literal[
                'highest_confidence',
                'high_confidence',
                'medium_confidence',
                'low_confidence',
                ] = 'high_confidence', 
            physical_interaction_score: bool = False, # currently this arg doesnt work for organisms other than human. can be fixed if necessary.
            ):
        """
        Wrapper function to download STITCH DTI data using pypath

        Args:
            organism: Name or NCBI Taxonomy IDs of organisms
                If None, DTI for all organisms will be downloaded.
            score_threshold: Minimum required interaction score. user can use
                pre-defined confidence limits or can define a custom value.
            physical_interaction_score: If True, returns physical interaction scores of interactions.

        """

        if organism is None:
            organism = string.string_species()

        organism = common.to_list(organism)
        
        print("Started downloading STITCH data")
        t0 = time()

        # map string ids to swissprot ids
        uniprot_to_string = uniprot.uniprot_data("database(STRING)", "*", True)
        self.string_to_uniprot = collections.defaultdict(list)
        for k,v in uniprot_to_string.items():
            for string_id in list(filter(None, v.split(";"))):
                self.string_to_uniprot[string_id.split(".")[1]].append(k)
        
        # mapping to convert pubchem ids to drugbank ids
        if self.unichem_external_fields_dict.get("pubchem", None):
            
            self.pubchem_to_drugbank = dict()
            for k, v in self.unichem_external_fields_dict["pubchem"].items():
                if len(v) > 1:
                    for value in list(v):
                        self.pubchem_to_drugbank[value] = k
                else:
                    self.pubchem_to_drugbank[list(v)[0]] = k
        else:
            
            self.pubchem_to_drugbank = dict()
            for k, v in unichem.unichem_mapping('drugbank', 'pubchem').items():
                if len(v) > 1:
                    for value in list(v):
                        self.pubchem_to_drugbank[value] = k
                else:
                    self.pubchem_to_drugbank[list(v)[0]] = k
        
        
        self.stitch_ints = []

        for tax in tqdm(organism):
            try:
                organism_stitch_ints = [
                    i for i in stitch.stitch_links_interactions(ncbi_tax_id=int(tax), score_threshold=score_threshold, physical_interaction_score= physical_interaction_score)
                    if i.partner_b in self.string_to_uniprot and i.partner_a in self.pubchem_to_drugbank] # filter with swissprot ids

                # TODO: later engage below print line to biocypher log 
                # print(f"Downloaded STITCH data with taxonomy id {str(tax)}")

                if organism_stitch_ints:
                    self.stitch_ints.extend(organism_stitch_ints)
            
            except TypeError: #'NoneType' object is not an iterator
                pass
                # TODO: later engage below print line to biocypher log 
                # print(f'Skipped tax id {tax}. This is most likely due to the empty file in database. Check the database file.')

        t1 = time()        
        print(f'STITCH data is downloaded in {round((t1-t0) / 60, 2)} mins')
        
    def process_stitch_dti_data(self):
        
        print("Started processing STITCH data")
        t0 = time()
        
        df_list = []

        for dti in self.stitch_ints:
            df_list.append((self.pubchem_to_drugbank[dti.partner_a], self.string_to_uniprot[dti.partner_b][0], dti.combined_score))


        stitch_dti_df = pd.DataFrame(df_list, columns=["drugbank_id", "uniprot_id", "stitch_combined_score"])

        stitch_dti_df.fillna(value=np.nan, inplace=True)

        # add source
        stitch_dti_df["source"] = "STITCH"
        
        # sort by stitch_combined_score
        stitch_dti_df.sort_values(by="stitch_combined_score", ignore_index=True, inplace=True, ascending=False)
        
        # remove duplicates
        self.stitch_dti_duplicate_removed_df = stitch_dti_df.groupby(["drugbank_id", "uniprot_id"], sort=False, as_index=False).aggregate({
                                                                                           "drugbank_id":"first",
                                                                                           "uniprot_id":"first",
                                                                                           "stitch_combined_score":self.get_median, 
                                                                                           "source":"first"}).replace("", np.nan)

        self.stitch_dti_duplicate_removed_df.fillna(value=np.nan, inplace=True)
        
        t1 = time()        
        print(f'STITCH data is processed in {round((t1-t0) / 60, 2)} mins')
        
    def merge_all_dtis(self):
        
        print("Started merging Drugbank and Chembl DTI data")
        t0 = time()

        # merge drugbank and chembl dti
        drugbank_plus_chembl_dti_df = self.drugbank_dti_duplicate_removed_df.merge(self.chembl_dti_duplicate_removed_df, how="outer", on=["uniprot_id", "drugbank_id"])
        
        # merge references
        drugbank_plus_chembl_dti_df["references"] = drugbank_plus_chembl_dti_df[["references_x", "references_y"]].apply(
        self.aggregate_column_level, axis=1)
        
        # merge sources 
        drugbank_plus_chembl_dti_df["source"] = drugbank_plus_chembl_dti_df[["source_x", "source_y"]].apply(
        self.merge_source_column, axis=1)
        
        # merge mechanism of action types
        drugbank_plus_chembl_dti_df["mechanism_of_action_type"] = drugbank_plus_chembl_dti_df[["mechanism_of_action_type_x", "mechanism_of_action_type_y"]].apply(lambda x: str(list(x.dropna())[0]).lower() if list(x.dropna()) else np.nan, axis=1)
        
        # drop redundant columns
        drugbank_plus_chembl_dti_df.drop(columns=["references_x", "references_y", "source_x", "source_y", 
                                          "mechanism_of_action_type_x", "mechanism_of_action_type_y", 
                                          ], inplace=True)
        
        t1 = time()        
        print(f'Drugbank and Chembl DTI data is merged in {round((t1-t0) / 60, 2)} mins')
        
        print("Started merging Drugbank+Chembl and Pharos DTI data")
        t0 = time()
        
        # merge drugbank+chembl and pharos dti
        drugbank_plus_chembl_plus_pharos_dti_df = drugbank_plus_chembl_dti_df.merge(self.pharos_dti_duplicate_removed_df, how="outer", on=["uniprot_id", "drugbank_id"])
        
        # merge references
        drugbank_plus_chembl_plus_pharos_dti_df["references"] = drugbank_plus_chembl_plus_pharos_dti_df[["references_x", "references_y"]].apply(
        self.aggregate_column_level, axis=1)
        
        # merge mechanism of action types
        drugbank_plus_chembl_plus_pharos_dti_df["mechanism_of_action_type"] = drugbank_plus_chembl_plus_pharos_dti_df[["mechanism_of_action_type_x", "mechanism_of_action_type_y"]].apply(
        lambda x: str(x.dropna().tolist()[0]).lower() if len(x.dropna().tolist())>0 else np.nan, axis=1)
        
        # merge activity value
        drugbank_plus_chembl_plus_pharos_dti_df["activity_value"] = drugbank_plus_chembl_plus_pharos_dti_df[["activity_value_x", "activity_value_y"]].apply(
        lambda x: x.dropna().tolist()[0] if len(x.dropna().tolist())>0 else np.nan, axis=1)
        
        # merge activity type
        drugbank_plus_chembl_plus_pharos_dti_df["activity_type"] = drugbank_plus_chembl_plus_pharos_dti_df[["activity_type_x", "activity_type_y"]].apply(
        lambda x: x.dropna().tolist()[0] if len(x.dropna().tolist())>0 else np.nan, axis=1)
        
        # merge sources
        drugbank_plus_chembl_plus_pharos_dti_df["source"] = drugbank_plus_chembl_plus_pharos_dti_df[["source_x", "source_y"]].apply(
        self.merge_source_column, axis=1)
        
        # drop redundant columns
        drugbank_plus_chembl_plus_pharos_dti_df.drop(columns=["source_x", "source_y", "references_x", "references_y",
                                                     "activity_value_x", "activity_value_y", "activity_type_x", "activity_type_y",
                                                      "mechanism_of_action_type_x", "mechanism_of_action_type_y"], inplace=True)
        
        t1 = time()        
        print(f'Drugbank+Chembl and Pharos DTI data is merged in {round((t1-t0) / 60, 2)} mins')
        
        print("Started merging Drugbank+Chembl+Pharos and Dgidb DTI data")
        t0 = time()
        
        # merge drugbank+chembl+pharos and dgidb
        drugbank_plus_chembl_plus_pharos_plus_dgidb_dti_df = drugbank_plus_chembl_plus_pharos_dti_df.merge(self.dgidb_dti_duplicate_removed_df, how="outer", on=["uniprot_id", "drugbank_id"])
        
        # merge references
        drugbank_plus_chembl_plus_pharos_plus_dgidb_dti_df["references"] = drugbank_plus_chembl_plus_pharos_plus_dgidb_dti_df[["references_x", "references_y"]].apply(self.aggregate_column_level, axis=1)
        
        # merge mechanism of action types
        drugbank_plus_chembl_plus_pharos_plus_dgidb_dti_df["mechanism_of_action_type"] = drugbank_plus_chembl_plus_pharos_plus_dgidb_dti_df[["mechanism_of_action_type_x", "mechanism_of_action_type_y"]].apply(
        lambda x: x.dropna().tolist()[0] if len(x.dropna().tolist())>0 else np.nan, axis=1)
        
        # merge sources
        drugbank_plus_chembl_plus_pharos_plus_dgidb_dti_df["source"] = drugbank_plus_chembl_plus_pharos_plus_dgidb_dti_df[["source_x", "source_y"]].apply(
        self.merge_source_column, axis=1)
        
        # drop redundant columns
        drugbank_plus_chembl_plus_pharos_plus_dgidb_dti_df.drop(columns=["source_x", "source_y", "references_x", "references_y", 
                                                                 "mechanism_of_action_type_x", "mechanism_of_action_type_y"], inplace=True)
        
        t1 = time()        
        print(f'Drugbank+Chembl+Pharos and Dgidb DTI data is merged in {round((t1-t0) / 60, 2)} mins')
        
        print("Started merging Drugbank+Chembl+Pharos+Dgidb and Stitch DTI data")
        t0 = time()
        
        # merge drugbank+chembl+pharos+dgidb and stitch
        drugbank_plus_chembl_plus_pharos_plus_dgidb_plus_stitch_dti_df = drugbank_plus_chembl_plus_pharos_plus_dgidb_dti_df.merge(self.stitch_dti_duplicate_removed_df, how="outer", on=["uniprot_id", "drugbank_id"])
        
        # merge sources
        drugbank_plus_chembl_plus_pharos_plus_dgidb_plus_stitch_dti_df["source"] = drugbank_plus_chembl_plus_pharos_plus_dgidb_plus_stitch_dti_df[["source_x", "source_y"]].apply(
        self.merge_source_column, axis=1)
        
        # drop redundant columns
        drugbank_plus_chembl_plus_pharos_plus_dgidb_plus_stitch_dti_df.drop(columns=["source_x", "source_y"], inplace=True)
        
        t1 = time()        
        print(f'Drugbank+Chembl+Pharos+Dgidb and Stitch DTI data is merged in {round((t1-t0) / 60, 2)} mins')
        
        print("Started merging Drugbank+Chembl+Pharos+Dgidb+Stitch and Kegg DTI data")
        t0 = time()
        
        drugbank_plus_chembl_plus_pharos_plus_dgidb_plus_stitch_plus_kegg_dti_df = drugbank_plus_chembl_plus_pharos_plus_dgidb_plus_stitch_dti_df.merge(self.kegg_dti_df, how="outer", on=["uniprot_id", "drugbank_id"])
        
        # merge sources
        drugbank_plus_chembl_plus_pharos_plus_dgidb_plus_stitch_plus_kegg_dti_df["source"] = drugbank_plus_chembl_plus_pharos_plus_dgidb_plus_stitch_plus_kegg_dti_df[["source_x", "source_y"]].apply(
        self.merge_source_column, axis=1)
        
        # drop redundant columns
        drugbank_plus_chembl_plus_pharos_plus_dgidb_plus_stitch_plus_kegg_dti_df.drop(columns=["source_x", "source_y"], inplace=True)
        
        t1 = time()        
        print(f'Drugbank+Chembl+Pharos+Dgidb+Stitch and Kegg DTI data is merged in {round((t1-t0) / 60, 2)} mins')
        
        # create self object from final dataframe
        self.all_dti_df = drugbank_plus_chembl_plus_pharos_plus_dgidb_plus_stitch_plus_kegg_dti_df
        
        
    def merge_all_ddis(self):
        
        print("Started merging Kegg and DDInter DDI data")
        t0 = time()

        # merge kegg and ddinter ddi data
        kegg_plus_ddinter_ddi_df = self.kegg_ddi_duplicate_removed_df.merge(self.ddinter_duplicate_removed_df, how="outer", on=["drug1", "drug2"])
        
        # merge source columns
        kegg_plus_ddinter_ddi_df["source"] = kegg_plus_ddinter_ddi_df[["source_x", "source_y"]].apply(self.merge_source_column, axis=1)
        
        # merge interaction type columns
        kegg_plus_ddinter_ddi_df["interaction_type"] = kegg_plus_ddinter_ddi_df[["interaction_type_x", "interaction_type_y"]].apply(self.aggregate_column_level, axis=1)
        
        # drop redundant columns
        kegg_plus_ddinter_ddi_df.drop(columns=["source_x", "source_y", "interaction_type_x", "interaction_type_y"], inplace=True)
        
        t1 = time()        
        print(f'Kegg and DDInter DDI data is merged in {round((t1-t0) / 60, 2)} mins')
        
        self.all_ddi_df = kegg_plus_ddinter_ddi_df
        
    def get_drug_nodes(self):
        """
        Merges drug node information from different sources. 
        """

        self.node_list = []

        print('Started writing drug nodes')
        for k, v in tqdm(self.drugbank_drugs.items()):
            drug_id = normalize_curie('drugbank:' + k)
            
            props = {}
            for prop_key, prop_value in v.items():
                if prop_value:
                    props[prop_key.replace(" ", "_").lower()] = prop_value
                    

            self.node_list.append((drug_id, 'drug', props))


    def get_dti_edges(self, early_stopping=500):

        self.dti_edge_list = []
        
        print('Started writing DTI edges')
        for index, row in tqdm(self.all_dti_df.iterrows(), total=self.all_dti_df.shape[0]):
            
            _dict = row.to_dict()
            source = normalize_curie('drugbank:' + _dict["drugbank_id"])
            target = normalize_curie('uniprot:' +_dict["uniprot_id"])

            del _dict["drugbank_id"], _dict["uniprot_id"]

            props = dict()
            for k, v in _dict.items():
                if str(v) != "nan":
                    if isinstance(v, str) and "|" in v:
                        props[str(k).replace(" ","_").lower()] = v.split("|")
                    else:
                        props[str(k).replace(" ","_").lower()] = v


            self.dti_edge_list.append((None, source, target, "targets", props))
            
            if early_stopping and (index+1) == early_stopping:
                break

    def get_dgi_edges(self, early_stopping=500):
        
        self.dgi_edge_list = []
        
        print('Started writing DGI edges')
        for index, row in tqdm(self.ctd_cgi_duplicate_removed_df.iterrows(), total=self.ctd_cgi_duplicate_removed_df.shape[0]):
            _dict = row.to_dict()
            source = normalize_curie('drugbank:' + _dict["drugbank_id"])
            target = normalize_curie('ncbigene:' + _dict["entrez_id"])
            label = "_".join(["drug", _dict["action_type"], "gene"])

            del _dict["drugbank_id"], _dict["entrez_id"], _dict["action_type"]
            
            props = dict()
            for k, v in _dict.items():
                if str(v) != "nan":
                    if isinstance(v, str) and "|" in v:
                        props[str(k).replace(" ","_").lower()] = v.split("|")
                    else:
                        props[str(k).replace(" ","_").lower()] = v


            self.dgi_edge_list.append((None, source, target, label, props))
            
            if early_stopping and (index+1) == early_stopping:
                break

    def get_ddi_edges(self, early_stopping=500):
        self.ddi_edge_list = []
        
        print('Started writing DGI edges')
        
        for index, row in tqdm(self.all_ddi_df.iterrows(), total=self.all_ddi_df.shape[0]):
            _dict = row.to_dict()
            
            source = normalize_curie('drugbank:' + _dict["drug1"])
            target = normalize_curie('drugbank:' + _dict["drug2"])
            
            del _dict["drug1"], _dict["drug2"]
            
            props = dict()
            for k, v in _dict.items():
                if str(v) != "nan":
                    if isinstance(v, str) and "|" in v:
                        props[str(k).replace(" ","_").lower()] = v.split("|")
                    else:
                        props[str(k).replace(" ","_").lower()] = v


            self.ddi_edge_list.append((None, source, target, "drug_interacts_with_drug", props))
            
            if early_stopping and (index+1) == early_stopping:
                break
                
