
from __future__ import annotations

from pypath.share import curl, settings, common
from pypath.inputs import drugbank, drugcentral, stitch, string, uniprot
from contextlib import ExitStack
from typing import Literal
from bioregistry import normalize_curie
from tqdm import tqdm
from time import time
import collections

class Drug:
    """
    WARNING: Drugbank, Drug Central and STITCH database download urls contain
        version number/date. Please update this urls (or request an update) in 
        resources.urls module of the pypath library before using this script
        in order to access the latest versions of data.

    Class that downloads drug data using pypath and reformats it to be ready
    for import into a BioCypher database.
    """

    def __init__(self, ):
        
        self.edge_list = []


    def download_drug_data(
        self, user, passwd, cache=False, debug=False, retries=3, 
        ):

        """
        Wrapper function to download drug data from various databases using pypath.
        Args:
            user: drugbank username
            passwd: drugbank password
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


            print('Started downloading and processing Drugbank data for drug nodes')
            # Drugbank Object
            drugbank_data = drugbank.DrugbankFull(user = user, passwd = passwd)
            
            # NODES
            
            self.drugbank_drugs = {}
            drugbank_drugs_external_ids = drugbank.drugbank_drugs(user = user, passwd = passwd)
            drugbank_drugs_external_ids = {drug.drugbank: drug._asdict() for drug in drugbank_drugs_external_ids}
            drugbank_external_fields = ['drugbank', 'kegg_compound', 'kegg_drug', 'pubchem_cid', 'pubchem_sid', 'chebi', 'chembl', 'pharmgkb', 'het']            
            drugbank_drugs_detailed = drugbank_data.drugbank_drugs_full(fields = ['type', 'cas_number', 'name', 'groups', 'general_references', 'atc_codes', 'drug_interactions'])

            all_fields = list(drugbank_drugs_detailed[0]._fields) + drugbank_external_fields
            for drug in drugbank_drugs_detailed:
                drugbank_id = drug.drugbank_id
                temp_dict = (
                    drug._asdict() | drugbank_drugs_external_ids[drugbank_id]
                        if drugbank_id in drugbank_drugs_external_ids else
                    drug._asdict()
                )

                self.drugbank_drugs[drugbank_id] = {f: temp_dict.get(f, None) for f in all_fields}

                del self.drugbank_drugs[drugbank_id]['drugbank_id']
                del self.drugbank_drugs[drugbank_id]['drugbank'] 

            print('Started downloading and processing Drug Central data for drug nodes')
            self.drugcentral_drugs = {}
            drugcentral_drugs_init = drugcentral.drugcentral_drugs()
            self.cas_to_drugbank = {value['cas_number']: key for key, value in self.drugbank_drugs.items()}
            for drug in drugcentral_drugs_init:
                if drug.cas in self.cas_to_drugbank:
                    drugbank_id = self.cas_to_drugbank[drug.cas]
                    self.drugcentral_drugs[drugbank_id] = drug._asdict()

            # EDGES
            
            # DRUG-TARGET
            print('Downloading DTI data: Drugbank')
            self.drugbank_dti = drugbank_data.drugbank_targets_full(fields=['drugbank_id', 'id', 'actions', 'references', 'known_action', 'polypeptide',])
            print('Downloading DTI data: Drug Central')
            self.drugcentral_dti = drugcentral.drugcentral_interactions()
            # print('Downloading DTI data: KEGG')
            # self.kegg_dti = 

            # for stitch we need to iterate through organisms in string db, this will take ~40 hours
            self.download_stitch_dti_data()
            # pubchem to drugbank id mapping, stitch returns pubchem CIDs as drug IDs
            self.db_pubchem_mapping = drugbank.drugbank_mapping(id_type='pubchem_compound', target_id_type='drugbank', user=user, passwd=passwd)


    def download_stitch_dti_data(
            self, 
            organism: str | list = None, 
            score_threshold: int | Literal[
                'highest_confidence',
                'high_confidence',
                'medium_confidence',
                'low_confidence',
                ] = 'high_confidence', 
            physical_interaction_score: bool = False # currently this arg doesnt work for organisms other than human. can be fixed if necessary.
            ):
        """
        Wrapper function to download STITCH DTI data using pypath
            
        It returns PubChem CIDs as drug IDs and STRING IDs as protein IDs.
        """

        if organism is None:
            organism = string.string_species()

        organism = common.to_list(organism)

        # map string ids to swissprot ids
        uniprot_to_string = uniprot.uniprot_data("database(STRING)", "*", True)
        self.string_to_uniprot = collections.defaultdict(list)
        for k,v in uniprot_to_string.items():
            for string_id in list(filter(None, v.split(";"))):
                self.string_to_uniprot[string_id.split(".")[1]].append(k)


        self.stitch_ints = []
        print("Started downloading STITCH data")
        t0 = time()

        for tax in tqdm(organism):
            try:
                organism_stitch_ints = [
                    i for i in stitch.stitch_links_interactions(ncbi_tax_id=int(tax), score_threshold=score_threshold, physical_interaction_score= physical_interaction_score)
                    if i.partner_b in self.string_to_uniprot] # filter with swissprot ids

                # TODO: later engage below print line to biocypher log 
                # print(f"Downloaded STITCH data with taxonomy id {str(tax)}")

                if organism_stitch_ints:
                    self.stitch_ints.extend(organism_stitch_ints)
            
            except TypeError: #'NoneType' object is not an iterator
                pass
                # TODO: later engage below print line to biocypher log 
                # print(f'Skipped tax id {tax}. This is most likely due to the empty file in database. Check the database file anyway.')

        t1 = time()
        print(f'STITCH data is downloaded in {round((t1-t0) / 60, 2)} mins')


    def get_drug_nodes(self):
        """
        Merges drug node information from different sources. 
        """

        self.node_list = []

        print('Started writing drug nodes')
        for k,v in tqdm(self.drugbank_drugs.items()):
            drug_id = normalize_curie('drugbank:' + k)
            props = v

            # add mapped drug central id if there is any
            
            props['drugcentral'] = (
                self.drugcentral_drugs[k]['drugcentral']
                    if k in self.drugcentral_drugs else
                None
                )

            del props['drug_interactions']

            self.node_list.append((drug_id, 'drug', props))

    def get_dti_edges(self):
        
        # drugbank
        print('Started writing Drugbank DTI')
        self.drugbank_dti_list = [] #remove this after final adjustments
        for dti in tqdm(self.drugbank_dti):
            if dti.polypeptide: # filter out targets without uniprot id
                if type(dti.polypeptide) == list:
                    for polypeptide in dti.polypeptide:
                        if polypeptide[1] == 'Swiss-Prot': # filter out Trembl IDs 
                            drug_id = normalize_curie('drugbank:' + dti.drugbank_id)
                            target_id = normalize_curie('uniprot:' + polypeptide[0])
                            props = dti._asdict()
                            props['source'] = 'DrugBank'
                            del props['drugbank_id']
                            del props['polypeptide']
                            self.drugbank_dti_list.append((drug_id, target_id, 'Targets', props))
                            self.edge_list.append((None, drug_id, target_id, 'Targets', props))
                
                elif type(dti.polypeptide) == tuple:
                    if dti.polypeptide[1] == 'Swiss-Prot':
                        drug_id = normalize_curie('drugbank:' + dti.drugbank_id)
                        target_id = normalize_curie('uniprot:' + dti.polypeptide[0])
                        props = dti._asdict()
                        props['source'] = 'DrugBank'
                        del props['drugbank_id']
                        del props['polypeptide']
                        self.drugbank_dti_list.append((drug_id, target_id, 'Targets', props))
                        self.edge_list.append((None, drug_id, target_id, 'Targets', props))

        # drug central
        print('Started writing Drug Central DTI')
        self.drugcentral_dti_list = []
        for dti in tqdm(self.drugcentral_dti):
            if dti.drug:
                if dti.drug.cas in self.cas_to_drugbank:
                    drugbank_id = self.cas_to_drugbank[dti.drug.cas]
                    drug_id = normalize_curie('drugbank:' + drugbank_id)
                    target_id = normalize_curie('uniprot:' + dti.uniprot)

                    drugcentral_dti_fields = ['target_type', 'act_value', 'act_type', 'relation', 'effect', 'tdl']
                    props = {key: value for key, value in dti._asdict().items() if key in drugcentral_dti_fields}
                    props['source'] = 'Drug Central'
                    self.drugcentral_dti_list.append((drug_id, target_id, 'Targets', props))
                    self.edge_list.append((None, drug_id, target_id, 'Targets', props))


        # stitch
        print('Started writing STITCH DTI')
        self.stitch_dti_list = [] #remove this after final adjustments

        for dti in tqdm(self.stitch_ints):
                
            drug_id = self.db_pubchem_mapping.get(dti.partner_a)
            if drug_id:
                drug_id = normalize_curie('drugbank:' + list(drug_id)[0])
                # convert string to uniprot id 
                target_id = normalize_curie('uniprot:' + self.string_to_uniprot[dti.partner_b][0])

                self.stitch_dti_list.append((drug_id, target_id, 'Targets', {'combined_score': dti.combined_score}))
                self.edge_list.append((None, drug_id, target_id, 'Targets', {'combined_score': dti.combined_score}))


    def get_ddi_edges(self):

        # drugbank
        self.drugbank_ddi_list = [] #remove this after final adjustments
        for k,v in tqdm(self.drugbank_drugs.items()):
            if v['drug_interactions']:
                for target in v['drug_interactions']:
                    source_drug_id = normalize_curie('drugbank:' + k)
                    target_drug_id = normalize_curie('drugbank:' + target)
                    self.drugbank_ddi_list.append((source_drug_id, target_drug_id, 'Interacts_with', {'source': 'DrugBank'}))
                    self.edge_list.append((None, source_drug_id, target_drug_id, 'Interacts_with', {'source': 'DrugBank'}))
                    
        print('Drugbank DDI are written')
        