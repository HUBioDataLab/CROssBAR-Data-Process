
from __future__ import annotations

from pypath.share import curl, settings, common
from pypath.inputs import drugbank, drugcentral, stitch, string, uniprot, dgidb, pharos
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

    def __init__(self, drugbank_user, drugbank_passwd):
        """
        Args
            drugbank_user: drugbank username
            drugbank_passwd: drugbank password
        """

        self.edge_list = []
        self.user = drugbank_user
        self.passwd = drugbank_passwd


    def download_drug_data(
        self, cache=False, debug=False, retries=3, 
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

            self.download_drugbank_data(self.user, self.passwd)
            self.drugcentral_drugbank_mapping(self.user, self.passwd)
            self.download_drugcentral_data() # TODO: we can add organism filter for DTI
            self.download_dgidb_data(self.user, self.passwd) # TODO: we can add organism filter for DTI
            self.download_pharos_data()

            # for stitch we need to iterate through organisms in string db, this will take ~40 hours
            self.download_stitch_dti_data(self.user, self.passwd) # TODO: we can add organism filter, score_threshold and physical_interaction_score args


    def download_drugbank_data(self, user, passwd):

        """
        Wrapper function to download DrugBank drug entries, DTI and DDI data using pypath
        Args
            user: drugbank username
            passwd: drugbank password
        """
        
        print('Downloading drug node data: Drugbank')
        # Drugbank Object
        self.drugbank_data = drugbank.DrugbankFull(user = user, passwd = passwd)

        # node data
        drugbank_drugs_external_ids = drugbank.drugbank_drugs(user = user, passwd = passwd)
        drugbank_drugs_detailed = self.drugbank_data.drugbank_drugs_full(fields = ['type', 'cas_number', 'name', 'groups', 'general_references', 'atc_codes', 'drug_interactions'])

        self.drugbank_drugs = {}
        drugbank_external_fields = ['drugbank', 'kegg_compound', 'kegg_drug', 'pubchem_cid', 'pubchem_sid', 'chebi', 'chembl', 'pharmgkb', 'het']            
        all_fields = list(drugbank_drugs_detailed[0]._fields) + drugbank_external_fields
        
        drugbank_drugs_external_ids = {drug.drugbank: drug._asdict() for drug in drugbank_drugs_external_ids}
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

        # edge data
        print('Downloading DTI data: Drugbank')
        self.drugbank_dti = self.drugbank_data.drugbank_targets_full(fields=['drugbank_id', 'id', 'actions', 'references', 'known_action', 'polypeptide',])


    def drugcentral_drugbank_mapping(self, user, passwd):
        # cas_to_drugbank = drugbank.drugbank_mapping(id_type='cas', target_id_type='drugbank', user=user, passwd=passwd)
        cas_to_drugbank = {value['cas_number']: key for key, value in self.drugbank_drugs.items()}
        drugcentral_to_cas = drugcentral.drugcentral_mapping(id_type='drugcentral', target_id_type='cas')
        self.drugcentral_to_drugbank = collections.defaultdict(list)
        for k, v in drugcentral_to_cas.items():
            if list(v)[0] in cas_to_drugbank:
                self.drugcentral_to_drugbank[k] = cas_to_drugbank[list(v)[0]]


    def download_drugcentral_data(
        self,
        organism: str | int | None = None):

        """
        Wrapper function to download Drug Central drug entries and DTI data using pypath
        Args
            organism: 
        """

        print('Downloading drug node data: Drug Central')
        # node data
        drugcentral_drugs_init = drugcentral.drugcentral_drugs()
        self.drugcentral_drugs = {}
        for drug in drugcentral_drugs_init:
            if drug.drugcentral in self.drugcentral_to_drugbank:
                drugbank_id = self.drugcentral_to_drugbank[drug.drugcentral]
                self.drugcentral_drugs[drugbank_id] = drug._asdict()

        print('Downloading DTI data: Drug Central')

        # edge data
        self.drugcentral_dti = drugcentral.drugcentral_interactions(organism)
    

    def download_dgidb_data(self, user, passwd):

        """
        Wrapper function to download DGIdb DTI data using pypath
        
        Args
            user: drugbank username
            passwd: drugbank password

        It returns ChEMBL IDs as drug IDs and Entrez Gene IDs as protein IDs.
        """

        # edge data
        print('Downloading DTI data: DGIdb')
        self.dgidb_dti = dgidb.dgidb_interactions()

        # map chembl ids to drugbank ids
        self.db_chembl_mapping = drugbank.drugbank_mapping(id_type='chembl', target_id_type='drugbank', user=user, passwd=passwd)
        
        # map entrez gene ids to swissprot ids
        uniprot_to_entrez = uniprot.uniprot_data("database(GeneID)", "*", True)
        self.entrez_to_uniprot = collections.defaultdict(list)
        for k,v in uniprot_to_entrez.items():
            for entrez_id in list(filter(None, v.split(";"))):
                self.entrez_to_uniprot[entrez_id].append(k)


    def download_kegg_data(self):
        pass

    
    def download_pharos_data(self):

        print('Downloading DTI data: Pharos')
        self.pharos_dti = pharos.pharos_targets(ligands=True)


    def download_stitch_dti_data(
            self, 
            user,
            passwd,
            organism: str | list = None, 
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

        Args
            user: drugbank username
            passwd: drugbank password
            organism: Name or NCBI Taxonomy IDs of organisms
            score_threshold: Minimum required interaction score. user can use
                pre-defined confidence limits or can define a custom value.
            physical_interaction_score: If True, returns physical interaction scores of interactions.

        Returns PubChem CIDs as drug IDs and STRING IDs as protein IDs.
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

        # mapping to convert pubchem ids to drugbank ids while writing dti edges
        self.db_pubchem_mapping = drugbank.drugbank_mapping(id_type='pubchem_compound', target_id_type='drugbank', user=user, passwd=passwd)
        
        print(f'STITCH data is downloaded in {round((t1-t0) / 60, 2)} mins')


    def get_drug_nodes(self):
        """
        Merges drug node information from different sources. 
        """

        self.node_list = []

        print('Started writing drug nodes')
        for k,v in tqdm(self.drugbank_drugs.items()):
            drug_id = normalize_curie('drugbank:' + k)
            props = {key: None if value == '' else value for key, value in v.items()}

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
        self.drugcentral_dti_list = [] #remove this after final adjustments
        for dti in tqdm(self.drugcentral_dti):
            if dti.drug:
                if dti.drug.drugcentral in self.drugcentral_to_drugbank:
                    drugbank_id = self.drugcentral_to_drugbank[dti.drug.drugcentral]
                    drug_id = normalize_curie('drugbank:' + drugbank_id)
                    target_id = normalize_curie('uniprot:' + dti.uniprot)

                    drugcentral_dti_fields = ['target_type', 'act_value', 'act_type', 'relation', 'effect', 'tdl']
                    props = {key: value for key, value in dti._asdict().items() if key in drugcentral_dti_fields}
                    props['source'] = 'Drug Central'

                    self.drugcentral_dti_list.append((drug_id, target_id, 'Targets', props))
                    self.edge_list.append((None, drug_id, target_id, 'Targets', props))

        # dgidb
        print('Started writing DGIdb DTI')
        self.dgidb_dti_list = [] #remove this after final adjustments
        for dti in tqdm(self.dgidb_dti):
            if dti.drug_chembl and 'chembl' in dti.drug_chembl:
                chembl_id = dti.drug_chembl.split('chembl:')[1]
                entrez_id = dti.entrez
                drugbank_id = self.db_chembl_mapping.get(chembl_id)
                uniprot_id = self.entrez_to_uniprot.get(entrez_id)
                if drugbank_id and uniprot_id:
                    drug_id = normalize_curie('drugbank:' + list(drugbank_id)[0])
                    target_id = normalize_curie('uniprot:' + list(uniprot_id)[0])

                    dgidb_dti_fields = ['type', 'score', 'pmid']
                    props = {key: value for key, value in dti._asdict().items() if key in dgidb_dti_fields}
                    props['source'] = 'DGIdb'

                    self.dgidb_dti_list.append((drug_id, target_id, 'Targets', props))
                    self.edge_list.append((None, drug_id, target_id, 'Targets', props))

        # pharos
        print('Started writing Pharos DTI')
        self.pharos_dti_list = [] #remove this after final adjustments
        for dti in tqdm(self.pharos_dti):
            if dti['ligands']:
                uniprot_id = dti['uniprot']
                for ligand in dti['ligands']:
                    drugcentral_id = [synonym['value'] for synonym in ligand['synonyms'] if synonym['name'] == 'DrugCentral'][0]
                    drugbank_id = self.drugcentral_to_drugbank.get(drugcentral_id)
                    if drugbank_id:
                        drug_id = normalize_curie('drugbank:' + drugbank_id)
                        target_id = normalize_curie('uniprot:' + uniprot_id)

                        activity_list = []
                        no_change_keys = ['actid', 'type', 'moa', 'value']
                        for activity in ligand['activities']:
                            temp_activity_dict = {k: activity.get(k, None) for k in no_change_keys}
                            pubmed_ids = (
                                [i['pmid'] for i in activity['pubs']]
                                    if activity['pubs'] else
                                None
                            )
                            temp_activity_dict['pubs'] = pubmed_ids
                            activity_list.append(temp_activity_dict)

                        self.pharos_dti_list.append((drug_id, target_id, 'Targets', {'activities': activity_list}))
                        self.edge_list.append((None, drug_id, target_id, 'Targets', {'activities': activity_list}))
            
        # stitch
        print('Started writing STITCH DTI')
        self.stitch_dti_list = [] #remove this after final adjustments

        for dti in tqdm(self.stitch_ints):
                
            drugbank_id = self.db_pubchem_mapping.get(dti.partner_a)
            uniprot_id = self.string_to_uniprot.get(dti.partner_b)
            if drugbank_id and uniprot_id:
                drug_id = normalize_curie('drugbank:' + list(drugbank_id)[0])
                target_id = normalize_curie('uniprot:' + list(uniprot_id)[0])
                self.stitch_dti_list.append((drug_id, target_id, 'Targets', {'combined_score': dti.combined_score, 'source': 'STITCH'}))
                self.edge_list.append((None, drug_id, target_id, 'Targets', {'combined_score': dti.combined_score, 'source': 'STITCH'}))


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
        