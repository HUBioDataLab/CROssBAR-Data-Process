from pypath.share import curl, settings
from pypath.inputs import drugbank, drugcentral
from contextlib import ExitStack
from bioregistry import normalize_curie
from tqdm import tqdm

class Drug:
    """
    Class that downloads drug data using pypath and reformats it to be ready
    for import into a BioCypher database.
    """

    def __init__(self, ):

        pass


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

            # Drugbank Object
            drugbank_data = drugbank.DrugbankFull(user = user, passwd = passwd)
            
            # NODES
            self.drugbank_drugs = {}

            drugbank_drugs_external_ids = drugbank.drugbank_drugs(user = user, passwd = passwd)
            drugbank_drugs_external_ids = {drug.drugbank: drug._asdict() for drug in drugbank_drugs_external_ids}
            drugbank_external_fields = ['drugbank', 'kegg_compound', 'kegg_drug', 'pubchem_cid', 'pubchem_sid', 'chebi', 'chembl', 'pharmgkb', 'het']

            for k, v in drugbank_drugs_external_ids.items():
                drugbank_drugs_external_ids[k] = {key: value for key, value in v.items() if key in drugbank_external_fields}
            
            drugbank_drugs_detailed = drugbank_data.drugbank_drugs_full(fields = ['type', 'cas_number', 'name', 'groups', 'general_references', 'atc_codes', 'drug_interactions'])

            for drug in drugbank_drugs_detailed:
                drugbank_id = drug._asdict()['drugbank_id']
                self.drugbank_drugs[drugbank_id] = (
                    drug._asdict() | drugbank_drugs_external_ids[drugbank_id]
                        if drugbank_id in drugbank_drugs_external_ids else
                    drug._asdict()
                )

                del self.drugbank_drugs[drugbank_id]['drugbank_id']

            # self.drugcentral_drugs = drugcentral.drugcentral_drugs()

            # EDGES
            
            # DRUG-TARGET
            self.drugbank_dti = drugbank_data.drugbank_targets_full(fields=['drugbank_id', 'id', 'actions', 'references', 'known_action', 'polypeptide',])


    def get_drug_nodes(self):
        """
        Merges drug node information from different sources. 
        """

        self.node_list = []
        self.drugbank_node_list = [] #remove this after final adjustments

        for k,v in tqdm(self.drugbank_drugs.items()):
            drug_id = normalize_curie('drugbank:' + k)
            props = v
            del props['drug_interactions']
            self.drugbank_node_list.append((drug_id, props))
            self.node_list.append((drug_id, 'drug', props))


    def get_drug_edges(self):
        self.edge_list = []

        # DRUG-TARGET

        # drugbank
        self.drugbank_dti_list = [] #remove this after final adjustments
        for dti in tqdm(self.drugbank_dti):
            if dti.polypeptide: # filter out targets without uniprot id
                if type(dti.polypeptide) == list:
                    for polypeptide in dti.polypeptide:
                        if polypeptide[1] == 'Swiss-Prot': # filter out Trembl IDs 
                            drug_id = normalize_curie('drugbank:' + dti.drugbank_id)
                            target_id = normalize_curie('uniprot:' + polypeptide[0])
                            props = dti._asdict()
                            del props['drugbank_id']
                            del props['polypeptide']
                            self.drugbank_dti_list.append((drug_id, target_id, 'Targets', props))
                            self.edge_list.append((None, drug_id, target_id, 'Targets', props))
                
                elif type(dti.polypeptide) == tuple:
                    if dti.polypeptide[1] == 'Swiss-Prot':
                        drug_id = normalize_curie('drugbank:' + dti.drugbank_id)
                        target_id = normalize_curie('uniprot:' + dti.polypeptide[0])
                        props = dti._asdict()
                        del props['drugbank_id']
                        del props['polypeptide']
                        self.drugbank_dti_list.append((drug_id, target_id, 'Targets', props))
                        self.edge_list.append((None, drug_id, target_id, 'Targets', props))

                else:
                    print(type(dti.polypeptide))
                    break
                
        print('Drugbank DTI are written')

        # DRUG-DRUG

        # drugbank
        self.drugbank_ddi_list = [] #remove this after final adjustments
        for k,v in tqdm(self.drugbank_drugs.items()):
            if v['drug_interactions']:
                for target in v['drug_interactions']:
                    source_drug_id = normalize_curie('drugbank:' + k)
                    target_drug_id = normalize_curie('drugbank:' + target)
                    self.drugbank_ddi_list.append((source_drug_id, target_drug_id, 'Interacts_with'))
                    self.edge_list.append((None, source_drug_id, target_drug_id, 'Interacts_with'))
        print('Drugbank DDI are written')
        