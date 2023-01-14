from __future__ import annotations

from pypath.share import curl, settings
from pypath.inputs import chembl, uniprot
from contextlib import ExitStack
from bioregistry import normalize_curie

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
    
    def get_compound_nodes(self):
        """
        Reformats compound node data to be ready for import into a BioCypher database.
        """
        print('Writing compound nodes...')
        self.nodes = []
        for compound in tqdm(self.compounds):
            if compound.chembl:
                compound_id = normalize_curie('chembl:' + compound.chembl)
                props = {
                    'chembl_id': compound.chembl,
                    'type': compound.type,
                    'full_mwt': compound.full_mwt,
                    'species': compound.species,
                    'heavy_atoms': compound.heavy_atoms,
                    'alogp': compound.alogp,
                }
                self.nodes.append((compound_id, 'Compound', props))
    
    def get_compound_target_edges(self):
        """
        Reformats compound-target edge data to be ready for import into a BioCypher database.
        """

        print('Writing compound-target edges...')
        self.edges = []
        self.compound_target_edges = []
        for act in tqdm(self.chembl_acts):
            # ids
            compound_id = normalize_curie('chembl:' + act.chembl)
            target_uniprot_id = self.target_dict.get(act.target_chembl)
            if target_uniprot_id:
                target_id = normalize_curie('uniprot:' + target_uniprot_id)

                pubmed_id = self.document_to_pubmed.get(act.document)
                
                assay = self.assay_dict.get(act.assay_chembl)
                if assay:
                    # add only binding assays
                    assay_type = assay.assay_type
                    assay_conf = assay.confidence_score

                    props = {
                        'pchembl': act.pchembl,
                        'standard_value': act.standard_value,
                        'standard_type': act.standard_type,
                        'assay_chembl' : act.assay_chembl,
                        'assay_type' : assay_type,
                        'assay_conf' : assay_conf,
                        'pubmed_id' : pubmed_id,
                    }

                    self.compound_target_edges.append((None, compound_id, target_id, 'Targets', props))
                    self.edges.append((None, compound_id, target_id, 'Targets', props))
