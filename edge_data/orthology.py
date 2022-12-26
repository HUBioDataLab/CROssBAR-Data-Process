from __future__ import annotations

from pypath.share import curl, settings
from pypath.inputs import oma, uniprot, pharos
from pypath.utils import taxonomy
from contextlib import ExitStack
from bioregistry import normalize_curie
from time import time

from biocypher._logger import logger
from tqdm import tqdm


OMA_ORGANISMS = {
    9606,
    10090,
    3702,
    10116,
    559292,
    9913,
    1264690,
    83333,
    6239,
    1423,
    39947,
    44689,
    7227,
    8355,
    7955,
    9031,
    1773,
    9601,
}.union(set(taxonomy.taxids.keys()))


class Orthology:
    """
    Class that downloads orthology data using pypath and reformats it to be ready
    for import into a BioCypher database.
    """

    def __init__(self):
        pass

    
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
    

    def download_oma_data(self, tax=OMA_ORGANISMS):
        """
        Downloads orthology data from OMA.

        Args
            tax: list of taxids to download data for.
        """

        self.entry_name_to_uniprot = uniprot.uniprot_data(field= 'entry name', reviewed = True, organism= '*')
        self.entry_name_to_uniprot = {v:k for k,v in self.entry_name_to_uniprot.items()}

        logger.debug("Started downloading OMA orthology data")
        t0 = time()
        
        self.oma_orthology = []

        for t in tqdm(tax):
            tax_orthology = oma.oma_orthologs(organism_a = 9606, organism_b = t)
            tax_orthology = [i for i in tax_orthology if i.id_a in self.entry_name_to_uniprot and i.id_b in self.entry_name_to_uniprot]
            self.oma_orthology.extend(tax_orthology)

        t1 = time()
        logger.info(f'OMA orthology data is downloaded in {round((t1-t0) / 60, 2)} mins')

    
    def download_pharos_data(self):
        """
        Downloads orthology data from Pharos.
        """

        logger.debug("Started downloading Pharos orthology data")
        t0 = time()

        uniprot_to_entrez = uniprot.uniprot_data(field= 'database(GeneID)', reviewed = True, organism= '*')
        self.entrez_to_uniprot = {}
        for k,v in uniprot_to_entrez.items():
            for entrez in v.split(';'):
                if entrez != '':
                    self.entrez_to_uniprot[entrez] = k
   
        pharos_orthology_init = pharos.pharos_targets(orthologs=True)
        self.pharos_orthology = []
        for protein in pharos_orthology_init:
            if protein['orthologs']:
                for ortholog in protein['orthologs']:
                    if str(ortholog['geneid']) in self.entrez_to_uniprot:
                        self.pharos_orthology.append(
                            {
                                'id_a': protein['uniprot'],
                                'id_b': self.entrez_to_uniprot[str(ortholog['geneid'])],
                                'db_id_b' : ortholog['dbid'],
                            })

        t1 = time()
        logger.info(f'Pharos orthology data is downloaded in {round((t1-t0) / 60, 2)} mins')

    
    def get_orthology_edges(self):
        """
        Reformats orthology data to be ready for import into a BioCypher database.
        """

        self.edge_list = []

        self.oma_orthology_edges = []
        logger.info("Preparing OMA orthology edges.")
        for orthology in self.oma_orthology:
            id_a = normalize_curie("uniprot:" + self.entry_name_to_uniprot[orthology.id_a])
            id_b = normalize_curie("uniprot:" + self.entry_name_to_uniprot[orthology.id_b])
            props = {
                'rel_type': orthology.rel_type,
                'score': orthology.score,
            }

            self.oma_orthology_edges.append((None, id_a, id_b, 'Is_Ortholog_To', props))
            self.edge_list.append((None, id_a, id_b, 'Is_Ortholog_To', props))

        self.pharos_orthology_edges = []
        logger.info("Preparing Pharos orthology edges.")
        for orthology in self.pharos_orthology:
            id_a = normalize_curie("uniprot:" + orthology['id_a'])
            id_b = normalize_curie("uniprot:" + orthology['id_b'])

            self.pharos_orthology_edges.append((None, id_a, id_b, 'Is_Ortholog_To', {'db_id_b': orthology['db_id_b']}))
            self.edge_list.append((None, id_a, id_b, 'Is_Ortholog_To', {'db_id_b': orthology['db_id_b']}))
