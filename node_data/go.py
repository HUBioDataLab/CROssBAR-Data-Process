from pypath.share import curl, settings
from pypath.inputs import interpro, uniprot
from pypath.inputs import go as go_input
from pypath.utils import go as go_util
from contextlib import ExitStack
from bioregistry import normalize_curie

from time import time

from tqdm.notebook import tqdm

from biocypher._logger import logger

logger.debug(f"Loading module {__name__}.")

class GO:
    """
    Class that downloads Gene Ontology data using pypath and reformats it to be ready
    for import into a BioCypher database.
    """

    def __init__(self, ):
        pass


    def download_go_data(
        self, organism = 9606, cache=False, debug=False, retries=6,
    ):
        """
        Wrapper function to download Gene Ontology data using pypath; used to access
        settings.
        Args:
            organism: ncbi tax id or known name of organism of interest
            cache: if True, it uses the cached version of the data, otherwise
            forces download.
            debug: if True, turns on debug mode in pypath.
            retries: number of retries in case of download error.
        """

        # stack pypath context managers
        with ExitStack() as stack:

            stack.enter_context(settings.context(retries=retries))

            if debug:
                stack.enter_context(curl.debug_on())

            if not cache:
                stack.enter_context(curl.cache_off())
            
            logger.debug("Started downloading Gene Ontology entry data")
            t0 = time()
            
            self.go_ontology = go_util.GeneOntology()
            
            t1 = time()
            logger.info(f'Gene Ontology entry data is downloaded in {round((t1-t0) / 60, 2)} mins')
            
            logger.debug("Started downloading Gene Ontology annotation data")
            t0 = time()
            
            self.go_annots = go_input.go_annotations_all(organism = organism, fields= ['qualifier', 'go_id', 'reference', 'evidence_code']) # returns dict of uniprot ids as keys and go term annotations as values
            
            t1 = time()
            logger.info(f'Gene Ontology annotation data is downloaded in {round((t1-t0) / 60, 2)} mins')
            
            logger.debug("Started downloading Interpro2go data")
            t0 = time()            
            
            self.interpro2go = interpro.interpro_xrefs(db_type = 'go') # returns dict of interpro ids as keys and go term annotations as values
            t1 = time()
            logger.info(f'Interpro2go data is downloaded in {round((t1-t0) / 60, 2)} mins')
            
            self.swissprots = list(uniprot._all_uniprots(organism = organism, swissprot =True))

    def get_go_nodes(self):
        """
        Downloads Gene Ontology term data from ebi.ac.uk/QuickGO through pypath.
        """
        logger.info("Preparing nodes.")
        
        # create list of nodes
        self.node_list = []
        
        aspect_to_node_label_dict = {"C":"cellular component", "P":"biological process", "F":"molecular function"}
        for go_term in tqdm(self.go_ontology.name.keys()): #keys in self.go_ontology.name is current go ids in the ontology
            
            go_id = normalize_curie("go:" + go_term)

            node_props = {}
            node_props['name'] = self.go_ontology.name[go_term]
            label = aspect_to_node_label_dict[self.go_ontology.aspect[go_term]]

            self.node_list.append((go_id, label, node_props))

    
    def get_go_edges(self):
        """
        Downloads various edges of GO term nodes through pypath.
        """
        logger.info("Preparing Protein-GO edges.")
        
        self.edge_list = []

        # PROTEIN-GO EDGES
        # Protein annotation data from ebi.ac.uk/GOA
        self.protein_to_go_edges = []
        for k, v in tqdm(self.go_annots.items()):
            if k in self.swissprots:
                protein_id = normalize_curie("uniprot:" + k)
                for annotation in list(v):
                    # filter electronic annotations and qualifiers with "NOT" and the ones that not in go ontology
                    if annotation.go_id in self.go_ontology.aspect.keys() and annotation.evidence_code != 'IEA' and not str(annotation.qualifier).startswith("NOT"): 
                        go_id = normalize_curie("go:" + annotation.go_id)
                        edge_label = str(annotation.qualifier).title()
                        props = {
                            'reference': annotation.reference,
                            'evidence_code': annotation.evidence_code,
                        }
                        self.protein_to_go_edges.append((None, protein_id, go_id, edge_label, props)) # TODO delete this row after checking data and keep only self.edge_list.append() line
                        self.edge_list.append((None, protein_id, go_id, edge_label, props))

        logger.info("Preparing GO-GO edges.")
        
        # GO-GO EDGES
        self.go_to_go_edges = []
        for k, v in tqdm(self.go_ontology.ancestors.items()):
            source_go_id = normalize_curie("go:" + k)
            
            for ancestor in list(v):
                target_go_id = normalize_curie("go:" + ancestor[0])
                edge_label = str(ancestor[1]).title()
                self.go_to_go_edges.append((None, source_go_id, target_go_id, edge_label)) # TODO delete this row after checking data and keep only self.edge_list.append() line
                self.edge_list.append((None, source_go_id, target_go_id, edge_label, {}))

        
        logger.info("Preparing Domain-GO edges.")
        
        # DOMAIN-GO EDGES    
        domain_function_label_dict = {
            'P': 'Involved_in',
            'F': 'Enables',
            'C': 'Located_in',
        }

        self.domain_to_go_edges = []
        for k, v in tqdm(self.interpro2go.items()):
            if v:
                for go_term in v:
                    aspect = self.go_ontology.aspect.get(go_term)
                    edge_label = domain_function_label_dict.get(aspect)
                    interpro_id = normalize_curie("interpro:" + k)
                    go_id = normalize_curie("go:" + go_term)
                    self.domain_to_go_edges.append((None, interpro_id, go_id, edge_label)) # TODO delete this row after checking data and keep only self.edge_list.append() line
                    self.edge_list.append((None, interpro_id, go_id, edge_label, {}))
