from pypath.share import curl, settings
from pypath.inputs import interpro, uniprot
from pypath.inputs import go as go_input
from pypath.utils import go as go_util
from contextlib import ExitStack
from bioregistry import normalize_curie


class GO:
    """
    Class that downloads Gene Ontology data using pypath and reformats it to be ready
    for import into a BioCypher database.
    """

    def __init__(self, ):
        pass


    def download_go_data(
        self, organism = '*', cache=False, debug=False, retries=3,
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

            self.go_ontology = go_util.GeneOntology()
            self.go_annots = go_input.go_annotations_all(organism = organism, fields= ['qualifier', 'go_id', 'reference', 'evidence_code']) # returns dict of uniprot ids as keys and go term annotations as values
            self.interpro2go = interpro.interpro_xrefs(db_type = 'go') # returns dict of interpro ids as keys and go term annotations as values
            self.swissprots = list(uniprot._all_uniprots(organism = organism, swissprot =True))

    def get_go_nodes(self):
        """
        Downloads Gene Ontology term data from ebi.ac.uk/QuickGO through pypath.
        """

        # create list of nodes
        self.node_list = []

        for go_term in self.go_ontology.name.keys(): #keys in self.go_ontology.name is current go ids in the ontology
            
            go_id = normalize_curie("go:" + go_term)

            node_props = {}
            node_props['name'] = self.go_ontology.name[go_term]
            node_props['aspect'] = self.go_ontology.aspect[go_term]

            self.node_list.append((go_id, "go", node_props))

    
    def get_go_edges(self):
        """
        Downloads various edges of GO term nodes through pypath.
        """

        self.edge_list = []

        # PROTEIN-FUNCTION EDGES
        # Protein annotation data from ebi.ac.uk/GOA
        self.protein_function_edges = []
        for k, v in self.go_annots.items():
            if k in self.swissprots:
                protein_id = normalize_curie("uniprot:" + k)
                for annotation in list(v):
                    if annotation.evidence_code != 'IEA': # filter electronic annotations
                        go_id = normalize_curie("go:" + annotation.go_id)
                        edge_label = annotation.qualifier
                        props = {
                            'reference': annotation.reference,
                            'evidence_code': annotation.evidence_code,
                        }
                        self.protein_function_edges.append((protein_id, go_id, edge_label, props)) # TODO delete this row after checking data and keep only self.edge_list.append() line
                        self.edge_list.append((None, protein_id, go_id, edge_label, props))


        # FUNCTION-FUNCTION EDGES
        self.function_function_edges = []
        for k, v in self.go_ontology.ancestors.items():
            source_go_id = normalize_curie("go:" + k)
            
            for ancestor in list(v):
                target_go_id = normalize_curie("go:" + ancestor[0])
                edge_label = ancestor[1]
                self.function_function_edges.append((source_go_id, target_go_id, edge_label)) # TODO delete this row after checking data and keep only self.edge_list.append() line
                self.edge_list.append((None, source_go_id, target_go_id, edge_label, {}))


        # DOMAIN-FUNCTION EDGES
        go_aspect_dict = self.go_ontology.aspect

        domain_function_label_dict = {
            'P': 'involved_in',
            'F': 'enables',
            'C': 'located_in',
        }

        self.domain_function_edges = []
        for k, v in self.interpro2go.items():
            if v:
                for go_term in v:
                    aspect = go_aspect_dict.get(go_term)
                    edge_label = domain_function_label_dict.get(aspect)
                    interpro_id = normalize_curie("interpro:" + k)
                    go_id = normalize_curie("go:" + go_term)
                    self.domain_function_edges.append((interpro_id, go_id, edge_label)) # TODO delete this row after checking data and keep only self.edge_list.append() line
                    self.edge_list.append((None, interpro_id, go_id, edge_label, {}))
