from pypath.share import curl, settings
from pypath.inputs import interpro
from pypath.utils import go
from contextlib import ExitStack
from bioregistry import normalize_curie


class InterPro:
    """
    Class that downloads InterPro data using pypath and reformats it to be ready
    for import into a BioCypher database.
    """

    def __init__(self):
        pass



    def download_domain_data(
        self, cache=False, debug=False, retries=3,
    ):
        """
        Wrapper function to download InterPro domain data using pypath; used to access
        settings.
        Args:
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

            self.get_interpro_nodes()
            self.get_interpro_edges(page_size = 150) # WARNING: decrease page_size parameter if there is a curl error about timeout in the pypath_log

    def get_interpro_nodes(self):
        """
        Downloads InterPro data from ebi.ac.uk/interpro/ through pypath.
        """

        # create list of nodes
        self.node_list = []

        interpro_entries = interpro.interpro_entries() # returns a list of namedtuples
        interpro_structural_xrefs = interpro.interpro_xrefs(db_type = 'structural')
        interpro_external_xrefs = interpro.interpro_xrefs(db_type = 'external')

        # reformat list of namedtuples so that it will be ready to be imported into BioCypher 
        for entry in interpro_entries:

            interpro_id = normalize_curie("interpro:" + entry.interpro_id)
            node_props = entry._asdict()
            node_props['structural_xrefs'] = interpro_structural_xrefs[entry.interpro_id]
            node_props['external_xrefs'] = interpro_external_xrefs[entry.interpro_id]

            del node_props['interpro_id']

            self.node_list.append((interpro_id, "interpro", node_props))


    def get_interpro_edges(self, page_size = 150):
        """
        Downloads InterPro data from ebi.ac.uk/interpro/ through pypath.
        """

        # create list of edges
        self.edge_list = []

        # DOMAIN-PROTEIN EDGES
        # WARNING: decrease page_size parameter if there is a curl error about timeout in the pypath_log
        interpro_annotations = interpro.interpro_annotations(page_size = page_size, reviewed = True, tax_id = '')
        
        self.protein_domain_edges = []
        for k, v in interpro_annotations.items():
            for annotation in v:
                interpro_id = normalize_curie("interpro:" + annotation.interpro_id)
                uniprot_id = normalize_curie("uniprot:" + k)
                edge_props = {}
                edge_props['locations'] = annotation.locations # TODO annotation.locations is a set. should we convert it to list?
                
                
                self.protein_domain_edges.append((uniprot_id, interpro_id, 'has', {'locations': annotation.locations})) # TODO delete this row after checking data and keep only self.edge_list.append() row
                self.edge_list.append((None, uniprot_id, interpro_id, 'has', {'locations': annotation.locations}))

        # DOMAIN-FUNCTION EDGES
        interpro2go = interpro.interpro_xrefs(db_type = 'go') # returns dict of interpro ids as keys and go term annotations as values
        goa = go.GOAnnotation()
        _go = goa.ontology
        go_aspect_dict = _go.aspect

        domain_function_label_dict = {
            'P': 'involved_in',
            'F': 'enables',
            'C': 'located_in',
        }

        self.domain_function_edges = []
        for k, v in interpro2go.items():
            if v:
                for go_term in v:
                    aspect = go_aspect_dict.get(go_term)
                    edge_label = domain_function_label_dict.get(aspect)
                    interpro_id = normalize_curie("interpro:" + k)
                    go_id = normalize_curie("go:" + go_term)
                    self.domain_function_edges.append((interpro_id, go_id, edge_label)) # TODO delete this row after checking data and keep only self.edge_list.append() row
                    self.edge_list.append((None, interpro_id, go_id, edge_label, {}))

    


        
    


