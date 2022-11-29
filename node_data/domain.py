from pypath.share import curl, settings
from pypath.inputs import interpro
from pypath.utils import go
from contextlib import ExitStack
from bioregistry import normalize_curie
from tqdm import tqdm

class InterPro:
    """
    Class that downloads InterPro data using pypath and reformats it to be ready
    for import into a BioCypher database.
    """

    def __init__(self, cache=False, debug=False, page_size=150, retries=6):
        """        
        Args:
            cache: if True, it uses the cached version of the data, otherwise
            forces download.
            debug: if True, turns on debug mode in pypath.
            retries: number of retries in case of download error.
        """
        self.cache = cache
        self.debug = debug
        self.retries = retries
        self.page_size = page_size


    def download_domain_node_data(self):
        """
        Downloads domain node data from Interpro
        """

        # stack pypath context managers
        with ExitStack() as stack:

            stack.enter_context(settings.context(retries=self.retries))

            if self.debug:
                stack.enter_context(curl.debug_on())

            if not self.cache:
                stack.enter_context(curl.cache_off())
                
            self.interpro_entries = interpro.interpro_entries() # returns a list of namedtuples
            self.interpro_structural_xrefs = interpro.interpro_xrefs(db_type = 'structural')
            self.interpro_external_xrefs = interpro.interpro_xrefs(db_type = 'external')
            
    def download_domain_edge_data(self, cache=False, debug=False, retries=3):
        """
        Downloads Uniprot annotation data from Interpro
        """
        
        with ExitStack() as stack:

            stack.enter_context(settings.context(retries=self.retries))

            if self.debug:
                stack.enter_context(curl.debug_on())

            if not self.cache:
                stack.enter_context(curl.cache_off())
            
            # WARNING: decrease page_size parameter if there is a curl error about timeout in the pypath_log
            self.interpro_annotations = interpro.interpro_annotations(page_size = self.page_size, reviewed = True, tax_id = '')       

    def get_interpro_nodes(self):
        """
        Downloads InterPro data from ebi.ac.uk/interpro/ through pypath.
        """

        # create list of nodes
        self.node_list = []

        # selected attributes
        primary_attributes = ["protein_count", "name", "type", "parent_list", "child_list"]
        external_attributes = ["pfam", "ec", "pdb"]

        for entry in tqdm(self.interpro_entries):
            props = dict()
            interpro_props = entry._asdict()
            domain_id = normalize_curie("interpro:" + entry.interpro_id)

            for element in primary_attributes:
                if interpro_props.get(element):
                    if element == "protein_count":
                        props[element] = int(interpro_props.get(element))
                    else:
                        if isinstance(interpro_props.get(element), list) and len(interpro_props.get(element)) == 1:                        
                            props[element] = interpro_props.get(element)[0]
                        else:
                            props[element] = interpro_props.get(element)

            for element in external_attributes:
                if element == "pfam":
                    if interpro_props.get('member_list').get(str(element).upper()):
                        if len(interpro_props.get('member_list').get(str(element).upper())) == 1:
                            props[element] = interpro_props.get('member_list').get(str(element).upper())[0]
                        else:                        
                            props[element] = interpro_props.get('member_list').get(str(element).upper())
                elif element == "ec":
                    if self.interpro_external_xrefs.get(entry.interpro_id).get(str(element).upper()):
                        if len(self.interpro_external_xrefs.get(entry.interpro_id).get(str(element).upper())) == 1:
                            props[element] = self.interpro_external_xrefs.get(entry.interpro_id).get(str(element).upper())[0]
                        else:
                            props[element] = self.interpro_external_xrefs.get(entry.interpro_id).get(str(element).upper())
                else:
                    if self.interpro_structural_xrefs.get(entry.interpro_id).get(str(element).upper()):
                        if len(self.interpro_structural_xrefs.get(entry.interpro_id).get(str(element).upper())) == 1:
                            props[element] = self.interpro_structural_xrefs.get(entry.interpro_id).get(str(element).upper())[0]
                        else:
                            props[element] = self.interpro_structural_xrefs.get(entry.interpro_id).get(str(element).upper())
                
            self.node_list.append((domain_id, "domain", props))

    def get_interpro_edges(self):
        # create list of edges
        self.edge_list = []

        # DOMAIN-PROTEIN EDGES
        for k, v in tqdm(interpro_annotations.items()):            
            # k -> uniprot id
            for annotation in v:
                interpro_id = normalize_curie("interpro:" + annotation.interpro_id)
                uniprot_id = normalize_curie("uniprot:" + k)                
                self.edge_list.append((None, uniprot_id, interpro_id, "Has", {"locations":list(annotation.locations)}))
