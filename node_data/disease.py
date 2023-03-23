
from __future__ import annotations

from pypath.share import curl, settings, common
from pypath.inputs import drugbank, drugcentral, stitch, string, uniprot, dgidb, pharos, ctdbase
from contextlib import ExitStack
from typing import Literal
from bioregistry import normalize_curie
from tqdm import tqdm
from time import time
import collections

class Disease:
    """
    Class that downloads disease data using pypath and reformats it to be ready
    for import into a BioCypher database.
    """

    def __init__(self):
        """
        Args
            some_arg: Some argument
        """

        self.edge_list = []


    def download_disease_data(
        self, cache=False, debug=False, retries=3, 
        ):

        """
        Wrapper function to download disease data from various databases using pypath.

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

            # TODO self.download_xxx_data()

    def download_disease_template(self):

        """
        Wrapper function to download Template disease entries, DTI and DDI data using pypath
        Args
            some_arg: Some argument
        """

        # node data
        print('Downloading disease node data: Template')
        # TODO self.template_diseases = template_diseases()

        # edge data
        print('Downloading DTI data: Template')
        # TODO self.template_dti = template_dti()


    def get_disease_nodes(self):
        """
        Merges disease node information from different sources. 
        """

        self.node_list = []

        print('Started writing drug nodes')
        # TODO Do some stuff

    def get_dti_edges(self):
        
        # template_1
        print('Started writing Template DTI')

        # template_2
        print('Started writing Template_2 DTI')