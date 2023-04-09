from __future__ import annotations

from pypath.share import curl, settings, common
from pypath.inputs import (
    pathophenodb,
    ctdbase,
    clinvar,
    disgenet,
    kegg_api,
    pharos,
    chembl,
)
from contextlib import ExitStack
from typing import Literal
from bioregistry import normalize_curie
from tqdm import tqdm
from time import time
import collections

kegg_org = "hsa"


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

        self.node_list = []
        self.edge_list = []
        self.download_all_data()

    def download_all_data(self):
        self.download_disease_data()
        # self.download_kegg_disease_nodes()
        self.download_efo_mondo_disease_nodes()
        self.download_disease_is_a_disease_efo_edges()
        self.download_disease_is_a_disease_mondo_edges()
        self.download_disgenet_disease_disease_edges()
        self.download_clinvar_gene_disease_edges()
        self.download_disgenet_gene_disease_edges()
        # self.download_kegg_gene_disease_edges()
        self.download_jensenlab_diseases_gene_disease_edges()
        # self.download_ctdbase_disease_gene_edges()
        # self.download_ctdbase_disease_chemical_edges()
        self.download_opentargets_disease_protein_edges()
        self.download_pharos_disease_protein_edges()
        self.download_pathopheno_disease_organism_edges()
        # self.download_ctdbase_disease_pathway_edges()
        # self.download_kegg_disease_drug_edges()
        self.download_chembl_disease_drug_edges()

    def download_disease_data(
        self,
        cache=False,
        debug=False,
        retries=3,
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

    def download_disease_template(self):
        """
        Wrapper function to download Template disease entries, DTI and DDI data using pypath
        Args
            some_arg: Some argument
        """

        # node data
        print("Downloading disease node data: Template")
        # TODO self.template_diseases = template_diseases()

        # edge data
        print("Downloading DTI data: Template")
        # TODO self.template_dti = template_dti()

    def get_disease_nodes(self):
        """
        Merges disease node information from different sources.
        """

        self.node_list = []

        print("Started writing disease nodes")

    def download_kegg_disease_nodes(self):
        disease_ids = kegg_api._Disease()
        pass

    def download_efo_mondo_disease_nodes(self):
        # TODO
        pass

    def download_disease_is_a_disease_efo_edges(self):
        # TODO
        pass

    def download_disease_is_a_disease_mondo_edges(self):
        # TODO
        pass

    def download_disgenet_disease_disease_edges(self):

        if not hasattr(self, 'disgenet_api'):
            self.disgenet_api = disgenet.DisgenetApi()

        if (not hasattr(self, 'disgenet_disease_ids')) or not self.disgenet_disease_ids:
            self.disgenet_disease_ids = disgenet.disease_id_mappings().keys()
            print(self.disgenet_disease_ids)
            print(type(self.disgenet_disease_ids))

        self.disgenet_dda_gene_edges = []
        self.disgenet_dda_variant_edges = []

        for disease_id in self.disgenet_disease_ids:
            try:
                self.disgenet_dda_gene_edges.extend(
                    self.disgenet_api.get_ddas_that_share_genes(disease_id)
                )
                self.disgenet_dda_variant_edges.extend(
                    self.disgenet_api.get_ddas_that_share_variants(disease_id)
                )
            except TypeError:
                print(f'{disease_id} not available')

    def download_clinvar_gene_disease_edges(self):
        print("Started writing ClinVAR GD")
        self.clinvar_gd_edges = clinvar.clinvar_raw()

    def download_disgenet_gene_disease_edges(self):

        if not hasattr(self, 'disgenet_api'):
            self.disgenet_api = disgenet.DisgenetApi()
        
        if (not hasattr(self, 'disgenet_disease_ids')) or not self.disgenet_disease_ids:
            self.disgenet_disease_ids = disgenet.disease_id_mappings().keys()
            print(self.disgenet_disease_ids)
            print(type(self.disgenet_disease_ids))

        self.disgenet_gda_edges = []

        for disease_id in self.disgenet_disease_ids:

            try:
                self.disgenet_gda_edges.extend(
                    self.disgenet_api.get_gdas_by_diseases(disease_id)
                )

            except TypeError:
                print(f'{disease_id} not available')

    def download_kegg_gene_disease_edges(self):
        self.kegg_disease_to_gene = kegg_api.disease_to_gene(kegg_org)
        self.kegg_gene_to_disease = kegg_api.gene_to_disease(kegg_org)

    def download_jensenlab_diseases_gene_disease_edges(self):
        
        # TODO Not implemented to pypath
        self.jensenlabs_knowledge_filtered = diseases.knowledge_filtered()
        self.jensenlabs_experimental_filtered = diseases.experimental_filtered()

    def download_ctdbase_disease_gene_edges(self):
        print("Started writing CTDBase DG")
        self.ctdbase_dp_edges = ctdbase.ctdbase_relations("gene_disease")

    def download_ctdbase_disease_chemical_edges(self):
        print("Started writing CTDBase DC")
        self.dp_edges = ctdbase.ctdbase_relations("chemical_disease")
        print(self.dp_edges)

    def download_opentargets_disease_protein_edges(self):
        # TODO Open Targets is not merged yet
        self.opentargets_indirect = opentargets.overall_indirect_score()
        self.opentargets_direct = opentargets.overall_direct_score()

    def download_pharos_disease_protein_edges(self):
        self.pharos_edges = pharos.pharos_targets(diseases=True)

    def download_pathopheno_disease_organism_edges(self):
        print("Started writing PathoPheno DO")
        self.pathopheno_do_edges = pathophenodb.disease_pathogen_interactions()

    def download_ctdbase_disease_pathway_edges(self):
        print("Started writing CTDBase DP")
        self.ctdbase_dp_edges = ctdbase.ctdbase_relations("disease_pathway")

    def download_kegg_disease_drug_edges(self):
        # FIXME KEGG
        self.kegg_disease_to_drug = kegg_api.disease_to_drug()
        self.kegg_drug_to_disease = kegg_api.drug_to_disease()
        pass

    def download_chembl_disease_drug_edges(self):
        self.chembl_disease_drug = chembl.chembl_drug_indications()


Disease()
