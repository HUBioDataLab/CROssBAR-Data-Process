class PathwayEdgeType(Enum):
    PROTEIN_TO_PATHWAY = auto()
    REACTOME_HIERARCHICAL_RELATIONS = auto()
    DRUG_TO_PATHWAY = auto()
    DISEASE_TO_PATHWAY = auto()
    PATHWAY_TO_PATHWAY = auto()
    

class Pathway:
    def __init__(self, drugbank_user, drugbank_passwd, edge_types: Union[list[PathwayEdgeType], None] = None,
                add_prefix = True, kegg_organism="hsa"):
        
        self.drugbank_user = drugbank_user
        self.drugbank_passwd = drugbank_passwd
        self.add_prefix = add_prefix
        
        if not kegg_organism:
            self.kegg_organism = list(kegg_local._Organism()._data.keys())
        else:
            self.kegg_organism = self.ensure_iterable(kegg_organism)
        
        self.set_edge_types(edge_types=edge_types)
    
    def download_pathway_data(
        self,
        cache=False,
        debug=False,
        retries=3,
    ):
        """
        Wrapper function to download pathway data from various databases using pypath.
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
                
            
            self.download_reactome_data()
            
            self.download_kegg_data()
            
            self.download_ctd_data()
            
            self.download_compath_data()
    
    def download_reactome_data(self):
        print("Started downloading Reactome data")
        t0 = time()
        
        self.reactome_pathways = reactome.reactome_data("pathways")
        
        if PathwayEdgeType.REACTOME_HIERARCHICAL_RELATIONS in self.edge_types:            
            self.reactome_hierarchial_relations = reactome.reactome_data("pathway_relations")
        
        if PathwayEdgeType.PROTEIN_TO_PATHWAY in self.edge_types:            
            self.reactome_uniprot_pathway = reactome.reactome_data("uniprot2all")
            
        if PathwayEdgeType.DRUG_TO_PATHWAY in self.edge_types:
            self.reactome_chebi_pathway = reactome.reactome_data("chebi2all")
            self.chebi_to_drugbank = {list(v)[0]:k for k, v in unichem.unichem_mapping("drugbank", "chebi").items()}
        
        t1 = time()
        print(f"Reactome data is downloaded in {round((t1-t0) / 60, 2)} mins")
        
    def download_kegg_data(self):
        print("Started downloading KEGG data")
        t0 = time()
        
        self.kegg_pathway_abbv_organism_name_dict = {k:v[1] for k, v in kegg_local._Organism()._data.items()}
        
        self.kegg_pathways = []
        for org in tqdm(self.kegg_organism):
            try:
                self.kegg_pathways.extend(kegg_local._kegg_list("pathway", org=org))
            except:
                print(f"Error occured in {org} organism data")
        
        if PathwayEdgeType.PROTEIN_TO_PATHWAY in self.edge_types:
            self.kegg_gene_to_pathway = {}
            for org in tqdm(self.kegg_organism):                
                self.kegg_gene_to_pathway = self.kegg_gene_to_pathway | kegg_local.gene_to_pathway(org=org)
                
            self.kegg_to_uniprot = {v.strip(";").split(";")[0]:k for k, v in uniprot.uniprot_data("database(KEGG)", 9606, True).items()}
            
        if PathwayEdgeType.DRUG_TO_PATHWAY in self.edge_types:
            self.kegg_drug_to_pathway = kegg_local.drug_to_pathway()
            drugbank_data = drugbank.DrugbankFull(user = self.drugbank_user, passwd = self.drugbank_passwd)
            drugbank_drugs_external_ids = drugbank_data.drugbank_external_ids_full()
            self.kegg_drug_to_drugbank = {v.get("KEGG Drug"):k for k, v in drugbank_drugs_external_ids.items() if v.get("KEGG Drug")}
            
        if PathwayEdgeType.DISEASE_TO_PATHWAY in self.edge_types:
            
            self.kegg_disease_to_pathway = kegg_local.disease_to_pathway()
            
            self.kegg_diseases_mappings = {}
            for dis in self.kegg_disease_to_pathway.keys():
                result = kegg_local.get_diseases(dis)
                self.kegg_diseases_mappings[dis] = result[0].db_links
                
        t1 = time()
        print(f"KEGG data is downloaded in {round((t1-t0) / 60, 2)} mins")
                
    def download_ctd_data(self):
        print("Started downloading CTD data")
        t0 = time()
        
        if PathwayEdgeType.DISEASE_TO_PATHWAY in self.edge_types:
            self.ctd_disease_pathway = ctdbase.ctdbase_relations("disease_pathway")
            
        t1 = time()
        print(f"CTD data is downloaded in {round((t1-t0) / 60, 2)} mins")
        
    def download_compath_data(self):
        print("Started downloading Compath data")
        t0 = time()
        
        if PathwayEdgeType.PATHWAY_TO_PATHWAY in self.edge_types:            
            self.compath_pathway_pathway = compath.get_all_mappings(source_db=None, target_db=None)
        
        t1 = time()
        print(f"CTD data is downloaded in {round((t1-t0) / 60, 2)} mins")
        
    def process_reactome_protein_pathway(self):
        if not hasattr(self, "reactome_uniprot_pathway"):
            self.download_reactome_data()
            
        print("Started processing Reactome protein-pathway data")
        t0 = time()
        
        df_list = []
        for pp in self.reactome_uniprot_pathway:
            if pp.evidence_code not in ["IEA"]:                
                df_list.append((pp.uniprot_id, pp.pathway_id))
            
        df = pd.DataFrame(df_list, columns=["uniprot_id", "pathway_id"])
        
        df.drop_duplicates(ignore_index=True, inplace=True)
        
        df["source"] = "Reactome"
        
        t1 = time()
        print(f"Reactome protein-pathway data is processed in {round((t1-t0) / 60, 2)} mins")
        
        return df
    
    def process_kegg_protein_pathway(self):
        if not hasattr(self, "kegg_gene_to_pathway"):
            self.download_kegg_data()
            
        print("Started processing KEGG protein-pathway data")
        t0 = time()
        
        df_list = []
        for gene, pathways in self.kegg_gene_to_pathway.items():
            if not gene.startswith("org") and self.kegg_to_uniprot.get(gene):
                organism_prefix = gene.split(":")[0].strip()
                for pathway in pathways.PathwayEntries:
                    df_list.append((self.kegg_to_uniprot[gene], pathway.pathway_id.replace("map", organism_prefix)))
                    
        df = pd.DataFrame(df_list, columns=["uniprot_id", "pathway_id"])
        
        df.drop_duplicates(ignore_index=True, inplace=True)
        
        df["source"] = "KEGG"
        
        t1 = time()
        print(f"KEGG protein-pathway data is processed in {round((t1-t0) / 60, 2)} mins")
        
        return df
    
    def process_reactome_drug_pathway(self):
        if not hasattr(self, "reactome_chebi_pathway"):
            self.download_reactome_data()
            
        print("Started processing Reactome drug-pathway data")
        t0 = time()
        
        df_list = []
        for cp in self.reactome_chebi_pathway:
            if cp.evidence_code not in ["IEA"] and self.chebi_to_drugbank.get(cp.chebi_id):
                df_list.append((self.chebi_to_drugbank[cp.chebi_id], cp.pathway_id))
                
        df = pd.DataFrame(df_list, columns=["drug_id", "pathway_id"])
        
        df.drop_duplicates(ignore_index=True, inplace=True)
        
        df["source"] = "Reactome"
        
        t1 = time()
        print(f"Reactome drug-pathway data is processed in {round((t1-t0) / 60, 2)} mins")
        
        return df
    
    def process_kegg_drug_pathway(self):        
        if not hasattr(self, "kegg_drug_to_pathway"):
            self.download_kegg_data()
            
        print("Started processing KEGG drug-pathway data")
        t0 = time()
        
        df_list = []
        for drug, pathways in self.kegg_drug_to_pathway.items():
            if self.kegg_drug_to_drugbank.get(drug):
                for pathway in pathways.PathwayEntries:
                    df_list.append((self.kegg_drug_to_drugbank[drug], pathway.pathway_id.replace("map", "hsa"),))
                    
        df = pd.DataFrame(df_list, columns=["drug_id", "pathway_id"])
        
        df.drop_duplicates(ignore_index=True, inplace=True)
        
        df["source"] = "KEGG"
        
        t1 = time()
        print(f"KEGG drug-pathway data is processed in {round((t1-t0) / 60, 2)} mins")
        
        return df
    
    def process_kegg_disease_pathway(self):
        if not hasattr(self, "kegg_disease_to_pathway"):
            self.download_kegg_data()        
        if not hasattr(self, "mondo_mappings"):
            self.prepare_mondo_mappings()
            
        print("Started processing KEGG disease-pathway data")
        t0 = time()
        
        kegg_dbs_to_mondo_dbs = {"MeSH":"MESH", "OMIM":"OMIM", "ICD-10":"ICD10CM",}
        
        df_list = []
        for disease, pathways in self.kegg_disease_to_pathway.items():
            found = False
            disease_id = None
            for db in kegg_dbs_to_mondo_dbs.keys():
                if found:
                    break
                if self.kegg_diseases_mappings[disease].get(db):
                    for ref in self.ensure_iterable(self.kegg_diseases_mappings[disease][db]):
                        if found:
                            break

                        if self.mondo_mappings[kegg_dbs_to_mondo_dbs[db]].get(ref):
                            disease_id = self.mondo_mappings[kegg_dbs_to_mondo_dbs[db]].get(ref)

            if disease_id:            
                for pathway in pathways.PathwayEntries:
                    df_list.append((disease_id, pathway.pathway_id.replace("map", "hsa"),))
                    
        
        df = pd.DataFrame(df_list, columns=["disease_id", "pathway_id"])
        
        df.drop_duplicates(ignore_index=True, inplace=True)
        
        df["source"] = "KEGG"
        
        t1 = time()
        print(f"KEGG disease-pathway data is processed in {round((t1-t0) / 60, 2)} mins")
        
        return df
    
    def process_ctd_disease_pathway(self):
        if not hasattr(self, "ctd_disease_pathway"):
            self.download_ctd_data()        
        if not hasattr(self, "mondo_mappings"):
            self.prepare_mondo_mappings()
       
        kegg_pathways_checker_list = [i[0] for i in kegg_local._kegg_list("pathway", org="hsa")]
            
        print("Started processing CTD disease-pathway data")
        t0 = time()
        
        df_list = []
        for dp in self.ctd_disease_pathway:
            disease_db = dp.DiseaseID.split(":")[0]
            pathway_db = dp.PathwayID.split(":")[0]
            if self.mondo_mappings[disease_db].get(dp.DiseaseID.split(":")[1]):
                if pathway_db == "KEGG":
                    if dp.PathwayID.split(":")[1].replace("_","") in kegg_pathways_checker_list:
                        df_list.append((self.mondo_mappings[disease_db].get(dp.DiseaseID.split(":")[1]),
                                       dp.PathwayID.split(":")[1].replace("_",""),))
                else:
                    df_list.append((self.mondo_mappings[disease_db].get(dp.DiseaseID.split(":")[1]),
                                       dp.PathwayID.split(":")[1],))
                    
        df = pd.DataFrame(df_list, columns=["disease_id", "pathway_id"])
        
        df.drop_duplicates(ignore_index=True, inplace=True)
        
        df["source"] = "CTD"
        
        t1 = time()
        print(f"CTD disease-pathway data is processed in {round((t1-t0) / 60, 2)} mins")
        
        return df
    
    def merge_protein_pathway_data(self):
        
        kegg_df = self.process_kegg_protein_pathway()
        
        reactome_df = self.process_reactome_protein_pathway()
        
        print("Started merging protein-pathway data")
        t0 = time()
        
        merged_df = pd.concat([kegg_df, reactome_df], ignore_index=True)
            
        t1 = time()
        print(f"protein-pathway edge data is merged in {round((t1-t0) / 60, 2)} mins")
        
        return merged_df
    
    def merge_drug_pathway_data(self):
        
        kegg_df = self.process_kegg_drug_pathway()
        
        reactome_df = self.process_reactome_drug_pathway()
        
        print("Started merging drug-pathway data")
        t0 = time()
        
        merged_df = pd.concat([kegg_df, reactome_df], ignore_index=True)
        
        t1 = time()
        print(f"Drug-pathway edge data is merged in {round((t1-t0) / 60, 2)} mins")
        
        return merged_df
    
    def merge_disease_pathway_data(self):
        
        kegg_df = self.process_kegg_disease_pathway()
        
        ctd_df = self.process_ctd_disease_pathway()
        
        print("Started merging disease-pathway data")
        t0 = time()
        
        merged_df = pd.merge(kegg_df, ctd_df, how="outer", on=["disease_id", "pathway_id"])
        
        merged_df["source"] = merged_df[["source_x", "source_y"]].apply(self.merge_source_column, axis=1)
        
        merged_df.drop(columns=["source_x", "source_y"], inplace=True)
        
        t1 = time()
        print(f"Drug-pathway edge data is merged in {round((t1-t0) / 60, 2)} mins")
        
        return merged_df
        
    def get_nodes(self, label="pathway"):
        if not hasattr(self, "reactome_pathways"):
            self.download_reactome_data()
        if not hasattr(self, "kegg_pathways"):
            self.download_kegg_data()
            
        print("Started writing pathway nodes")
            
        node_list = []
        
        for p in tqdm(self.reactome_pathways):
            pathway_id = self.add_prefix_to_id(prefix="reactome", identifier=p.pathway_id)                
            node_list.append((pathway_id, label, {"name":p.pathway_name, "organism":p.organism}))
            
        for p in tqdm(self.kegg_pathways):
            pathway_id = self.add_prefix_to_id(prefix="kegg.pathway", identifier=p[0])
            node_list.append((pathway_id, label, {"name":p[1].split("-")[0].strip(), "organism":self.kegg_pathway_abbv_organism_name_dict.get(p[0][:3]) }))            
            
        return node_list
    
    def get_edges(self):
        
        print("Started writing all pathway edges")
        edge_list = []
        
        edge_list.extend(self.get_protein_pathway_edges())
        
        edge_list.extend(self.get_drug_pathway_edges())
        
        edge_list.extend(self.get_disease_pathway_edges())
        
        edge_list.extend(self.get_pathway_pathway_edges())
        
        edge_list.extend(self.get_reactome_hierarchical_edges())
        
        return edge_list
        
    def get_protein_pathway_edges(self, label="protein_take_part_in_pathway"):
        
        protein_pathway_edges_df = self.merge_protein_pathway_data()
        
        print("Started writing protein-pathway edges")
        
        edge_list = []
        for index, row in tqdm(protein_pathway_edges_df.iterrows(), total=protein_pathway_edges_df.shape[0]):
            _dict = row.to_dict()
            
            if _dict["source"] == "Reactome":                
                pathway_id = self.add_prefix_to_id(prefix="reactome", identifier=_dict["pathway_id"])
            else:
                pathway_id = self.add_prefix_to_id(prefix="kegg.pathway", identifier=_dict["pathway_id"])
            
            uniprot_id = self.add_prefix_to_id(prefix="uniprot", identifier=_dict["uniprot_id"])
            
            del _dict["uniprot_id"], _dict["pathway_id"]
            
            props = {}
            for k, v in _dict.items():
                if str(v) != "nan":
                    if isinstance(v, str) and "|" in v:
                        props[k] = v.split("|")
                    else:
                        props[k] = v
            
            edge_list.append((None, uniprot_id, pathway_id, label, props))
            
        
        return edge_list
    
    def get_drug_pathway_edges(self, label="drug_has_target_in_pathway"):
        
        drug_pathway_edges_df = self.merge_drug_pathway_data()
        
        print("Started writing drug-pathway edges")
        
        edge_list = []
        for index, row in tqdm(drug_pathway_edges_df.iterrows(), total=drug_pathway_edges_df.shape[0]):
            _dict = row.to_dict()
            
            if _dict["source"] == "Reactome":                
                pathway_id = self.add_prefix_to_id(prefix="reactome", identifier=_dict["pathway_id"])
            else:
                pathway_id = self.add_prefix_to_id(prefix="kegg.pathway", identifier=_dict["pathway_id"])
                
            drug_id = self.add_prefix_to_id(prefix="drugbank", identifier=_dict["drug_id"])
            
            del _dict["drug_id"], _dict["pathway_id"]
            
            props = {}
            for k, v in _dict.items():
                if str(v) != "nan":
                    if isinstance(v, str) and "|" in v:
                        props[k] = v.split("|")
                    else:
                        props[k] = v
            
            edge_list.append((None, drug_id, pathway_id, label, props))
            
        return edge_list
    
    def get_disease_pathway_edges(self, label="disease_modulates_pathway"):
        
        disease_pathway_edges_df = self.merge_disease_pathway_data()
        
        print("Started writing disease-pathway edges")
        
        edge_list = []
        for index, row in tqdm(disease_pathway_edges_df.iterrows(), total=disease_pathway_edges_df.shape[0]):
            _dict = row.to_dict()
            
            if _dict["pathway_id"].startswith("R-"):
                pathway_id = self.add_prefix_to_id(prefix="reactome", identifier=_dict["pathway_id"])
            else:
                pathway_id = self.add_prefix_to_id(prefix="kegg.pathway", identifier=_dict["pathway_id"])
                
            disease_id = self.add_prefix_to_id(prefix="MONDO", identifier=_dict["disease_id"])
            
            del _dict["disease_id"], _dict["pathway_id"]
            
            props = {}
            for k, v in _dict.items():
                if str(v) != "nan":
                    if isinstance(v, str) and "|" in v:
                        props[k] = v.split("|")
                    else:
                        props[k] = v
            
            edge_list.append((None, disease_id, pathway_id, label, props))
            
        return edge_list
    
    def get_pathway_pathway_edges(self):
        if not hasattr(self, "compath_pathway_pathway"):
            self.download_compath_data()
            
        print("Started writing pathway-pathway edges")
        
        edge_list = []
        for pp in tqdm(self.compath_pathway_pathway):
            if pp.source_db in ["kegg", "reactome"] and pp.target_db in ["kegg", "reactome"]:
                if pp.relation == "isPartOf":
                    label = "pathway_is_part_of_pathway"
                elif pp.relation == "equivalentTo":
                    label = "pathway_is_equivalent_to_pathway"
                    
                if pp.id_1.startswith("R-"):
                    pathway_id1 = self.add_prefix_to_id(prefix="reactome", identifier=pp.id_1)
                else:
                    pathway_id1 = self.add_prefix_to_id(prefix="kegg.pathway", identifier=pp.id_1)
                    
                if pp.id_2.startswith("R-"):
                    pathway_id2 = self.add_prefix_to_id(prefix="reactome", identifier=pp.id_2)
                else:
                    pathway_id2 = self.add_prefix_to_id(prefix="kegg.pathway", identifier=pp.id_2)
                    
                edge_list.append((None, pathway_id1, pathway_id2, label, {}))
                
        return edge_list
    
    def get_reactome_hierarchical_edges(self, label="pathway_participates_pathway"):
        if not hasattr(self, "reactome_hierarchial_relations"):
            self.download_reactome_data()
        
        print("Started writing reactome hierarchial edges")
        
        edge_list = []
        for pp in tqdm(self.reactome_hierarchial_relations):
            parent_id = self.add_prefix_to_id(prefix="reactome", identifier=pp.parent)
            child_id = self.add_prefix_to_id(prefix="reactome", identifier=pp.child)
            
            edge_list.append((None, child_id, parent_id, label, {}))
            
        return edge_list
    
    def prepare_mondo_mappings(self):        
        mondo = ontology.ontology(ontology="mondo", fields=["is_obsolete", "obo_xref"])
        
        self.mondo_mappings = collections.defaultdict(dict)
        mapping_db_list = ["MESH", "OMIM", "ICD10CM"]

        for term in mondo:
            if not term.is_obsolete and term.obo_id and "MONDO" in term.obo_id and term.obo_xref:
                for xref in term.obo_xref:
                    if xref.get("database") in mapping_db_list:                
                        db = mapping_db_list[mapping_db_list.index(xref.get("database"))]
                        self.mondo_mappings[db][xref["id"]] = term.obo_id
        
    def add_prefix_to_id(self, prefix=None, identifier=None, sep=":") -> str:
        """
        Adds prefix to database id
        """
        if self.add_prefix and identifier:
            return normalize_curie(prefix + sep + str(identifier))
        
        return identifier
    
    def merge_source_column(self, element, joiner="|"):
        
        _list = []
        for e in list(element.dropna().values):
            if joiner in e:
                for i in e.split(joiner):
                    _list.append(i)
            else:
                _list.append(e)
                
        return joiner.join(list(dict.fromkeys(_list).keys()))
    
    def ensure_iterable(self, element):
        if isinstance(element, (list, tuple, set)):
            return element
        else:
            return [element]
    
    def set_edge_types(self, edge_types):
        if edge_types:
            self.edge_types = edge_types
        else:
            self.edge_types = [field for field in PathwayEdgeType]
