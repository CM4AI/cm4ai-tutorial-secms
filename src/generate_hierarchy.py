import os
import pandas as pd 

from cellmaps_utils.provenance import ProvenanceUtil
from cellmaps_generate_hierarchy.ppi import CosineSimilarityPPIGenerator
from cellmaps_generate_hierarchy.hierarchy import CDAPSHiDeFHierarchyGenerator
from cellmaps_generate_hierarchy.maturehierarchy import HiDeFHierarchyRefiner
from cellmaps_generate_hierarchy.hcx import HCXFromCDAPSCXHierarchy
from cellmaps_generate_hierarchy.runner import CellmapsGenerateHierarchy

from load_raw_data import load_secms_data, get_disambiguated_dfs, get_formatted_df


cell_line = "MDA-MB-468"
category = "CTRL"
PPI_CUTOFFS = [0.001, 0.002, 0.003]

base_path = "data/raw"
embedding_rocrate_path = 'data/results/embedding'
hierarchy_path = 'data/results/hierarchydir'

os.makedirs(embedding_rocrate_path, exist_ok=True)

secms_data = load_secms_data(base_path)
secms_data_clean, secms_data_ambig = get_disambiguated_dfs(secms_data)
secms_formatted = get_formatted_df(secms_data_clean, cell_line, category)

embedding_file_path = os.path.join(embedding_rocrate_path, 'coembedding_emd.tsv')

prov_utils = ProvenanceUtil(raise_on_error=True)
if not os.path.exists(embedding_file_path):
    secms_formatted.to_csv(embedding_file_path, sep='\t', index=False)
    # Generate provenance for Embedding File
    crate_details = f'SEC-MS Embedding Data {cell_line} {category}'
    prov_utils.register_rocrate(embedding_rocrate_path,
                                name=crate_details,
                                organization_name='CM4AI',
                                project_name=crate_details,
                                description=crate_details,
                                keywords=['secms', cell_line, category, 'embedding'])
    
# Generator that creates edge lists used as input for HiDeF
ppigen = CosineSimilarityPPIGenerator(embeddingdirs=[embedding_rocrate_path], cutoffs=PPI_CUTOFFS)

# Refiner that performs some cleanup of the hierarchy
refiner = HiDeFHierarchyRefiner(provenance_utils=prov_utils)

# Converter that converts CDAPS CX hierarchy into HCX format
converter = HCXFromCDAPSCXHierarchy()

# Creates hierarchy generator
hiergen = CDAPSHiDeFHierarchyGenerator(refiner=refiner,
                                       hcxconverter=converter,
                                       provenance_utils=prov_utils)

# Constructor of the object that takes all the objects created above to make a hierarchy
x = CellmapsGenerateHierarchy(outdir=hierarchy_path,
                              inputdirs=embedding_rocrate_path,
                              ppigen=ppigen,hiergen=hiergen,
                              provenance_utils=prov_utils)

# Runs the hierarchy generation
x.run()