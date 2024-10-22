###### read velocyto loom files and perform scvelo analysis ########

import scvelo as scv
import anndata

# Load the preprocessed anndata object containing metadata
adata = anndata.read_h5ad('SSRI_anndata.h5ad')

# Load the velocyto loom files for both conditions
loom_file_condition1 = 'P29_sertraline.loom'
loom_file_condition2 = 'P29_ssri_control.loom'

# Read the loom files into AnnData objects
adata_condition1 = scv.read(loom_file_condition1, cache=True)
adata_condition2 = scv.read(loom_file_condition2, cache=True)

# Merge the metadata from the preprocessed anndata object into the loom data
adata_condition1.obs = adata.obs
adata_condition2.obs = adata.obs

# Concatenate the two conditions into a single AnnData object
adata_combined = adata_condition1.concatenate(adata_condition2, batch_key='condition', batch_categories=['condition1', 'condition2'])

# Perform scvelo analysis
scv.pp.filter_and_normalize(adata_combined)
scv.pp.moments(adata_combined)
scv.tl.velocity(adata_combined)
scv.tl.velocity_graph(adata_combined)


# Plot the velocity embedding
scv.pl.velocity_embedding_stream(adata_combined, basis='umap', color='condition')
