import scanpy as sc
import matplotlib.pyplot as plt
import pandas as pd
import json
import os
import h5py
from cell2location.utils import select_slide
adata=sc.read_h5ad('./Cell2location_outputs/sp.h5ad')
del(adata.obsm['means_cell_abundance_w_sf'])
del(adata.obsm['q05_cell_abundance_w_sf'])
del(adata.obsm['q95_cell_abundance_w_sf'])
del(adata.obsm['stds_cell_abundance_w_sf'])


adata
# Write 10X h5 function cited from https://github.com/scverse/anndata/issues/595
def write_10X_h5(adata, file):
    if '.h5' not in file: file = f'{file}.h5'
    def int_max(x):
        return int(max(np.floor(len(str(int(max(x)))) / 4), 1) * 4)
    def str_max(x):
        return max([len(i) for i in x])

    w = h5py.File(file, 'w')
    grp = w.create_group("matrix")
    grp.create_dataset("barcodes", data=np.array(adata.obs_names, dtype=f'|S{str_max(adata.obs_names)}'))
    grp.create_dataset("data", data=np.array(adata.X.data, dtype=f'<i{int_max(adata.X.data)}'))
    ftrs = grp.create_group("features")
    # this group will lack the following keys:
    # '_all_tag_keys', 'feature_type', 'genome', 'id', 'name', 'pattern', 'read', 'sequence'
    ftrs.create_dataset("feature_type", data=np.array(adata.var.feature_types, dtype=f'|S{str_max(adata.var.feature_types)}'))
    ftrs.create_dataset("genome", data=np.array(adata.var.genome, dtype=f'|S{str_max(adata.var.genome)}'))
    ftrs.create_dataset("id", data=np.array(adata.var.SYMBOL, dtype=f'|S{str_max(adata.var.SYMBOL)}'))
    ftrs.create_dataset("name", data=np.array(adata.var.SYMBOL, dtype=f'|S{str_max(adata.var.SYMBOL)}'))
    grp.create_dataset("indices", data=np.array(adata.X.indices, dtype=f'<i{int_max(adata.X.indices)}'))
    grp.create_dataset("indptr", data=np.array(adata.X.indptr, dtype=f'<i{int_max(adata.X.indptr)}'))
    grp.create_dataset("shape", data=np.array(list(adata.X.shape)[::-1], dtype=f'<i{int_max(adata.X.shape)}'))

def writeVisium(adata_vis,spName):
  slide = select_slide(adata_vis,spName)
  slide=slide.raw.to_adata()
  new_dir = spName
  os.makedirs(new_dir, exist_ok=True)
  os.chdir(new_dir)
  write_10X_h5(slide,"filtered_feature_bc_matrix.h5")
  os.makedirs("spatial", exist_ok=True)
  os.chdir("spatial")
  slide = select_slide(slide, spName)
  plt.imsave('tissue_hires_image.png', slide.uns['spatial'][spName]['images']['hires'])
  plt.imsave('tissue_lowres_image.png', slide.uns['spatial'][spName]['images']['lowres'])
  with open("scalefactors_json.json", 'w') as file:
      json.dump(slide.uns['spatial'][spName]['scalefactors'], file)
  df1_selected = slide.obs.iloc[:, :3]
  Obsm=pd.DataFrame(slide.obsm['spatial'])
  Obsm[[1, 0]] = Obsm[[0, 1]]
  Obsm.index=slide.obs.index
  df_horizontal = pd.concat([df1_selected,Obsm ], axis=1)
  df_horizontal.to_csv('tissue_positions_list.csv',header=False)
  os.chdir("..")
  os.chdir("..")


writeVisium(adata,"WSA_LngSP9258464")
