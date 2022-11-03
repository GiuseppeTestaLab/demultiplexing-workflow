

#https://github.com/KenLauLab/dropkick/blob/master/dropkick_tutorial.ipynb tutorial

import scanpy as sc; sc.set_figure_params(color_map="viridis", frameon=False)
import dropkick as dk
import os
import pandas as pd


if not os.path.exists(snakemake.params[1]):
	os.mkdir(snakemake.params[1])

print(snakemake.params[1])
print(snakemake.input[0])
print(snakemake.input[1])
print(snakemake.params[0])


barcodeMask = pd.read_csv(snakemake.input[1]).iloc[:,0].to_list()
adata = sc.read_10x_mtx(snakemake.input[0])
adata =  adata[barcodeMask]

adata.var_names_make_unique()


dk.dropkick(adata, n_jobs=snakemake.params[0])


adata.obs.loc[adata.obs.dropkick_score < .3,"dropkick_label_03"] = "LowQuality"
adata.obs.loc[adata.obs.dropkick_score >= .3,"dropkick_label_03"] = "GoodQuality"


adata.obs["dropkick_label_03"].to_csv(snakemake.output[0], index=True, index_label = "barcode", sep = "\t")




