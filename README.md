# VisiumH5AD2Seurat
Convert spatial transcriptomics H5AD file into 10X original data format

# I hope this can be helpful for those who encounter similar difficulties!

## Purpose

Both R environment and Python environment are important platforms for downstream analysis of high-dimensional sequencing data, with popular packages including Seurat in R and scanpy in Python. There are already many interactive solutions for converting Anndata H5AD format from single-cell sequencing into Seurat, such as [zellkonverter](https://github.com/theislab/zellkonverter). However, for Visium spatial transcriptomics, there is currently a lack of effective conversion solutions. Maintainers from the Satija Lab have expressed [disdain](https://github.com/satijalab/seurat/issues/8191) for this issue, which is truly absurd. I can see efforts from Scverse team to be compatible with the Seurat environment, like zellkonverter, but Satija Lab seems lazy about this matter, which is really disappointing.

Original h5ad data set is from [locationLung](https://cellgeni.cog.sanger.ac.uk/5-locations-lung/Cell2location_outputs.zip). The raw data has been processed via scanpy and following cell2location, but I intend to use Seurat and [semla](https://ludvigla.github.io/semla/) for subsequent analysis.

