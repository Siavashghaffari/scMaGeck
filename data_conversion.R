library("Seurat")

Name<-'sublib4'
path<-'/gstore/scratch/u/ghaffars/scMaGeck'
path_dir <- file.path(path, Name)
setwd(path_dir)
data <- Read10X(data.dir = paste0("matrix_files_", Name))
metadata <- read.csv(paste0("metadata_", Name, ".csv"))
rownames(metadata) <-metadata$X
options(Seurat.object.assay.version = "v3")
seurat_hdf5 <- CreateSeuratObject(counts = data, meta.data=metadata)


saveRDS(seurat_hdf5, file = "raw_qc.RDS")
