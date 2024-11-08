library(scMAGeCK)
library(Seurat)

### Please Give this Project a name and assign the input and output directory
Name<-'sublib4'

path<-'/gstore/scratch/u/ghaffars/scMaGeck'  # input Directory

out_dir <- file.path("/gstore/project/crc_recursion_2/scMaGeck/", Name) # output Directory
##################################


# Set the path directories
path_dir <- file.path(path, Name)
setwd(path_dir)

# Create the output directory
if (file.exists(out_dir)){
} else {
  dir.create(out_dir)
}

# set the BARCODE and RDS file path 
BARCODE = paste0("barcode_rec_",Name,".txt")

bc_frame=read.table(BARCODE,header = T,as.is = T)

## RDS can be a Seurat object or local RDS file path that contains the scRNA-seq dataset
RDS = "raw_qc.RDS"
rds_object=readRDS(RDS)

## Run all the genes: gene by gene
gene_list<-grep('NTC',unique(rds_object$gene_symbol),value = T, invert = T)
for (gene_list1 in gene_list){
    eff_object <- scmageck_eff_estimate(rds_object, bc_frame, perturb_gene=gene_list1, 
                                     non_target_ctrl = 'NTC')

    eff_estimat=eff_object$eff_matrix
    rds_subset=eff_object$rds
    eff_traget=eff_object$target_gene_search_result
    
    write.csv(eff_traget[[gene_list1]]$deframe, file=paste0(out_dir,"/DF_", gene_list1,".csv"))
    df<- data.frame(eff_estimat)
    write.csv(df, file=paste0(out_dir,"/PS_", gene_list1,".csv"))}


