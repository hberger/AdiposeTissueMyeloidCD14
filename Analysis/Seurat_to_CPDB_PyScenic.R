library(Seurat)
#library(sceasy)

load("Results/data_storage/2021-05-31/BCN_NAFLD_MoMF_DC_only.combined_cluster_annotated.RData")

# source: https://www.cellphonedb.org/faq-and-troubleshooting
# take raw data and normalise it
count_raw <- scData.combined@assays$RNA@counts
count_norm <- sweep(count_raw, 2, colSums(count_raw)*1/10000, "/" )
#count_norm <- apply(count_raw, 2, function(x) (x/sum(x))*10000)
#count_norm <- scData.combined@assays$RNA@data
write.table(count_norm, 'Results/data_storage/2021-05-31/BCN_MoMF_DC_only_norm_counts_CPDB.txt', sep='\t', quote=F, col.names=NA)

# generating meta file
meta_data <- scData.combined@meta.data[,c('orig_cell_id','cluster_label'), drop=F]   #####  cluster is the user’s specific cluster column

write.table(meta_data, 'Results/data_storage/2021-05-31/BCN_MoMF_DC_only_meta_CPDB.txt', sep='\t', quote=F, row.names=F)
          
