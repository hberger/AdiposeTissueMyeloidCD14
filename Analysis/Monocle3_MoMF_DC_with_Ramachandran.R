library(SeuratWrappers)
library(Seurat)
library(monocle3)

result_folder = paste("./Results",format(Sys.time(), "%Y-%m-%d"),sep="/")
if (!file.exists(result_folder)) dir.create(result_folder, recursive=T)

load("./Results/data_storage/2022-01-31/scData_NAFLD_Ramachandran_complete_MoMF.DC.cluster_annotated.RData")

scData_cd = as.cell_data_set(scData_MoMF.DC)
scData_cd = cluster_cells(scData_cd, reduction_method="UMAP")
scData_cd = learn_graph(scData_cd, use_partition=TRUE)

root_cells = WhichCells(scData_MoMF.DC, expression = (seurat_clusters==1)) # cla Mo

scData_cd = order_cells(scData_cd, reduction_method = "UMAP", root_cells = root_cells)

pdf(file=file.path(result_folder, "Monocle3_MoMF_DC_with_Ramachandran.pdf"), height=6, width=6)
DimPlot(scData_MoMF.DC, label = T, label.size = 2) + NoLegend()
plot_cells(cds = scData_cd, color_cells_by = "pseudotime", show_trajectory_graph = TRUE)
dev.off()

