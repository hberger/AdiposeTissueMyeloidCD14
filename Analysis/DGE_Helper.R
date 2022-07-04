compute_DGE <- function(data, target_suffix, control_suffix, group_parameter, min_cells_per_group = 5) {
  all_celltypes = unique(data$cluster_label)
  
  dge_results = list()
  
  cellcounts_per_group = table(data$cluster_label, data@meta.data[[group_parameter]])
  sel_celltypes = rownames(cellcounts_per_group[apply(cellcounts_per_group, 1, min)>min_cells_per_group,])
  
  for (tt in sel_celltypes) {
    tmp = FindMarkers(data , ident.1 = paste0(tt,target_suffix), ident.2 = paste0(tt,control_suffix), verbose = FALSE, assay = "RNA")
    tmp$celltype = tt
    tmp$GeneSymbol = rownames(tmp)
    dge_results[[tt]] = tmp
  }
  
  return(dge_results)
}

import_sigs <- function() {
  d1 = as.data.frame(read_excel("../Literature/Hildreth_NatImmunology_Cluster_Markers.xlsx", sheet = "Cluster_Markers_CD45PosNeg", skip = 2))
  d2 = as.data.frame(read_excel("../Literature/Hildreth_NatImmunology_Cluster_Markers.xlsx", sheet = "Cluster_Markers_Cellype_NK_ILCs", skip = 2))
  d3 = as.data.frame(read_excel("../Literature/Hildreth_NatImmunology_Cluster_Markers.xlsx", sheet = "Cluster_Markers_Celltype_Myeloi", skip = 2))
  
  d1 = subset(d1, abs(avg_logFC) > 0.5)
  d2 = subset(d2, abs(avg_logFC) > 0.5)
  d3 = subset(d3, abs(avg_logFC) > 0.5)
  
  d_all = rbind(d1, d2, d3)[, c("cluster","gene")]
  colnames(d_all) = c("gs_name","GeneSymbol")
  d_all$gs_name = paste0("Hildgreth_", d_all$gs_name)
  
  gs_to_entrez = select(org.Hs.eg.db, keys=unique(d_all$GeneSymbol), keytype="SYMBOL",columns="ENTREZID")
  gs_to_entrez = gs_to_entrez[!duplicated(gs_to_entrez$SYMBOL),]
  rownames(gs_to_entrez) = gs_to_entrez$SYMBOL

  angiogenesis_IPA=unlist(strsplit("ADAM9,ADM,AKR1B1,APOE,ATF3,AXL,C1QA,CCL2,CD28,CD36,CD63,CD81,CDK6,CDKN1A,CHD7,CRYAB,CSF1R,CTSB,CXCL12,DAAM1,DAB2,EGFL7,EGR1,EIF4A3,EMP1,EMP2,ENG,EPAS1,FGFR1,GAS6,HBEGF,HLA-DQB1,HMOX1,HRH1,HTRA1,ICAM1,IGF1,IGFBP4,IL18,JUN,KCNMA1,KLF2,LMNA,LYVE1,MERTK,MYC,NR4A1,NRP1,NRP2,PDGFB,PDGFC,PFKFB3,PLEKHG5,PLPP3,PLXDC1,PRDM1,PRKACB,RCAN1,RHOB,RND3,SAA1,SCARB1,SPRED1,STAB1,STARD13,TGFBI,THBD,TIMP3,TNFRSF25,VCAM1,VEGFA,WLS",","))
  endothelial_chemotaxis_top = unlist(strsplit("LGMN,NRP1,HSPB1,VEGFA,FGFR1,NR4A1,PLEKHG5",","))
  d_all = rbind(d_all, data.frame(gs_name = "Angiogenesis_IPA", GeneSymbol=angiogenesis_IPA))
  d_all = rbind(d_all, data.frame(gs_name = "Endothelial_chemotaxis_top", GeneSymbol=endothelial_chemotaxis_top))
  
  d_all$entrez_gene = gs_to_entrez[d_all$GeneSymbol,"ENTREZID"]
  d_all[["GeneSymbol"]] <- NULL
  
  return(d_all)
  #write.table(hildgreth_ref_top_lists_gs, file="hildgreth.txt", row.names=F, sep="\t", quote=F)
  
}

plot_top_genes <- function(data, dge_results_all, title) {
  for (n in sort(unique(data$cluster_label))) {
    tmp = subset(dge_results_all, celltype==n)
    if(nrow(tmp)<1) next
    tmp_up = subset(tmp, avg_log2FC > 0)
    tmp_up$direction = "up"
    tmp_down = subset(tmp, avg_log2FC < 0)
    tmp_down$direction = "down"
    
    final_top_genes = rbind(head(tmp_up,10), head(tmp_down,10))
    final_top_genes = final_top_genes[order(final_top_genes$avg_log2FC),]
    final_top_genes$gene = factor(final_top_genes$GeneSymbol, levels = final_top_genes$GeneSymbol)
    p = ggplot(final_top_genes, aes(y=avg_log2FC, x=gene, fill=direction)) + geom_bar(stat="identity") + ggtitle(paste0(n, ", ", title)) + coord_flip() + theme(axis.text.y =  element_text(size=12) )
    print(p)
    
  }
}

enrichGSEA <- function(genes, collection, organism) {
  m_t2g <- msigdbr(species = organism, category = collection)[, c("gs_name", "entrez_gene")]
  m_t2g$gs_name = ifelse(nchar(m_t2g$gs_name)>80, paste0(substr(m_t2g$gs_name,1,80),"~"), m_t2g$gs_name)
  enricher(genes, TERM2GENE=m_t2g) 
}

enrichMAFRES <- function(genes) {
  m_t2g <- maf_ref_top_lists_gs
  m_t2g$gs_name = ifelse(nchar(m_t2g$gs_name)>80, paste0(substr(m_t2g$gs_name,1,80),"~"), m_t2g$gs_name)
  enricher(genes, TERM2GENE=m_t2g) 
}

enrichSignatures <- function(genes) {
  m_t2g <- signatures_gs
  m_t2g$gs_name = ifelse(nchar(m_t2g$gs_name)>80, paste0(substr(m_t2g$gs_name,1,80),"~"), m_t2g$gs_name)
  enricher(genes, TERM2GENE=m_t2g) 
}

compute_functional_enrichment <- function(data, curr_dge_tab, title) {
  require(clusterProfiler)
  require(org.Hs.eg.db)
  require(msigdbr)
  require(ReactomePA)
    
  curr_dge_tab$direction = with(curr_dge_tab, ifelse(avg_log2FC > 0, "Up", ifelse(avg_log2FC < 0, "Down","nc")) )
  curr_dge_tab$gene_group = with(curr_dge_tab, paste0(celltype,"_",direction))
  
  tmp = AnnotationDbi::select(org.Hs.eg.db, keys=rownames(sel_dds@assays$RNA), keytype = "SYMBOL", columns="ENTREZID")

  gs2e = tapply(tmp$ENTREZID, tmp$SYMBOL, paste, collapse=",")
  e2gs = tapply(tmp$SYMBOL, tmp$ENTREZID, paste, collapse=",")
  
  genes_by_class = split(gs2e[curr_dge_tab$GeneSymbol], factor(curr_dge_tab$gene_group))
  genes_by_class = lapply(genes_by_class, function(x) unique(x[!is.na(x)]))
  
  
  all_enrichments = new.env() # futures can only be assigned to an environment, not a list

  all_enrichments[["GO_BP"]] %<-% compareCluster(geneCluster = genes_by_class, fun = "enrichGO", OrgDb="org.Hs.eg.db", ont="BP")
  #all_enrichments[["KEGG"]] %<-% compareCluster(geneCluster = genes_by_class, fun = "enrichKEGG", organism = "hsa", use_internal_data=T)
  all_enrichments[["REACTOME"]] %<-% compareCluster(geneCluster = genes_by_class, fun = "ReactomePA::enrichPathway", organism = "human")
  all_enrichments[["GSEA_H"]] %<-% compareCluster(geneCluster = genes_by_class, fun = "enrichGSEA", organism = "Homo sapiens", collection="H")
  all_enrichments[["GSEA_C2"]] %<-% compareCluster(geneCluster = genes_by_class, fun = "enrichGSEA", organism = "Homo sapiens", collection="C2")
  all_enrichments[["GSEA_C3"]] %<-% compareCluster(geneCluster = genes_by_class, fun = "enrichGSEA", organism = "Homo sapiens", collection="C3")
  all_enrichments[["GSEA_C7"]] %<-% compareCluster(geneCluster = genes_by_class, fun = "enrichGSEA", organism = "Homo sapiens", collection="C7")
  #all_enrichments[["MAF_RES"]] %<-% compareCluster(geneCluster = genes_by_class, fun = "enrichMAFRES")
  all_enrichments[["Signatures"]] %<-% compareCluster(geneCluster = genes_by_class, fun = "enrichSignatures")
  
  all_enrichments = as.list(all_enrichments)
    
  # write to table
  for (n in names(all_enrichments)) {
    tmp = all_enrichments[[n]]
    tmp@compareClusterResult$database = n
    all_enrichments[[n]] = tmp
  }
  
  HM_cluster_enrichment = do.call(rbind, lapply(all_enrichments, function(x) x@compareClusterResult) )
  
  HM_cluster_enrichment$GeneSymbols = Map(function(x) {sapply(strsplit(x,"/"), function(y) sort(paste(e2gs[y], collapse=",")))}, HM_cluster_enrichment$geneID )
  cc = colnames(HM_cluster_enrichment)
  sel_cols = c("Cluster","database")
  sel_cols = append(sel_cols, cc[!cc %in% sel_cols & !cc =="geneID"])
  fwrite(HM_cluster_enrichment[, sel_cols], file=file.path(result_folder,paste0("DGE_",title,"_enrichments_all.txt")), sep="\t", quote=F)

  return(all_enrichments)
}
