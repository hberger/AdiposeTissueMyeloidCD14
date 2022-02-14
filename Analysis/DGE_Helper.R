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
  all_enrichments[["MAF_RES"]] %<-% compareCluster(geneCluster = genes_by_class, fun = "enrichMAFRES")
  
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
