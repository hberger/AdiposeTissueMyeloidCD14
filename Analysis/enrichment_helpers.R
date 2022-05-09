require(org.Hs.eg.db)
require(data.table)

set_gene_list <- function(genes, ref_folder) {
  g2e <<- AnnotationDbi::select(org.Hs.eg.db, keys=genes, keytype = "SYMBOL",columns = "ENTREZID")
  g2e <<- g2e[!duplicated(g2e$SYMBOL),]
  rownames(g2e) <<- g2e$SYMBOL
  
  homologene_human_mouse = fread(file.path(ref_folder, "HomologousGenes/HomoloGene/build68/homologene_9606_10090.txt"))
  all_entrez_ids = data.frame(EntrezID=unique(g2e$ENTREZID))
  all_symbols = data.frame(GeneSymbol=unique(g2e$SYMBOL))
  a1 = merge(all_entrez_ids, homologene_human_mouse[,c("entrez_1","symbol_1","entrez_2","symbol_2"), with=F], by.x="EntrezID", by.y="entrez_1", all.x=T, sort=F)
  
  human_to_mouse <<- a1
  human_to_mouse <<- subset(human_to_mouse, !is.na(entrez_2) & !is.na(symbol_2) & !symbol_2 == "" & !is.na(EntrezID))
  rownames(human_to_mouse) <<- human_to_mouse$EntrezID
  
  all_symbols = data.frame(GeneSymbol=unique(g2e$SYMBOL))
  a1 = merge(all_symbols, homologene_human_mouse[,c("entrez_1","symbol_1","symbol_2"), with=F], by.x="GeneSymbol", by.y="symbol_1", all.x=T, sort=F)
  
  mouse_to_human <<- a1
  mouse_to_human <<- subset(mouse_to_human, !is.na(GeneSymbol) & !is.na(symbol_2) & !symbol_2 == "")
  rownames(mouse_to_human) <<- mouse_to_human$symbol_2
  
  
}


load_sigs <- function() {

  d1 = as.data.frame(read_excel("../Literature/Hildreth_NatImmunology_Cluster_Markers.xlsx", sheet = "Cluster_Markers_CD45PosNeg", skip = 2))
  d2 = as.data.frame(read_excel("../Literature/Hildreth_NatImmunology_Cluster_Markers.xlsx", sheet = "Cluster_Markers_Cellype_NK_ILCs", skip = 2))
  d3 = as.data.frame(read_excel("../Literature/Hildreth_NatImmunology_Cluster_Markers.xlsx", sheet = "Cluster_Markers_Celltype_Myeloi", skip = 2))
  
  d1 = subset(d1, abs(avg_logFC) > 0.5)
  d2 = subset(d2, abs(avg_logFC) > 0.5)
  d3 = subset(d3, abs(avg_logFC) > 0.5)
  
  d_all = rbind(d1, d2, d3)[, c("cluster","gene")]
  colnames(d_all) = c("gs_name","GeneSymbol")
  
  # gs_to_entrez = select(org.Hs.eg.db, keys=unique(d_all$GeneSymbol), keytype="SYMBOL",columns="ENTREZID")
  # gs_to_entrez = gs_to_entrez[!duplicated(gs_to_entrez$SYMBOL),]
  # rownames(gs_to_entrez) = gs_to_entrez$SYMBOL
  
  d_all$entrez_gene = g2e[d_all$GeneSymbol,"ENTREZID"]
  d_all[["GeneSymbol"]] <- NULL
  d_all$gs_name = paste0("Hildgreth_", d_all$gs_name)
  
  #hildgreth_ref_top_lists_gs = d_all
  #write.table(hildgreth_ref_top_lists_gs, file="hildgreth.txt", row.names=F, sep="\t", quote=F)

  # sc_signature_tab = data.frame(gs_name = rep(names(sc_signatures), times=unlist(lapply(sc_signatures, length)) ), entrez_gene = g2e[unlist(sc_signatures),"ENTREZID"], stringsAsFactors = F )
  
  sc_signature_tab = d_all
  sc_signature_tab = subset(sc_signature_tab, !is.na(entrez_gene))
  
  return(sc_signature_tab)
}

load_maf_ressource <- function() {
  # differentially expressed genes from macrophage stimulations in GSE47189
  maf_resource_env = new.env()
  load("../../DataSets/GSE47189_MacrophageResource/DGE_results_Macrophages.Rdata", envir = maf_resource_env)
  
  all_maf_dge = get("all_results", maf_resource_env)
  

  maf_res_top_lists = list()
  
  # select differentially expressed genes and convert to mouse Entrez IDs
  for (n in names(all_maf_dge)) {
    n_new = gsub('\UFFB5',"^", n)
    tmp = subset(all_maf_dge[[n]], logFC > 1 & adj.P.Val < 0.05 & !is.na(ENTREZID))
    if(nrow(tmp)>0) maf_res_top_lists[[paste0(n,"_UP")]] = data.frame(gs_name = paste0(n_new,"_UP"), entrez_gene = tmp$ENTREZID, stringsAsFactors = F)
    tmp = subset(all_maf_dge[[n]], logFC < -1 & adj.P.Val < 0.05 & !is.na(ENTREZID))
    if(nrow(tmp)>0) maf_res_top_lists[[paste0(n,"_DOWN")]] = data.frame(gs_name = paste0(n_new,"_DOWN"), entrez_gene = tmp$ENTREZID, stringsAsFactors = F)
  }
  
  maf_ref_top_lists_gs = do.call(rbind, maf_res_top_lists)
  rownames(maf_ref_top_lists_gs) <- NULL
  
  rm(all_maf_dge)
  gc()

  return(maf_ref_top_lists_gs)  
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

enrichScSigs <- function(genes) {
  m_t2g <- sc_signature_tab
  m_t2g$gs_name = ifelse(nchar(m_t2g$gs_name)>80, paste0(substr(m_t2g$gs_name,1,80),"~"), m_t2g$gs_name)
  enricher(genes, TERM2GENE=m_t2g)
}


compute_functional_enrichment <- function(sc_data, curr_dge_tab, title) {
  
  curr_dge_tab$direction = with(curr_dge_tab, ifelse(avg_log2FC > 0, "Up", ifelse(avg_log2FC < 0, "Down","nc")) )
  curr_dge_tab$gene_group = with(curr_dge_tab, paste0(celltype,"_",direction))
  
  tmp = AnnotationDbi::select(org.Hs.eg.db, keys=rownames(sc_data@assays$RNA), keytype = "SYMBOL", columns="ENTREZID")
  
  gs2e = tapply(tmp$ENTREZID, tmp$SYMBOL, paste, collapse=",")
  e2gs = tapply(tmp$SYMBOL, tmp$ENTREZID, paste, collapse=",")
  
  genes_by_class = split(gs2e[curr_dge_tab$GeneSymbol], factor(curr_dge_tab$gene_group))
  genes_by_class = lapply(genes_by_class, function(x) unique(x[!is.na(x)]))
  
  
  all_enrichments = new.env() # futures can only be assigned to an environment, not a list
  
  # future has problems identifying the enrich... functions from preexisting packages as global variabales (at least on Windows) - lets define the full path for them
  all_enrichments[["GO_BP"]] %<-% compareCluster(geneCluster = genes_by_class, fun = enrichGO, OrgDb="org.Hs.eg.db", ont="BP")
  all_enrichments[["KEGG"]] %<-% compareCluster(geneCluster = genes_by_class, fun = enrichKEGG, organism = "hsa", use_internal_data=T)
  all_enrichments[["REACTOME"]] %<-% compareCluster(geneCluster = genes_by_class, fun = enrichPathway, organism = "human")
  all_enrichments[["GSEA_H"]] %<-% compareCluster(geneCluster = genes_by_class, fun = enrichGSEA, organism = "Homo sapiens", collection="H")
  all_enrichments[["GSEA_C2"]] %<-% compareCluster(geneCluster = genes_by_class, fun = enrichGSEA, organism = "Homo sapiens", collection="C2")
  all_enrichments[["GSEA_C3"]] %<-% compareCluster(geneCluster = genes_by_class, fun = enrichGSEA, organism = "Homo sapiens", collection="C3")
  all_enrichments[["GSEA_C8"]] %<-% compareCluster(geneCluster = genes_by_class, fun = enrichGSEA, organism = "Homo sapiens", collection="C8")
  #all_enrichments[["MAF_RES"]] %<-% compareCluster(geneCluster = genes_by_class, fun = enrichMAFRES)
  all_enrichments[["scSignatures_Lit"]] %<-% compareCluster(geneCluster = genes_by_class, fun = enrichScSigs)  
  
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
  
  tmp <- group_split(HM_cluster_enrichment %>% group_by(Cluster))
  tmp2 <- lapply(tmp, function(x) x[order(x$database, x$pvalue), c("Cluster","database","Description","GeneRatio","BgRatio","pvalue","p.adjust","qvalue","Count","GeneSymbols")] )
  tmp2 <- lapply(tmp2, function(x) {x$GeneSymbols = unlist(lapply(x$GeneSymbols, paste, collapse=", ")) ; return(x) } )
  names(tmp2) = unlist(lapply(tmp2, function(x) unique(make.names(as.character(x$Cluster)))))
  
  write_xlsx(tmp2, path=file.path(result_folder,paste0("DGE_",title,"_enrichments_all.xlsx")))
  
  
  return(all_enrichments)
}


