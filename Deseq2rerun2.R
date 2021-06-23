
library(tidyverse)
library(Rsamtools)
library(GenomicAlignments)
library(GenomicFeatures)
library(pheatmap)
library(RColorBrewer)
library(AnnotationDbi)
library(org.Hs.eg.db)
library(DESeq2)
library(ggplot2)
library(rafalib)
library(magrittr)

setwd('/data/ltt/projects/LiTing/GSE135815/results/')


# 1 building count matrix: SummarizedExperiment ----------------------------
#meta data: 
rnaseq_meta <- read_tsv('./deseq2_rerun2//rna_seq_meta.txt')
rnaseq_meta

bamfiles <- BamFileList(rnaseq_meta$bam_file, yieldSize = 2000000)
seqinfo(bamfiles)
#gtf file:
hg38_v37 <- makeTxDbFromGFF('/ref/annotation/gencode.v37.chr_patch_hapl_scaff.annotation.gtf')
exons_by_gene <- exonsBy(hg38_v37, by ='gene')

#calculate the raw counts for each sample: 
# Note that fragments will be counted only once to each gene, even if they overlap
#multiple exons of a gene which may themselves be overlapping.
#
se <- summarizeOverlaps(features=exons_by_gene, 
						reads=bamfiles,
						mode="Union",
						singleEnd=FALSE,
						ignore.strand=TRUE,
						fragments=TRUE )
# add colData info 
rowRanges(se)
colData(se) <-  DataFrame(rnaseq_meta)
colData(se)
colnames(se) <- rnaseq_meta$sample
assay(se)[1:2,]
save(se, file = './deseq2_rerun2/RNA_Seq_raw_counts_se.rda')
load('./deseq2_rerun2/RNA_Seq_raw_counts_se.rda')
rownames(se)
colnames(se)
rowRanges(se)
colData(se)
se[1:2,] %>% assay
#raw read counts for each sample:
assay(se) %>% colSums %>% divide_by(1e6)
colData(se)


# 2. generate dds dataset from se dataset with DESeq2 --------------------
#generate dds:
dds <- DESeqDataSet(se, design = ~treatment)
dds$treatment

#remove unexpressed reads:
dds <- dds[rowSums(counts(dds)) >1, ]
nrow(dds)


# add size factors:
dds <- estimateSizeFactors(dds)
sizeFactors(dds)
colSums(counts(dds))

mypar()
plot(sizeFactors(dds), colSums(counts(dds)))
abline(lm(colSums(counts(dds)) ~sizeFactors(dds) + 0))

# 3. PCA and sample clustering ------------------------------------------
# 3.1 reads normalization -----------------------------------------------

# log2(normalized_counts + 1) methods (normTransform):
log_norm_counts <- log2(counts(dds, normalized =T) + 1)
log_norm_counts %>% head
#or
log_norm <- normTransform(dds)
log_norm
log_norm %>% assay %>% head

# rlog methods:
rld <- rlog(dds, blind = F)
# vst methods:
vsd <- vst(dds, blind = F)


# 3.2 PCA and hierarchical clustering on normalized counts --------------

log_norm %>% plotPCA(intgroup = 'treatment', ntop =1000)
rld %>% plotPCA(intgroup = 'treatment', ntop =1000)
vsd %>% plotPCA(intgroup = 'treatment', ntop =1000)

#hierarchical clustering
mypar(1,2)
plot(hclust(dist(t(log_norm %>% assay))))
plot(hclust(dist(t(assay(rld)))))


# 4. differential expression analysis -------------------------------------
# 4.1 DEG analysis by DESeq() function -----------------
dds <- DESeq(dds)
dds$sizeFactor

#raw read counts and normalized read counts:
counts(dds)[1:2,]
counts(dds, normalized =T) %>% head

#differentical expression results:
#contract: a character vector with exactly three elements:
#	the name of a factor in the design formula.
#	the name of the numerator level for the fold change.
#	the name of the denominator level for the fold change.
#res <- results(dds, contrast = c('treatment', 'actd', 'ctrl')) 
#res
#mcols(res)
#summary(res)

fetch_DEG <- function(deseq_set, group_name, ctrl_group, expr_group, 
                      up_sig =F, down_sig =F, sig=F, row_name = 'ensembl_id'){
  #Usage:
  #	given deseq dataset and compared group names, return clean tibble with 3 
  #col:
  #	gene_symbol, log2FC, padj
  #parameters:
  #	sig: logical, set T to only report significant DGSs with filter condition:
  #		abs(log2FC) >=1 & padj <= 0.05
  #	up_sig: logical, set T to only report significantly up-regulated genes.
  #	down_sig: logical, set T to only reprot significantly down-regulated genes.
  
  res = results(deseq_set, contrast = c(group_name, expr_group, ctrl_group))
  tb_res = res %>% as.data.frame() %>% as_tibble(rownames = row_name)
  tb_res_clean = tb_res %>% extract(c(row_name, 'log2FoldChange', 'padj'))
  
  if(sig){
    tb_res_clean %<>% filter(abs(padj) <= 0.05 & abs(log2FoldChange) >= 1)
  }else if (up_sig){
    tb_res_clean %<>% filter(abs(padj) <= 0.05 & log2FoldChange >= 1)
  }else if (down_sig){
    tb_res_clean %<>% filter(abs(padj) <= 0.05 & log2FoldChange <= -1)
  }
  
  label <- sprintf('%s_vs_%s', expr_group, ctrl_group)
  colnames(tb_res_clean) <- c(row_name,
                              paste0('log2FC_', label),
                              paste0('padj_', label))
  
  return(tb_res_clean)
}

res <- fetch_DEG(deseq_set = dds, 
                 group_name = 'treatment', 
                 expr_group = 'patient', 
                 ctrl_group = 'ctrl')
res

res_sig <- fetch_DEG(deseq_set = dds, sig = T,
                     group_name = 'treatment', 
                     expr_group = 'patient', 
                     ctrl_group = 'ctrl')
res_sig


# 4.2 annotate genes with Symbol and Entrez_id ------------
library(AnnotationDbi)
library(org.Hs.eg.db)
res$symbol <- mapIds(org.Hs.eg.db,
                     keys= sub('\\.[0-9]*$', '', x = res$ensembl_id),
                     column="SYMBOL",
                     keytype="ENSEMBL",
                     multiVals="first")
res$entrez <- mapIds(org.Hs.eg.db,
                     keys=sub('\\.[0-9]*$', '', x = res$ensembl_id),
                     column="ENTREZID",
                     keytype="ENSEMBL",
                     multiVals="first")

is.na(res$symbol ) %>% table

res %<>% arrange(padj_patient_vs_ctrl)

# 4.3 FPKM/TPM calculation --------------------------------------------------------

#calculate by fpkm function:
gene_fpkm <- fpkm(dds, robust = F) %>% as_tibble(rownames = 'ensembl_id')

#calculate fpkm from raw reads counts:

gene_tp_length <- rowRanges(dds) %>% GenomicRanges::reduce() %>% width() %>% sum()
gene_fpkm_2 <- counts(dds) %>% as_tibble(rownames = 'ensembl_id') %>%
  mutate(gene_length = gene_tp_length) %>%
  gather(key = "sample", value = "raw_counts", 2:6) %>%
  group_by(sample) %>%
  mutate(fpkm = (raw_counts/gene_length ) / sum(raw_counts) *1e9) %>%
  dplyr::select(-raw_counts, -gene_length) %>%
  spread(sample, fpkm)

tail(gene_fpkm)
tail(gene_fpkm_2)

#calculate TPM from FPKM:
gene_tpm <- 
  gene_fpkm %>%
  gather(sample, fpkm, 2:6) %>%
  group_by(sample) %>%
  mutate(tpm = (fpkm / sum(fpkm))* 1e6 ) %>%
  dplyr::select(-fpkm) %>%
  spread(sample, tpm) 

#check:
gene_tpm %>% dplyr::select(-ensembl_id) %>% colSums()


# 4.4 output DEGs -------------------------------------------------------------
#merge deg results and gene tpm values
res_tpm <- left_join(res, gene_tpm)

res_tpm

res_tpm %<>% arrange(padj_patient_vs_ctrl)
res_tpm_up <- filter(res_tpm, padj_actd_vs_ctrl<=0.05 &  log2FC_actd_vs_ctrl >= 1)
res_tpm_down <- filter(res_tpm, padj_actd_vs_ctrl<=0.05 &  log2FC_actd_vs_ctrl <= -1)

 write_csv(res_tpm, './deseq2_rerun2/DEG_patients_vs_ctrl_5samples.csv')
# write_csv(res_tpm_up, './results/deseq2/up_DEG_ko_vs_wt_log2FC_1_pajd_0.05.csv')
# write_csv(res_tpm_down, './results/deseq2/down_DEG_ko_vs_wt_log2FC_1_pajd_0.05.csv')



# 4.5 output TPM data for GSEA analysis -----------------------------------
system('mkdir -p ./results/gsea')

#two files are needed to perform GSEA analysis:
#	expression file: .gct
#	group file: .cls
#generate gct file:
res_tpm_gct <- res_tpm %>% mutate(DESCRIPTION = 'na') %>%
  dplyr::select(symbol, DESCRIPTION, 
                hela_actd1,hela_actd2,hela_actd3,
                hela_ctrl1, hela_ctrl2, hela_ctrl3)

gct_file <- './results/gsea/gene_tpm.gct'
cat('#1.2\n', file = gct_file)
cat("27069	6\n", file = gct_file, append = T)
write_tsv(res_tpm_gct, path = gct_file, append = T, col_names = T)

#generate cls file:
cls_file <- './results/gsea/expr.cls'
cat('6	2	1\n', file = cls_file)
cat('#actd	ctrl\n', file = cls_file, append = T)
cat('actd	actd	actd	ctrl	ctrl	ctrl\n', file = cls_file, append = T)





# 5 visualization -------------------------------------------------------

# 5.1 visualization of the top20 significant DEGs -----------------------

heatmapOfTopGenes <- function(expr_df, row_anno = NA, col_anno = NA,  n=30,
                              cluster_rows =F, cluster_cols =F,
                              show_rownames =T,show_colnames =F, 
                              main_lab =NA){
  #Usage:
  #	pheatmap of values of top n genes from given expr dataset.\
  #Parameters:
  #	expr_df: 	matrix of selected gene values.
  #	row_anno: 	vector of gene names.
  #	col_anno: 	data frame of colData.
  #	n: 			number of genes.
  expr_mat = expr_df[1:n, ] %>% as.data.frame()
  if(!is.na(row_anno) ){
    rownames(expr_mat) <- row_anno[1:n]
  }
  return(	pheatmap(mat = expr_mat,
                   annotation_col =  col_anno,
                   cluster_rows = cluster_rows,
                   cluster_cols = cluster_cols,
                   show_colnames  = show_colnames,
                   show_rownames = show_rownames,
                   main = main_lab))
}


heatmapOfTopGenes(res_tpm[6:11] %>% log2(),
                  row_anno = res_tpm$symbol,
                  col_anno = colData(dds)['treatment'] %>% as.data.frame(),
                  show_colnames = T, n = 30,
                  main_lab = 'Top 30 DEGs')


#up-regulated:

res_tpm_up <- res_tpm_up %>% filter(!is.na(symbol))


heatmapOfTopGenes(res_tpm_up[6:11] %>% log2(),
                  row_anno = res_tpm_up$symbol,
                  col_anno = colData(dds)['treatment'] %>% as.data.frame(),
                  show_colnames = T, n = 30,
                  main_lab = 'Top 30 up-regulated DEGs with KANSL1 knock-out')

#down-regulated:

res_tpm_down <- res_tpm_down %>% filter(!is.na(symbol))


heatmapOfTopGenes(res_tpm_down[6:11] %>% log2(),
                  row_anno = res_tpm_down$symbol,
                  col_anno = colData(dds)['treatment'] %>% as.data.frame(),
                  show_colnames = T, n = 30,
                  main_lab = 'Top 30 down-regulated DEGs with KANSL1 knock-out')



# 5.2 visualization of the top10 GOs and KEGG pathways -------------------
# 5.2.1 david enrichment analysis ----------------------------------------
#Use david to manipulate enrichment analysis with up/down regulated gene lists

# 5.2.2 top 10 GOs -------------------------
up_gos <- read_tsv('./david/go/up-regulated-gene/up-regulated-chart-go.txt')

par(oma = c(4, 11,2,2), mar = c(2,11,2,1))
barplot( -log10(up_gos$PValue[10:1]), horiz = T, las = 2,
         sub = '-log10(p_value)',
         names.arg = up_gos$Term[10:1] %>% gsub('GO:[0-9].*~', '', x =.),
         main = 'Top 10 GO enriched by up-regulated genes')

# 5.2.3 top 10 KEGG pathways --------------------------------------------

up_kegg <- read_tsv('./david/kegg/up-regulated-genes/up-charts-kegg.txt')
barplot( -log10(up_kegg$PValue[10:1]), horiz = T, las = 2,
         sub = '-log10(p_value)',
         names.arg = up_kegg$Term[10:1] %>% gsub('hsa[0-9].*:', '', x =.),
         main = 'Top 10 kegg pathway enriched by up-regulated genes')



