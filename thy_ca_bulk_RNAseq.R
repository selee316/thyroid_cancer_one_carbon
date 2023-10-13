# library
if(T){
  library(tidyverse)
  library(readxl)
  library(ggpubr)
  library(ggrepel)
  library(GSVA)
  library(corrplot)
  library(DESeq2)
  library(fgsea)
  library(circlize)
  library(piano)
  library(GSA)
  library(ggsci)
  library(scales)
  library(mgsub)
  library(cowplot)
}

# design setting
theme_se = theme_bw() + theme(axis.text = element_text(size=18), axis.title = element_text(size=18), 
                              panel.border = element_blank(), axis.line = element_line(), plot.title = element_text(size=20, face="bold"))

if(F){
  ec_dt <- read_tsv('/path/to/thyroid_944s_merged_RSEM_expected_count.txt')
  ec_dt[1:5,1:5]
  ec_dt <- ec_dt[,c('gene_id', total_ids)]
  dim(ec_dt) #60669 828
  
  #load biomart data
  bm_dt <- read_tsv('/path/to/biomart_human_geneID_hgncsymbol_ensemblgenename_genetype.txt')
  bm_dt <- bm_dt %>% dplyr::rename(gene_id = `Gene stable ID`, gene_type = `Gene type`, gene_name = `Gene name`) %>% select(gene_id, gene_type, gene_name) %>% unique()
  bm_dt %>% group_by(gene_id) %>% dplyr::count() %>% arrange(desc(n)) #check duplicated ids
  
  #expression profile 
  dim(ec_dt) #60669 828
  ec_dt <- ec_dt[,c("gene_id", m_meta_dt$sample_id)]
  dim(ec_dt) #60669 828
  ec_dt <- ec_dt %>% separate(gene_id, c('gene_id','gene_subid'), sep='\\.')
  ec_dt %>% group_by(gene_id) %>% dplyr::count() %>% filter(n>1) #check duplicated ids 45rows
  ec_dt %>% filter(gene_id == 'ENSG00000002586')
  ec_dt$rowSum <- ec_dt %>% select(starts_with('TH_')) %>% rowSums()
  ec_dt$rowMean <- ec_dt %>% select(starts_with('TH_')) %>% rowMeans()
  ec_dt <- ec_dt %>% filter(rowSum > 0) #remove rows with all 0s
  ec_dt %>% group_by(gene_id) %>% dplyr::count() %>% filter(n>1) #recheck duplicated ids 0
  ec_dt <- left_join(ec_dt, bm_dt)
  ec_dt <- ec_dt %>% filter(gene_type=='protein_coding') #filter in protein-coding only
  nrow(ec_dt) #before duplicated row removal 19842
  dup_list <- ec_dt %>% group_by(gene_name) %>% dplyr::count() %>% filter(n>1) %>% pull(gene_name)
  length(dup_list) #No. of duplicated gene names 13
  inc_tbl <- tibble()
  for (t_gene in dup_list){
    inc_tbl <- bind_rows(inc_tbl, ec_dt %>% filter(gene_name == t_gene) %>% arrange(desc(rowMean)) %>% head(n=1))
  }
  ec_dt <- bind_rows(ec_dt %>% filter(!gene_name %in% dup_list), inc_tbl)
  nrow(ec_dt) #after duplicated row removal 19829
  if(F){
    ec_dt %>% saveRDS('/path/to/ec_dt.rds')
  }
  
  ec_mx <- ec_dt %>% select(gene_name, starts_with('TH_')) %>% column_to_rownames('gene_name') %>% as.matrix() 
  ec_mx <- round(ec_mx)  #rounding
}

#run DESeq2

if(F){
  
  col_dt <- m_meta_dt %>% select(sample_id, donor_group, TorN) %>% as.data.frame() %>% column_to_rownames('sample_id')
  all(colnames(ec_mx) == rownames(col_dt))
  col_dt <- col_dt[colnames(ec_mx),]
  all(colnames(ec_mx) == rownames(col_dt))
  
  dds <- DESeqDataSetFromMatrix(countData = ec_mx,
                                colData = col_dt,
                                design = ~ TorN)
  un_dds <- DESeq(dds)
  un_fpm <- fpm(un_dds, robust = T)
  un_fpm <- un_fpm[rowSums(un_fpm) > 0,]
}
if(F){
  saveRDS(un_dds, '/path/to/un_dds.rds')
  saveRDS(un_fpm, '/path/to/un_fpm.rds')
}

#batch effect corrrction using ComBat_seq
if(F){
  m_meta_dt <- read_tsv('/path/to/210313_Thyroid_ca_827s_meta_tbl.tsv', col_types = cols(Age = 'n'))
  # add batch on meta table
  m_meta_dt <- m_meta_dt %>% mutate(batch = ifelse(library_type == 'str' & tissue_buffer == 'PAXgene', 'str_PAX',
                                                   ifelse(library_type == 'str', 'str_others',
                                                          ifelse(library_type == 'acc' & tissue_buffer %in% c('PAXgene','RNAlater'), 'acc_PAX_later','acc_others'))))
  
  batch <- m_meta_dt %>% select(sample_id, batch) %>% column_to_rownames('sample_id') %>% .[colnames(ec_mx),]
  ec_bc_mx <- sva::ComBat_seq(ec_mx, batch=batch, group=NULL)
  if(F){
    saveRDS(ec_bc_mx, '/path/to/ec_bc_mx.rds')
  }
}

#Run DESeq2 with bc_ec_mx, T vs N
if(F){
  ec_bc_mx <- readRDS('/path/to/ec_bc_mx.rds')
  col_dt <- m_meta_dt %>% select(sample_id, tissue_group, TorN) %>% column_to_rownames('sample_id')
  all(colnames(ec_bc_mx) == rownames(col_dt))
  col_dt <- col_dt[colnames(ec_bc_mx),]
  all(colnames(ec_bc_mx) == rownames(col_dt))
  
  bc_dds <- DESeqDataSetFromMatrix(countData = ec_bc_mx,
                                   colData = col_dt,
                                   design = ~ TorN)
  library("BiocParallel")
  bc_dds <- DESeq(bc_dds, parallel = T, BPPARAM=MulticoreParam(8))
  bc_fpm <- fpm(bc_dds, robust = T)
  bc_fpm <- bc_fpm[rowSums(bc_fpm) > 0,]
  
  if(F){
    saveRDS(bc_dds, '/path/to/bc_dds.rds')
    saveRDS(bc_fpm, '/path/to/bc_fpm.rds')
  }
}

#make tpm table
tpm <- function(counts,len) {
  x <- counts/len
  return(t(t(x)*1e6/colSums(x)))
}

bc_dds <- readRDS('/path/to/bc_dds.rds')
ct_dt <- counts(bc_dds, normalized= T)
dim(ct_dt)
el_dt <- read_tsv('/path/to/merged_RSEM_effective_length.txt')
el_dt$mean_EL <- rowMeans(el_dt[2:ncol(el_dt)])
el_dt <- el_dt %>% select(gene_id, mean_EL) %>% separate(gene_id, c('gene_id', 'gene_subid'), sep='\\.')
ec_dt <- readRDS('/path/to/ec_dt.rds')
m_el_dt <- left_join(ec_dt %>% select(gene_id, gene_subid, gene_name), el_dt %>% select(gene_id, gene_subid, mean_EL)%>% unique())
m_el_dt %>% filter(gene_name == 'MTHFS')
el_mx <-  m_el_dt %>% select(gene_name, mean_EL) %>% column_to_rownames('gene_name') %>% as.matrix()
el_v <- el_mx[rownames(ct_dt),]
#remove zero length genes
zero_lengths <- names(el_v[el_v==0])
ct_dt <- ct_dt[!rownames(ct_dt) %in% zero_lengths,] 
el_v <- el_v[!names(el_v) %in% zero_lengths]
all(rownames(ct_dt)==names(el_v))
tpm_mx <- tpm(ct_dt, el_v)
l_tpm_mx <- log10(tpm_mx+0.01)


# load gene sets
c5_list <- GSEABase::getGmt('/path/to/c5.all.v7.5.1.symbols.gmt') %>% GSEABase::geneIds()
c2_list <- GSEABase::getGmt('/path/to/c2.all.v7.5.1.symbols.gmt') %>% GSEABase::geneIds()

gobp_names <- names(c5_list)[grepl('GOBP',names(c5_list))]
gobp_list <- c5_list[gobp_names]
kegg_names <- names(c2_list)[grepl('KEGG',names(c2_list))]
kegg_list <- c2_list[kegg_names]

oc_gene_list <- c('PHGDH','PSAT1','PSPH','SHMT1','SHMT2','MTHFD1','MTHFD2','MTHFD1L','MTHFD2L','TYMS', 'MTFMT', 'DHFR',
                  'MTHFR','ALDH1L1','ALDH1L2','SFXN1','SFXN2','SFXN3', 'SLC1A5', 'FOLH1', 'SLC1A4', 'FOLR1', 'FOLR2')
setdiff(oc_gene_list, rownames(l_tpm_mx))  

tds_genes = c('DIO1','DIO2','DUOX1','DUOX2','FOXE1','GLIS3','NKX2-1','PAX8','SLC26A4','SLC5A5','SLC5A8','TG','THRA','THRB','TPO','TSHR')

#TDS using gsva
if(F){
  library(GSVA)
  tds_list= list(TDSgenes = tds_genes)
  gsva_res <- gsva(l_tpm_mx, tds_list)
  gsva_res %>% saveRDS('/path/to/gsva_tds.rds')
}
tds_dt <- tibble(sample_id=names(gsva_res[1,]),tds=gsva_res[1,])


# run DESeq for pca plot
ec_bc_mx <- readRDS('/path/to/ec_bc_mx.rds')
tumor_meta_dt <- read_excel('G/path/to/primaryT_TDS_group_forSPSS.xlsx')

col_dt <- tumor_meta_dt %>% select(sample_id, tds_0_low_high, tissue_group2) %>% column_to_rownames('sample_id')
dim(col_dt) #369 2

all(colnames(ec_bc_mx) == rownames(col_dt))
s_mx <- ec_bc_mx[,rownames(col_dt)]
all(colnames(s_mx) == rownames(col_dt))
dim(s_mx)

dds <- DESeqDataSetFromMatrix(countData = s_mx,
                              colData = col_dt,
                              design = ~ tds_0_low_high)
dds <- DESeq(dds)
un_fpm <- fpm(dds, robust = T)
un_fpm <- un_fpm[rowSums(un_fpm) > 0,]

# pca plot
un_fpm.pca <- prcomp(t(un_fpm), center = T, scale = T)
str(un_fpm.pca)
tmp_dt <- un_fpm.pca$x %>% as.data.frame() %>% rownames_to_column('sample_id') %>% as_tibble()
tmp_dt <- left_join(tmp_dt, meta_dt %>% select(sample_id, tissue_group2, tds_0_low_high, tds))

show_col(pal_d3("category10")(10))
cancer_type_pal <- pal_d3("category10")(10)[c(3,2,6)]
names(cancer_type_pal) <- c('PTC_primary','ATC_primary','PD_primary')

ggplot(tmp_dt, aes(x=PC1, y=PC2, color=tds, shape=tissue_group2))+
  geom_point(size=4)+
  scale_color_viridis_c()+
  scale_shape_manual(values = c(17, 16, 0))+
  theme_se


# normal vs tumor - GSEA
# DEG
meta_dt <- read_tsv('/path/to/normal_PTC_PD_ATC_clinical_dt.tsv')
ec_bc_mx <- readRDS('/path/to/ec_bc_mx.rds')

group <- meta_dt %>% mutate(group = ifelse(tissue_group2 == 'normal_thyroid', 'Normal', 'Tumor'))

col_dt <- group %>%  select(sample_id, group, tissue_group2) %>% column_to_rownames('sample_id')
dim(col_dt) 

all(colnames(ec_bc_mx) == rownames(col_dt))
s_mx <- ec_bc_mx[,rownames(col_dt)]
all(colnames(s_mx) == rownames(col_dt))
dim(s_mx)

bc_dds <- DESeqDataSetFromMatrix(countData = s_mx,
                                 colData = col_dt,
                                 design = ~ group)
bc_dds <- DESeq(bc_dds)

res <- results(bc_dds, contrast = c('group','Tumor','Normal'))

mcols(res, use.names=TRUE)
m_res <- res %>% as.data.frame() %>% rownames_to_column('gene_name') %>% as_tibble() %>% arrange(desc(log2FoldChange))

tmp_dt <- m_res %>% mutate(nlpval = -log10(padj))
tmp_dt$nlpval[is.finite(tmp_dt$nlpval) == F] <- 300
tmp_dt <- tmp_dt %>% mutate(dist= sqrt(log2FoldChange^2 + nlpval ^2))
tmp_dt <- tmp_dt %>% mutate(dist = ifelse(log2FoldChange < 0, (-1)*dist, dist))

dds_stats <- tmp_dt$dist
names(dds_stats) <- tmp_dt$gene_name

gsea_res <- fgsea(pathways = kegg_list, stats = dds_stats, maxSize=1000)

gsea_m_res <- gsea_res %>% as_tibble()
gsea_m_res$rep_genes <- map_chr(gsea_m_res$leadingEdge, function(x) paste(x[1:5], collapse=','))

f_res <- gsea_m_res %>% filter(pval < 0.05) 

x1 <- f_res %>% arrange(NES) %>% head(n=10) %>% pull(pathway)
x2 <- f_res %>% arrange(NES) %>% tail(n=10) %>% pull(pathway)
x_order = c(x1,x2)

ggplot(f_res, aes(x=pathway, y=NES, fill= pval))+
  geom_bar(stat='identity')+
  scale_x_discrete(limits = x_order)+
  coord_flip()+
  theme_bw()



# TDS_low vs TDS_high - GSEA
# DEG
col_dt <- tumor_meta_dt %>%  select(sample_id, tds_0_low_high, tissue_group2) %>% column_to_rownames('sample_id')
dim(col_dt) 

all(colnames(ec_bc_mx) == rownames(col_dt))
s_mx <- ec_bc_mx[,rownames(col_dt)]
all(colnames(s_mx) == rownames(col_dt))
dim(s_mx)

bc_dds <- DESeqDataSetFromMatrix(countData = s_mx,
                                 colData = col_dt,
                                 design = ~ tds_0_low_high)
bc_dds <- DESeq(bc_dds)

res <- results(bc_dds, contrast = c('tds_0_low_high','low','high'))

mcols(res, use.names=TRUE)
m_res <- res %>% as.data.frame() %>% rownames_to_column('gene_name') %>% as_tibble() %>% arrange(desc(log2FoldChange))

tmp_dt <- m_res %>% mutate(nlpval = -log10(padj))
tmp_dt$nlpval[is.finite(tmp_dt$nlpval) == F] <- 300
tmp_dt <- tmp_dt %>% mutate(dist= sqrt(log2FoldChange^2 + nlpval ^2))
tmp_dt <- tmp_dt %>% mutate(dist = ifelse(log2FoldChange < 0, (-1)*dist, dist))

dds_stats <- tmp_dt$dist
names(dds_stats) <- tmp_dt$gene_name

gsea_res <- fgsea(pathways = kegg_list, stats = dds_stats, maxSize=1000)

gsea_m_res <- gsea_res %>% as_tibble()
gsea_m_res$rep_genes <- map_chr(gsea_m_res$leadingEdge, function(x) paste(x[1:5], collapse=','))

f_res <- gsea_m_res %>% filter(pval < 0.05) 

x1 <- f_res %>% arrange(NES) %>% head(n=10) %>% pull(pathway)
x2 <- f_res %>% arrange(NES) %>% tail(n=10) %>% pull(pathway)
x_order = c(x1,x2)

ggplot(f_res, aes(x=pathway, y=NES, fill= pval))+
  geom_bar(stat='identity')+
  scale_x_discrete(limits = x_order)+
  coord_flip()+
  theme_bw()



# normal vs tumor KEGG --> metabolic kegg genes volcano plot
gsea_m_res <- read_excel('/path/to/tumor_vs_normal_KEGG.xlsx')
f_gsea <- gsea_m_res %>% filter(pval < 0.05 & abs(NES) > 1.5) %>% arrange(desc(pval)) 
meta_kegg <- f_gsea$pathway[grepl('METABOLISM',f_gsea$pathway)]

c2_list <- GSEABase::getGmt('/path/to/c2.all.v7.5.1.symbols.gmt') %>% GSEABase::geneIds()
kegg_names <- names(c2_list)[grepl('KEGG',names(c2_list))]
kegg_list <- c2_list[kegg_names]

tumor_metab_kegg_genes <- unique(unlist(kegg_list[meta_kegg], use.name=F)) 

deg_res <- read_excel('/path/to/tumor_vs_normal_DEG.xlsx')
f_dt <- deg_res %>% filter(gene_name %in% tumor_metab_kegg_genes) 
genes <- f_dt %>% filter(padj<0.01) 

label_dt <- genes %>% filter(gene_name %in% label_list)

g1=c('Tumor')
g2=c('Normal')

ggplot(f_dt , aes(x=log2FoldChange, y=-log10(padj)))+
  geom_hline(yintercept = -log10(0.05), linetype = 'dashed', color='red')+
  geom_text_repel(data= genes, aes(label=gene_name))+
  geom_point()+
  theme_se+
  xlab(paste0(paste(g2,collapse=','),'<----->',paste(g1, collapse=',')))+
  theme(legend.position = 'none')



# TDS-low vs TDS-high KEGG --> metabolic kegg genes volcano plot
gsea_m_res <- read_excel('/path/to/TDS_0_DEG_primaryT_KEGG.xlsx')
f_gsea <- gsea_m_res %>% filter(pval < 0.05 & abs(NES) > 1.5) %>% arrange(desc(pval)) 
meta_kegg <- f_gsea$pathway[grepl('METABOLISM',f_gsea$pathway)]

c2_list <- GSEABase::getGmt('/path/to/c2.all.v7.5.1.symbols.gmt') %>% GSEABase::geneIds()
kegg_names <- names(c2_list)[grepl('KEGG',names(c2_list))]
kegg_list <- c2_list[kegg_names]

tds_metab_kegg_genes <- unique(unlist(kegg_list[meta_kegg], use.name=F)) 

deg_res <- read_excel('/path/to/TDS_0_DEG_primaryT.xlsx')
f_dt <- deg_res %>% filter(gene_name %in% tds_metab_kegg_genes)
genes <- f_dt %>% filter(padj<0.01)

label_dt <- genes %>% filter(gene_name %in% label_list)

g1=c('TDS_low')
g2=c('TDS_high')

ggplot(f_dt , aes(x=log2FoldChange, y=-log10(padj)))+
  geom_hline(yintercept = -log10(0.05), linetype = 'dashed', color='red')+
  geom_text_repel(data= genes, aes(label=gene_name))+
  geom_point()+
  theme_se+
  xlab(paste0(paste(g2,collapse=','),'<----->',paste(g1, collapse=',')))+
  theme(legend.position = 'none')


# double volcano plot
# union_list (TDS metabolic kegg, tumor metabolic kegg)
union_list <- union(tumor_metab_kegg_genes, tds_metab_kegg_genes)

# union_list in tumor DEG
tumor_DEG <- read_excel('/path/to/tumor_vs_normal_DEG.xlsx')
tumor_DEG_metab <- tumor_DEG %>% filter(gene_name %in% union_list)

# union_list in tds DEG
tds_DEG <- read_excel('/path/to//TDS_0_DEG_primaryT.xlsx')
tds_DEG_metab <- tds_DEG %>% filter(gene_name %in% union_list) 

#assign path
dt1 <- tumor_DEG_metab
dt2 <- tds_DEG_metab
nrow(dt1) 
nrow(dt2) 

#edit before merge
dt1 <- dt1 %>% mutate(tn_nlgpval = log10(pvalue)*(-1)) %>%
  mutate(tn_nlgpval = ifelse(log2FoldChange < 0, tn_nlgpval*-1, tn_nlgpval))
nrow(dt1) 
dt2 <- dt2 %>% mutate(hl_nlgpval = log10(pvalue)*-1) %>%
  mutate(hl_nlgpval = ifelse(log2FoldChange < 0, hl_nlgpval*-1, hl_nlgpval))
nrow(dt2) 

dt <- full_join(dt1 %>% select(gene_name, tn_nlgpval), dt2 %>% select(gene_name, hl_nlgpval))
nrow(dt) 
summary(dt$tn_nlgpval)
summary(dt$hl_nlgpval)

dt %>% arrange(desc(hl_nlgpval))

oc_gene_list <- c('PHGDH','PSAT1','PSPH','SHMT1','SHMT2','MTHFD1','MTHFD2','MTHFD1L','MTHFD2L','TYMS',
                  'MTHFR','ALDH1L1','ALDH1L2','SFXN1','SFXN2','SFXN3')

label1 <- dt %>% filter(tn_nlgpval > 30 & abs(hl_nlgpval) >10)
label2 <- label1 %>% filter(gene_name %in% oc_gene_list)

ggplot(dt, aes(x=tn_nlgpval, y=hl_nlgpval))+
  geom_point(color = "black", alpha = 0.2)+
  geom_point(data = label1, fill = "black", colour = "black")+
  geom_text_repel(data = label1, aes(label = gene_name), nudge_x=-5)+
  geom_text_repel(data = label2, aes(label = gene_name), col="red")+
  theme_se


############### primary tumor cohort #############
# complexheatmap
tumor_meta_dt <- read_excel('/path/to/primaryT_TDS_group_forSPSS.xlsx')
tumor_tpm <- read_tsv('/path/to/PTC_PD_ATC_tpm.tsv')

#remain only one sample per donor-tissue_group pair
ct_dt <- tumor_meta_dt %>% group_by(donor_id, tissue_group) %>% summarise(donor_sample_n=n())  
tmp_dt <- left_join(tumor_meta_dt, ct_dt)
tmp_dt2 <- tmp_dt %>% filter(donor_sample_n > 1)

multi_donor_ids <- tmp_dt2$donor_id %>% unique()
rep_ids <- tmp_dt2 %>% group_by(donor_id, tissue_group) %>% summarise(rep_sample_id = head(sample_id,n=1)) %>%
  pull(rep_sample_id)

tmp_dt <- bind_rows(tmp_dt %>% filter(!donor_id %in% multi_donor_ids),
                    tmp_dt %>% filter(sample_id %in% rep_ids))

nrow(tmp_dt)
tmp_dt
target_ids <- tmp_dt  %>% arrange(desc(tds))%>% pull(sample_id)

mx <- tumor_tpm %>% as.data.frame() %>% column_to_rownames('gene_name') %>% as.matrix()
l_mx <- log10(mx + 0.01)
tmp_mx <- l_mx[oc_gene_list, target_ids]

annot_cols <- tumor_meta_dt %>% select(sample_id,tissue_group2, tds) %>%
  dplyr::rename(`Tumor type`=tissue_group2, TDS=tds)

annot_cols <- annot_cols %>% column_to_rownames('sample_id')
annot_cols <- annot_cols[colnames(tmp_mx),]

show_col(pal_d3("category10")(10))
cancer_type_pal <- pal_d3("category10")(10)[c(3,2,6)]
names(cancer_type_pal) <- c('PTC_primary','ATC_primary','PD_primary')

top_anno <- ComplexHeatmap::HeatmapAnnotation(df = annot_cols, col = list(`Tumor type`=cancer_type_pal))

if(T){
  scaled_dt <- t(apply(tmp_mx, 1, scale))
  colnames(scaled_dt) <- colnames(tmp_mx)
}

col_fun = colorRamp2(c(-2,0, 2), c('blue', "white", 'red'))
chm <- ComplexHeatmap::Heatmap(scaled_dt, top_annotation = top_anno, show_row_names = T, 
                               show_column_names = F, 
                               clustering_distance_rows='euclidean', 
                               clustering_method_rows = 'complete', 
                               cluster_columns = F,
                               col=col_fun)
ComplexHeatmap::draw(chm, heatmap_legend_side = "left", annotation_legend_side = "left")


###### correlation
# correlation plot - TDS
cor_tbl <- tmp_mx %>% as.data.frame() %>% rownames_to_column('gene_name') %>% as_tibble()
cor_tbl <- cor_tbl %>% as.data.frame() %>% column_to_rownames('gene_name') %>% 
  t() %>% as.data.frame() %>% rownames_to_column('sample_id') %>% as_tibble()

cor_tbl <- left_join(cor_tbl, tumor_meta_dt %>% select(sample_id, tds, Tumor_size))

res_tbl <- tibble()
for (t_gene in oc_gene_list){
  print(t_gene)
  a <- cor_tbl[[t_gene]]
  b <- cor_tbl[["tds"]]
  res <- cor.test(a,b)
  tmp <-tibble(gene=t_gene, cor=res$estimate, pval=res$p.value)
  res_tbl <- bind_rows(res_tbl, tmp)
}

ggplot(res_tbl, aes(x=cor, y=-log10(pval)))+
  geom_point(size = 3)+
  geom_text_repel(aes(label=gene), size = 5)+
  ggtitle('TDS-1C correlation in our cohort')+
  theme_se + theme(legend.position = 'none')


# correlation plot - Tumor size
res_tbl <- tibble()
for (t_gene in oc_gene_list){
  print(t_gene)
  a <- cor_tbl[[t_gene]]
  b <- cor_tbl[["Tumor_size"]]
  res <- cor.test(a,b)
  tmp <-tibble(gene=t_gene, cor=res$estimate, pval=res$p.value)
  res_tbl <- bind_rows(res_tbl, tmp)
}

ggplot(res_tbl, aes(x=cor, y=-log10(pval)))+
  geom_point(size = 3)+
  geom_text_repel(aes(label=gene), size = 5)+
  ggtitle('Tumor size-1C correlation in our cohort')+
  theme_se + theme(legend.position = 'none')


# corr plot --> primary T
tmp_mx <- cor_tbl %>% select(-Tumor_size) %>% dplyr::rename(TDS = tds) %>% as.data.frame() %>% column_to_rownames('sample_id') %>% as.matrix()

res <- cor(tmp_mx)

col2 <- colorRampPalette(c("#053061", "#2166AC", "#4393C3", "#92C5DE", "#D1E5F0", "#FFFFFF", 
                           "#FDDBC7", "#F4A582", "#D6604D", "#B2182B", "#67001F"))


testRes = cor.mtest(tmp_mx, conf.level = 0.95) 

corrplot(res, type = "upper", order = "original",
         tl.col = "black", tl.srt = 45,
         p.mat = testRes$p, 
         sig.level = c(0.001, 0.01, 0.05), pch.cex = 0.9,
         insig = 'label_sig', pch.col = 'black', col=col2(200))


# TDS low vs high - bar graph
mx <- tumor_tpm %>% as.data.frame() %>% column_to_rownames('gene_name') %>% as.matrix()
l_mx <- log10(mx + 0.01)
l_dt <- l_mx %>% as.data.frame() %>% rownames_to_column('gene_name') %>% as_tibble()
l_dt <- l_dt %>% as.data.frame() %>% column_to_rownames('gene_name') %>% t() %>% as.data.frame() %>% rownames_to_column('sample_id') %>% as_tibble()

tbl <- left_join(tumor_meta_dt %>% select(sample_id, tds_0_low_high), l_dt %>% select(sample_id, oc_gene_list))
tbl %>% group_by(tds_0_low_high) %>% dplyr::count()

show_col(pal_d3("category20b")(20))
tds_pal <- pal_d3("category20b")(20)[c(9,18)]
names(tds_pal) <- c('high','low')

ggplot(tbl, aes(x=tds_0_low_high, y= SHMT2, fill=tds_0_low_high))+
  geom_boxplot(outlier.size=-1)+
  geom_jitter(size=3, alpha=0.2, height=0,width=0.1)+
  stat_compare_means(comparisons=list(c("high", "low")), method="t.test", aes(label=..p.signif..))+
  scale_x_discrete(limits = c("high", "low"))+
  theme_bw() + 
  theme(axis.text = element_text(size=20), axis.title = element_text(size=20), 
        panel.border = element_blank(), axis.line = element_line(), plot.title = element_text(size=22))+
  scale_fill_manual(values = tds_pal)+
  theme(legend.position = 'none')
