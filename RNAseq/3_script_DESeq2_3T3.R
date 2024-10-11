library(DESeq2)
library("biomaRt")

setwd("/home/bmartin/Desktop/RNAseq_analysis/")
cts <- as.matrix(read.csv("./count_matrix_3T3.csv", sep=",", row.names="gene_id",check.names=FALSE))
head(cts,2)
coldata <- read.csv("./metadata_3T3.csv", sep=",", row.names=1)
coldata$condition <- factor(coldata$condition)
dds <- DESeqDataSetFromMatrix(countData = cts, colData = coldata, design = ~condition)
keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep,]
dds$condition <- relevel(dds$condition, ref = "untreated")
dds <- DESeq(dds)
res <- results(dds)

### 1h vs untreated
res <- results(dds, name="condition_CHX_1h_vs_untreated")
res <- lfcShrink(dds, coef = 'condition_CHX_1h_vs_untreated', type = 'apeglm', res = res)

ens <- rownames(res)
enst <- gsub("\\.[0-9]*$","",ens)
mart <- useMart("ensembl", "mmusculus_gene_ensembl")
genes <- getBM(attributes=c("ensembl_transcript_id","external_gene_name","ensembl_gene_id"),
                           filters = "ensembl_gene_id",
                           values = enst,
                           mart = mart)

rownames(genes) <- genes$ensembl_transcript_id
names(genes) = c("transcript_id", "gene_name", "gene_id")
head(genes)
genes <- genes[match(enst, genes$gene_id),]
rownames(res) <- genes$gene_name

resOrdered <- res[order(res$pvalue),]
write.csv(as.data.frame(resOrdered), file = "results_3T3/DGE_1vs0_3T3_results.csv")

### 2h vs untreated
res <- results(dds, name="condition_CHX_2h_vs_untreated")
res <- lfcShrink(dds, coef = 'condition_CHX_2h_vs_untreated', type = 'apeglm', res = res)

ens <- rownames(res)
enst <- gsub("\\.[0-9]*$","",ens)

mart <- useMart("ensembl", "mmusculus_gene_ensembl")
genes <- getBM(attributes=c("ensembl_transcript_id","external_gene_name","ensembl_gene_id"),
                           filters = "ensembl_gene_id",
                           values = enst,
                           mart = mart)

rownames(genes) <- genes$ensembl_transcript_id
names(genes) = c("transcript_id", "gene_name", "gene_id")
genes <- genes[match(enst, genes$gene_id),]
rownames(res) <- genes$gene_name

resOrdered <- res[order(res$pvalue),]
write.csv(as.data.frame(resOrdered), file = "results_3T3/DGE_2vs0_3T3_results.csv")

### 5h vs untreated
res <- results(dds, name="condition_CHX_5h_vs_untreated")
res <- lfcShrink(dds, coef = 'condition_CHX_5h_vs_untreated', type = 'apeglm', res = res)

ens <- rownames(res)
enst <- gsub("\\.[0-9]*$","",ens)

mart <- useMart("ensembl", "mmusculus_gene_ensembl")
genes <- getBM(attributes=c("ensembl_transcript_id","external_gene_name","ensembl_gene_id"),
                           filters = "ensembl_gene_id",
                           values = enst,
                           mart = mart)

rownames(genes) <- genes$ensembl_transcript_id
names(genes) = c("transcript_id", "gene_name", "gene_id")
genes <- genes[match(enst, genes$gene_id),]
rownames(res) <- genes$gene_name

resOrdered <- res[order(res$pvalue),]
write.csv(as.data.frame(resOrdered), file = "results_3T3/DGE_5vs0_3T3_results.csv")

### 48h vs untreated
res <- results(dds, name="condition_CHX_48h_vs_untreated")
res <- lfcShrink(dds, coef = 'condition_CHX_48h_vs_untreated', type = 'apeglm', res = res)

ens <- rownames(res)
enst <- gsub("\\.[0-9]*$","",ens)

mart <- useMart("ensembl", "mmusculus_gene_ensembl")
genes <- getBM(attributes=c("ensembl_transcript_id","external_gene_name","ensembl_gene_id"),
                           filters = "ensembl_gene_id",
                           values = enst,
                           mart = mart)

rownames(genes) <- genes$ensembl_transcript_id
names(genes) = c("transcript_id", "gene_name", "gene_id")
genes <- genes[match(enst, genes$gene_id),]
rownames(res) <- genes$gene_name

resOrdered <- res[order(res$pvalue),]
write.csv(as.data.frame(resOrdered), file = "results_3T3/DGE_48vs0_3T3_results.csv")
