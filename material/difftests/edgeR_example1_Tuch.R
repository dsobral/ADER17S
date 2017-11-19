rawdata <- read.delim("edgeR_example1_Tuch.tab", check.names=FALSE, stringsAsFactors=FALSE)
library(edgeR)
y <- DGEList(counts=rawdata[,2:7], genes=rawdata[,1])

#Optional annotation part... important for the functional enrichment step...
#library(org.Hs.eg.db)
#idfound <- y$genes$RefSeqID %in% mappedRkeys(org.Hs.egREFSEQ)
#y <- y[idfound,]
#egREFSEQ <- toTable(org.Hs.egREFSEQ)
#m <- match(y$genes$RefSeqID, egREFSEQ$accession)
#y$genes$EntrezGene <- egREFSEQ$gene_id[m]
#egSYMBOL <- toTable(org.Hs.egSYMBOL)
#m <- match(y$genes$EntrezGene, egSYMBOL$gene_id)
#y$genes$Symbol <- egSYMBOL$symbol[m]

#Dealing with transcripts and genes...
#o <- order(rowSums(y$counts), decreasing=TRUE)
#y <- y[o,]
#d <- duplicated(y$genes$Symbol)
#y <- y[!d,]
#nrow(y)
#y$samples$lib.size <- colSums(y$counts)

y <- calcNormFactors(y)
y$samples
plotMDS(y)
Patient <- factor(c(8,8,33,33,51,51))
Tissue <- factor(c("N","T","N","T","N","T"))
data.frame(Sample=colnames(y),Patient,Tissue)
design <- model.matrix(~Patient+Tissue)
rownames(design) <- colnames(y)
y <- estimateDisp(y, design, robust=TRUE)
plotBCV(y)
fit <- glmFit(y, design)
lrt <- glmLRT(fit)
topTags(lrt)
o <- order(lrt$table$PValue)
cpm(y)[o[1:10],]
summary(de <- decideTestsDGE(lrt))

detags <- rownames(y)[as.logical(de)]
plotSmear(lrt, de.tags=detags)
abline(h=c(-1, 1), col="blue")

#this part requires the annotation part before...
#go <- goana(lrt)
#topGO(go, ont="BP", sort="Up", n=30)
