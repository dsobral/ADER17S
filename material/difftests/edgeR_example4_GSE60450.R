
library(edgeR)
rawdata <- read.delim("edgeR_example4_GSE60450_Lactation-GenewiseCounts.tab", header=TRUE)
row.names(rawdata)<-rawdata$EntrezGeneID
y <- DGEList(counts=rawdata[,3:14], genes=rawdata[,1:2])
y <- calcNormFactors(y)
y$samples
plotMDS(y,method = "logFC")
#This is NOT the same as the BCV associated to the logFC between groups
plotMDS(y,method = "bcv")
#If you want to know more...
?edgeR::plotMDS.DGEList

metadata <- read.delim("edgeR_example4_GSE60450_Lactation_metadata.tab", header=TRUE)

design <- model.matrix(~ CellType, data=metadata)
design
rownames(design) <- colnames(y)
y <- estimateDisp(y, design, robust=TRUE)
#These two lines are the "equivalent" to the exact test in the classical model
#But this is fitting a linear model to the data, iteratively improving it (fit)
fit <- glmFit(y, design)
lrt <- glmLRT(fit)

#Genes obtained here is equivalent to pairwise CellTypeL - CellTypeB != 0 (6 replicates each)
topgenes<-topTags(lrt, n=dim(rawdata)[[1]])
table(topgenes$table$FDR<0.05)



#Let's try again
y <- DGEList(counts=rawdata[,3:14], genes=rawdata[,1:2])
y <- calcNormFactors(y)
design <- model.matrix(~ 0 + CellType, data=metadata)
design
y <- estimateDisp(y, design, robust=TRUE)
fit <- glmFit(y, design)
lrt <- glmLRT(fit)

# This is equivalent to CellTypeL - 0 != 0 (6 replicates) - so almost ALL genes are different!
topgenes<-topTags(lrt, n=dim(rawdata)[[1]])
table(topgenes$table$FDR<0.05)



#Let's try again
y <- DGEList(counts=rawdata[,3:14], genes=rawdata[,1:2])
y <- calcNormFactors(y)
design <- model.matrix(~ 0 + CellType, data=metadata)
con <- makeContrasts(CellTypeL - CellTypeB, levels=design)
y <- estimateDisp(y, design, robust=TRUE, contrasts=con)
fit <- glmFit(y, design)
lrt <- glmLRT(fit, contrast = con)
topgenes<-topTags(lrt, n=dim(rawdata)[[1]])
table(topgenes$table$FDR<0.05)


#Let's try again
y <- DGEList(counts=rawdata[,3:14], genes=rawdata[,1:2])
y <- calcNormFactors(y)
design <- model.matrix(~ CellType + Status, data=metadata)



group <- factor(paste0(metadata$CellType, ".", metadata$Status))
y <- DGEList(counts=rawdata[,3:14], genes=rawdata[,1:2], group=group)
colnames(y) <- metadata$Sample
require(org.Mm.eg.db)
names <- select(org.Mm.eg.db,keys=rownames(y),columns="SYMBOL")
y$genes$name<-names$SYMBOL
keep <- rowSums(cpm(y) > 0.5) >= 2
y <- y[keep, , keep.lib.sizes=FALSE]


points <- c(0,1,2,15,16,17)
colors <- rep(c("blue", "darkgreen", "red"), 2)
plotMDS(y, col=colors[group], pch=points[group])
legend("topleft", legend=levels(group), pch=points, col=colors, ncol=2)

design <- model.matrix(~ 0 + group)
colnames(design) <- levels(group)
y <- estimateDisp(y, design, robust=TRUE)
plotBCV(y)

fit <- glmQLFit(y, design, robust=TRUE)
con <- makeContrasts(B.pregnant - B.lactate, levels=design)
qlf <- glmQLFTest(fit, contrast=con)
is.de <- decideTestsDGE(qlf, p.value=0.05)
plotSmear(qlf, de.tags=rownames(qlf)[is.de!=0])
#Deciding extra condition of logfc
tr <- glmTreat(fit, contrast=con, lfc=log2(1.2))
topTags(tr)

#ANOVA-Like for L samples 
con <- makeContrasts(L.PvsL = L.pregnant - L.lactate, L.VvsL = L.virgin - L.lactate, L.VvsP = L.virgin - L.pregnant, levels=design)
anov <- glmQLFTest(fit, contrast=con)
topTags(anov)

#Functional Enrichment
con <- makeContrasts(B.lactate - B.pregnant, levels=design)
qlf <- glmQLFTest(fit, contrast=con)
go <- goana(qlf, species = "Mm")
topGO(go, n=30)

#Looking for specific Pathway enrichment
library(GO.db)
cyt.go <- c("GO:0032465", "GO:0000281", "GO:0000920")
term <- select(GO.db, keys=cyt.go, columns="TERM")
Rkeys(org.Mm.egGO2ALLEGS) <- cyt.go

ind <- ids2indices(as.list(org.Mm.egGO2ALLEGS), fit$genes$EntrezGeneID)

con <- makeContrasts(B.virgin-B.lactate, levels=design)
#fr <- fry(y, index=ind, design=design, contrast=con)

res <- glmQLFTest(fit, contrast=con)
barcodeplot(res$table$logFC, ind[[1]], main=names(ind)[1])

