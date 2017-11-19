rawdata <- read.delim("edgeR_example4_GSE60450_Lactation-GenewiseCounts.tab", header=TRUE)
row.names(rawdata)<-rawdata$EntrezGeneID
metadata <- read.delim("edgeR_example4_GSE60450_Lactation_metadata.tab", header=TRUE)
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
fr <- fry(y, index=ind, design=design, contrast=con)

res <- glmQLFTest(fit, contrast=con)
barcodeplot(res$table$logFC, ind[[1]], main=names(ind)[1])

