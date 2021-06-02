BiocManager::install(version = "3.12")
BiocManager::install("tximport")
BiocManager::install("rhdf5")
BiocManager::install("BiocPkgTools")
BiocManager::install("ggplot")
install.packages("calibrate")
BiocManager::install("DESeq2")
BiocManager::install("topGO")
BiocManager::install("org.Dm.eg.db")
BiocManager::install("multtest")
BiocManager::install("Rgraphviz")
library(tximport)
library(rhdf5)
library(readr)
library(ggplot2)
library(calibrate)
library(DESeq2)
library(topGO)
library(org.Dm.eg.db)
library(genefilter)
library(multtest)
library(Rgraphviz)
setwd("C:/Users/luis_/Documents/Académico/Mastér en Bioinformatica y Bioestadística/TFM/Proyecto")
dir<-getwd()

BGER<-read.table("BGER-mod.fna",header=FALSE)
borrar<-seq(1,length(BGER$V1),2)
BGER<-BGER[borrar,]
write.table(BGER,"GeneNames.txt",sep = "\t",quote=FALSE,col.names=c("geneID"),row.names=FALSE)
BGER<-read.table("GeneNames.txt",header=FALSE)
GeneNames<-cbind(BGER,BGER[,1])
colnames(GeneNames)<-c("geneID","transcriptID")
head(GeneNames)
write.table(GeneNames, file ="TransciptNames.txt",sep="\t",quote = FALSE,row.names= FALSE)
Groups<-data.frame(muestra=c("ED0.1","ED0.2","ED1.1","ED1.2","ED2.1","ED2.2"),condicion=c("día-0","día-0","día-1","día-1","día-2","día-2"))
write.table(Groups, file ="Groups.txt",sep="\t",quote = FALSE,row.names= FALSE)
dir<-getwd()

archivos<-file.path(dir,"kallisto",Groups$muestra,"abundance.h5")
names(files)<-Groups$muestra
tx2gene<-read.table("TransciptNames.txt",header = TRUE)
txi<-tximport(archivos,type="kallisto",tx2gene=tx2gene,countsFromAbundance="no")
coldata<-data.frame(condition=factor(rep(c("día-0","día-1","día-2"),each=2)))
rownames(coldata)<-colnames(txi$counts)
ddsTxi<-DESeqDataSetFromTximport(txi,colData=coldata,design=~condition)

ddsTxi<-ddsTxi[rowSums(counts(ddsTxi))>30]

AED<-DESeq(ddsTxi)
result<-results(AED)
write.table(result,file="AED_BGER.lista",quote = FALSE, dec = ".", sep = "\t" )
genesUP<-result[which(result$padj<0.01 & result$log2FoldChange> 1),]
genesDOWN<-result[which(result$padj<0.01 & result$log2FoldChange< -1),]
write.table(genesUP, file = "GenesUP.lista", quote = FALSE, dec = ",", sep ="\t")
write.table(genesDOWN, file = "GenesDOWN.lista", quote = FALSE, dec = ",", sep ="\t")

plotMA(result, ylim=c(-5,5))
rld<-rlog(AED,blind=TRUE)
vsd<-vst(AED,blind=TRUE)
pcaData<-plotPCA(vsd,intgroup=c("condition"), returnData=TRUE)
percentVar<-round(100 * attr(pcaData, "percentVar"))
ggplot(pcaData, aes(PC1, PC2, color=condition)) +
  geom_point(size=3) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) + 
  coord_fixed()
with(result, plot(log2FoldChange, -log10(pvalue), pch=20, main="Volcano plot", xlim=c(-10,10)))
with(subset(result, padj<0.01), points(log2FoldChange, -log10(pvalue), pch=20, col="blue"))

orthologous<-read.table("1_2_1_orthologs.tab",header = FALSE,sep = "\t",row.names = NULL)
names(orthologous)<-c("Drosophila","B_germanica")
EXP<-read.table("AED_BGER-mod.list",header = TRUE,sep = "\t",dec = ",",row.names = NULL)
EXP<-data.frame(EXP$row.names,EXP$pvalue,row.names = NULL)
names(EXP)<-c("B_germanica", "pvalue")
orthologous_AED<-orthologous$Drosophila[match(EXP$B_germanica,orthologous$B_germanica)]
orthologous_AED<-orthologous_AED[!is.na(orthologous_AED)]
EXPpvalue<-orthologous$B_germanica[match(EXP$B_germanica,orthologous$B_germanica)]
EXPpvalue<-EXPpvalue[!is.na(EXPpvalue)]
EXPpvalue<-EXP$pvalue[match(EXPpvalue,EXP$B_germanica)]
orthologous_AED<-data.frame(orthologous_AED,EXPpvalue)
names(orthologous_AED)<-c("Drosophila","pvalue")
write.table(orthologous_AED,file="Drosophila_AED.list",sep = "\t",col.names=c("genes","pvalue"),row.names =FALSE)

Dlist<-as.numeric(orthologous_AED$pvalue)
names(Dlist)<-orthologous_AED$Drosophila
topDiffGenes <- function(allScore) {
  return(allScore < 0.01)
  }
GOdata<-new("topGOdata", description = "GO terms para B.germanica",ontology = "BP", allGenes = Dlist, geneSel = topDiffGenes, nodeSize = 10, annot = annFUN.org, ID = "ensembl", mapping = "org.Dm.eg")
resultFis<-runTest(GOdata, algorithm = "classic", statistic = "fisher")
resultKS.elim<-runTest(GOdata, algorithm = "elim", statistic = "ks")
allRes<-GenTable(GOdata, classicFisher = resultFis, elimKS = resultKS.elim, orderBy = "elimKS", ranksOf = "classicFisher", topNodes = 25)
write.table(allRes,file="GOresults.list",sep = "\t",col.names=names(allRes),row.names =FALSE)
showSigOfNodes(GOdata, score(resultFis), firstSigNodes = 19, useInfo = 'def')
showSigOfNodes(GOdata, score(resultKS.elim), firstSigNodes = 5, useInfo = 'def')

