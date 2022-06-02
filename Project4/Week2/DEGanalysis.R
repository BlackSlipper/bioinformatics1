################
# 2022-05-30   #
# Hanshin Shin #
# DEG analysis #
################

setRepositories(ind=1:8)
library(dplyr)
library(ggplot2)
library(data.table)
library(edgeR)
library(openxlsx)
library(biomaRt)
library(EDASeq)

#Loading featureCount data
getwd()
setwd("G:/hanshin/5.featureCounts/")
filelist <- list.files(pattern = ".txt$")

data <- fread(filelist[1])
colnames(data)[7] <- strsplit(strsplit(colnames(data)[7],"/")[[1]][6],"_")[[1]][1]
data <- data[,c(1,7)]
a
#save.image(file = "G:/hanshin/DEGanalysis.RData")
for(i in 2:length(filelist)){
  temp <- fread(filelist[i])
  colnames(temp)[7] <- strsplit(strsplit(colnames(temp)[7],"/")[[1]][6],"_")[[1]][1]
  temp <- temp[,7]
  data <- cbind(data,temp)
}
count <- data
rowname <- count[,1]
count = data.frame(data, row.names = rowname$Geneid)
count = count[,-1]

geneSymbol <- data$Geneid

#Loading metadata
file <- c()
for (c in 1:length(filelist)){
  file[c] <- fread(filelist[c])
}
metadata <- fread("G:/hanshin/SraRunTable.txt")

# a = LEC vector control 
# b = LEC V5/vIRF3 overexpression T test
# c = BEC vector control 
# d = BEC V5/vIRF3 overexpression T test


#Loading gene annotation

Annotation <- fread("G:/hanshin/mart_export.txt",sep=",")
Annotation <- as.data.frame(Annotation)

temp <- data.frame(geneSymbol)

temp <- merge(
  x = temp,
  y = Annotation,
  by.x = "geneSymbol",
  by.y = "Gene stable ID",
  sort = F
)
Annotation <- temp

all.equal(as.character(Annotation$geneSymbol),geneSymbol)  
data.frame(as.character(Annotation$geneSymbol),geneSymbol)  

# Removing unnecessary chromosome
indexRemoval <- c(grep("KQ",Annotation$`Chromosome/scaffold name`),grep("LSZ",Annotation$`Chromosome/scaffold name`),grep("LGE",Annotation$`Chromosome/scaffold name`),grep("MT",Annotation$`Chromosome/scaffold name`))

#Removing all zero counted genes

indexRemoval <- c(indexRemoval,(which(rowSums(count) == 0)))
count <- count[-indexRemoval,]
Annotation <- Annotation[-indexRemoval,]
geneSymbol <- intersect(rownames(count), geneSymbol)

dim(count)
dim(Annotation)
length(geneSymbol)

colSums(count)

ensembl_list <- Annotation$geneSymbol
a <- getGeneLengthAndGCContent(ensembl_list, "hsa")

#write.csv(count,"CountData.csv")
#Defining experimental variables


metadata_LEC <- metadata[c(1:6),]
count_LEC <- count[,-c(7:12)]



LEC <-factor(metadata_LEC$`genotype/variation`, levels = c("control","overexpression"))


design_LEC <- model.matrix( ~ LEC)

setwd("G:/hanshin/")


#1-2 TMM normalization(LEC_vector_control vs LEC_V5/IRF3 overexpression)
Y_LEC <- DGEList(counts =count_LEC, genes = geneSymbol)
Y_LEC <- calcNormFactors(Y_LEC, method = "TMM")
Y_LEC <- estimateDisp(Y_LEC, design_LEC)

fit_LEC <- glmFit(Y_LEC, design_LEC)
logCPM_LEC <- cpm(Y_LEC, normalized.lib.sizes = TRUE, log = T) #this is TMM table
###########DEG analysis######################################################################
###################################################################
Result_LEC <- glmLRT(fit_LEC, coef =2)
Result_LEC_table <- topTags(Result_LEC, n=dim(logCPM_LEC)[1], sort.by = "none")$table
View(Result_LEC_table)
#write.xlsx(Result_LEC_table, sheetName = "LEC control vs overexpression", file = "Pairwise comparison.xlsx")
#write.csv(Result_LEC_table,"LECvsLECoverexpression.csv")

print(sum(Result_LEC_table$PValue <0.01))
print(sum(Result_LEC_table$FDR < 0.01))
print(sum(Result_LEC_table$FDR < 0.05))
print(sum(Result_LEC_table$logFC>2))
print(sum(Result_LEC_table$logFC))
View(logCPM_LEC)
Result <- cbind(Result_LEC_table, logCPM_LEC, count)
head(Result)
######################DrawPlot###########################


CountData <- count
GeneIDs <- geneSymbol
Condition <- c("Control","Control","Control","Overexpression", "Overexpression", "Overexpression")
SampleID <- colnames(CountData)
MetaData <- data.frame(SampleID, Condition)
str(MetaData)
targets <- MetaData
targets$Condition <- factor(targets$Condition)
str(targets)
levels(targets$Condition)
Condition <- factor(COndition, levels =c("Control", "Overexpression"))

design <- model.matrix(~Condition, data =targets)
y <- DGEList(counts=CountData, gene=GeneIDs)
y <- calcNormFactors(y)
y <- estimateGLMRobustDisp(y, design)
plotBCV(y)

fit <- glmFit(y, design)
lrt <- glmLRT(fit,coef = 2)
Result_table <- topTags(lrt, n=dim(CountData)[1], sort.by="none")$table
head(Result_table)
nrow(CountData)[1]
print(sum(Result_table$PValue < 0.01))
print(sum(Result_table$FDR < 0.01))
print(sum(Result_table$FDR < 0.05))

TMM <- cpm(y, normalized.lib.sizes=TRUE,log=T)
head(TMM)
TMM_ColName <- paste("TMMValue_",gsub(".bam", "",colnames(TMM)),sep="")
colnames(TMM) <- TMM_ColName

head(CountData)
RawCount_ColName <- paste("RawCount_",gsub(".bam", "",colnames(CountData)),sep="")
colnames(CountData) <- RawCount_ColName

Result <- cbind(Result_table,TMM,CountData)
head(Result)
tail(Result)

write.table(Result, "Result.txt", sep="\t", quote=F, row.names=F)


### DrawPlot###
library(ggplot2)
Index <- as.numeric(order(Result$FDR)[1:10]) #Top 10 FDR genes

plot_data <- data.frame(Gene_expression = TMM[Index[1],],
                        Condition,
                        Gene_ID = rep(GeneIDs[Index[1]], length(Condition)))
View(TMM)
GeneIDs[Index[2]]

for(i in 2:10){
  temp <- data.frame(Gene_expression = TMM[Index[i],],
                     Condition,
                     Gene_ID =rep(GeneIDs[Index[i]], length(Condition)))
  
  plot_data <- rbind(plot_data, temp)
}


ggplot(data=plot_data, aes(x=Condition, y=Gene_expression)) +
  geom_boxplot() +
  geom_boxplot(aes(fill = Condition)) +
  facet_wrap(~ Gene_ID, nrow=3, ncol=4, scale = "free") +
  guides(fill=FALSE) +
  ylab("log2 TMM normalized values") +
  xlab("") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))


ggsave("TopPvalue_BoxPlot.tiff", width=8, height=6)


################## MDS plot #####################################
MDS_data <- plotMDS(Y_LEC,top=20000)
head(MDS_data)
plot_data <- data.frame(Condition = LEC, X=MDS_data$x, Y=MDS_data$y)
head(plot_data)
#####################Volcano plot #########################
plot(x= Result_LEC_table$logFC,
     y = -log10(Result_LEC_table$PValue),
     main = "Volcano Plot",
     xlab = "logFoldChange",
     ylab = "logP-Value",
     type ="p",
     col =adjustcolor("black", alpha=0.3),
     pch=1)

print(sum(Result_table$PValue < 0.01))
print(sum(Result_table$logFC < -3))
print(sum((-log10(Result_LEC_table$PValue)>2)& (Result_table$logFC > 3)))
