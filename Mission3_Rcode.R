#Set Environment
setRepositories(ind=1:7)
library(data.table)
library(DescTools)

#Load Data
setwd("C:/Users/User/Downloads/")
let7d_gene <- read.csv("let7d_pileup.csv")
let7g_gene <- read.csv("let7g_pileup.csv")
let7f_1_gene <- read.csv("let7f_1_pileup.csv")


#Reconstruct Data
data1 <- let7d_gene[,c(1,2,7)]
data2 <- let7g_gene[,c(1,2,7)]
data3 <- let7f_1_gene[,c(1,2,7)]

#1 각 position별로 base수를 셉니다.
baseCount <- c()
for (i in 1:nrow(data1)){
  baseCount <- c(baseCount, nchar(data1$matches[i]))
}
data1 <- cbind(data1, baseCount)

baseCount <- c()
for (i in 1:nrow(data2)){
  baseCount <- c(baseCount, nchar(data2$matches[i]))
}
data2 <- cbind(data2, baseCount)

baseCount <- c()
for (i in 1:nrow(data3)){
  baseCount <- c(baseCount, nchar(data3$matches[i]))
}
data3 <- cbind(data3, baseCount)

#2 각 position별로 shannon entropy를 계산합니다. DescTools 패키지를 이용

Entropylist <- c()
for (i in 1:nrow(data1)){
  baselist <- strsplit(data1$matches[i],"")
  mtx <- as.matrix(table(baselist))
  Entropy <- Entropy(mtx, base=2)
  Entropylist <- c(Entropylist, Entropy)
}
data1 <- cbind(data1, Entropylist)

Entropylist <- c()
for (i in 1:nrow(data2)){
  baselist <- strsplit(data2$matches[i],"")
  mtx <- as.matrix(table(baselist))
  Entropy <- Entropy(mtx, base=2)
  Entropylist <- c(Entropylist, Entropy)
}
data2 <- cbind(data2, Entropylist)

Entropylist <- c()
for (i in 1:nrow(data3)){
  baselist <- strsplit(data3$matches[i],"")
  mtx <- as.matrix(table(baselist))
  Entropy <- Entropy(mtx, base=2)
  Entropylist <- c(Entropylist, Entropy)
}
data3 <- cbind(data3, Entropylist)
View(data1)

#Bedgrpah format으로 출력하기

bedgraph_let7d <- data1[,c(1,2,2,5)]
colnames(bedgraph_let7d) <- c('chrom','chromStart','chromEnd','dataValue')

bedgraph_let7g <- data2[,c(1,2,2,5)]
colnames(bedgraph_let7g) <- c('chrom','chromStart','chromEnd','dataValue')

bedgraph_let7f_1 <- data3[,c(1,2,2,5)]
colnames(bedgraph_let7f_1) <- c('chrom','chromStart','chromEnd','dataValue')

write.table(bedgraph_let7d, 'bedgraph_let7d.csv', quote = F, sep =" ", row.names = F)
write.table(bedgraph_let7g, 'bedgraph_let7g.csv',quote = F, sep =" ", row.names = F)
write.table(bedgraph_let7f_1, 'bedgraph_let7f_1.csv', quote = F, sep =" ", row.names = F)

