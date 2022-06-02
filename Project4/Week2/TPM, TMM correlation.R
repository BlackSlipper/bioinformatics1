library(corrplot)
load("G:/hanshin/DEGanalysis.RData")

abundances <- cbind(count, genelength =a[,1]/1000)
abundances_1 <- cbind(count, genelength =a[,1]/1000)

#rawdata correlation plot
abundances_2 <- na.omit(abundances_1)
rawdata <- abundances_2[1:6]
N <- cor(rawdata)
corrplot(N, method ="shade")

#####TPM ###########
##################################genelength로 나누기
temp1<- abundances[1,]/abundances$genelength[1]
for (i in 1:nrow(abundances)){
  temp <- abundances[i,]/abundances$genelength[i]
  temp1 <- rbind(temp1,temp)
}
data <- temp1[-c(2),1:12] #중복된 행 제거
View(data)

##################################결측치 제거
data <- na.omit(data)

##################################1e6으로 나누기
normfactor <- colSums(data)/1000000

temp1 <- data[,1]/normfactor[1]
for (i in 1:ncol(data)){
  temp <- data[,i]/normfactor[i]
  temp1 <- cbind(temp1, temp)
}
tpm <- as.data.frame(temp1)
tpm <- tpm[,-1]
colnames(tpm) <- colnames(data)
rownames(tpm) <- rownames(data)
View(tpm)
##################################결측치, 0이 있는 행 제거
data2 <- na.omit(tpm)

##################################TPM Correlation plot 그리기
LECdata <- data2[1:6]
M <- cor(LECdata)
corrplot(M, method = "number")
corrplot(M, method="circle")
corrplot(M, method="square")
corrplot(M, method="ellipse")
corrplot(M, method ="shade")
corrplot(M, method="color")
corrplot(M, method="pie")

corrplot(M,
         method="shade",
         addshade="all",
         shade.col=NA,
         tl.col='red',
         tl.srt=30,
         diag=FALSE,
         addCoef.col="black",
         order= "FPC"
         #"FPC": First Principle Component,
         #"hclust": hierarchical clustering
         #"AOE": Angular Order of Eigenvectors
         )



###############
l2abund <- log2(LECdata +1)


hist(LECdata$SRR5887652, nclass=100)
hist(l2abund$SRR5887655, nclass=100)

TPM <- corrplot(cor(LECdata), method = "shade")
logTPM <- corrplot(cor(l2abund), method = "shade")

par(mfcol=c(1,2))
plot()
####TMM####################
Factors <- Y_LEC$samples$lib.size * Y_LEC$samples$norm.factors
tmmExp <- t(t(Y_LEC$counts)/Factors) * mean(Factors)
View(tmmExp)

par(mfcol=c(1,1))
corrplot(cor(tmmExp), method="shade")
save.image(file ="G:/hanshin/DEGanalysis.RData")
