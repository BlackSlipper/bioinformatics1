library(ggfortify)
library(pheatmap)
library(data.table)
library(dplyr)
library(plyr)
library(ggplot2)
library(FactoMineR)
library(factoextra)
library(corrplot)

pcaResult <- prcomp(t(tmmExp))
pcaResult1 <- prcomp(tmmExp)
autoplot(pcaResult)
autoplot(pcaResult1)
head(pcaResult$x)

biplot(pcaResult, scale = 0)

#calculate total variance explained by each principal component
var_explained <- pcaResult$sdev^2 /sum(pcaResult$sdev^2)

#create scree plot
qplot(c(1:6), var_explained) +
  geom_line() +
  xlab("Principal Component") +
  ylab("Variance Explained") +
  ggtitle("Scree Plot")+
  ylim(0,1)

#PCA() function
outputPCA <- PCA(tmmExp, graph =T)
print(outputPCA)

get_eigenvalue(outputPCA)
fviz_eig(outputPCA, addlabels=TRUE)
