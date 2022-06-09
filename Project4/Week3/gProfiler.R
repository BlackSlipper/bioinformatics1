################
# 2022-06-08   #
# Hanshin Shin #
# g:Profiler   #
################

setRepositories(ind=1:7)
library(gprofiler2)
library(enrichR)
library(feather)
#g:GOSt
#tempgeneIndex <- as.numeric(order(Result$FDR)[1:100]) #Top 100 FDR genes
#tempgeneIndex <- as.numeric(order(Result$PValue)[1:100]) #Top 100 p-value genes
tempgeneIndex <- which(Result_LEC_table$logFC >4 & Result_LEC_table$FDR<0.05) 




#print(sum(Result_table$PValue < 0.01))
print(sum(Result_LEC_table$logFC >4 & Result_LEC_table$FDR<0.05))
#print(sum(Result_LEC_table$logFC <-2 & Result_LEC_table$FDR<0.0001))


genelist <- c()
for(i in 1:length(tempgeneIndex)){
  
  genelist <- c(genelist, GeneIDs[tempgeneIndex[i]])
  print(i)
}

write.table(genelist, "genelist.csv", sep=" ")
gostres <- gost(query = genelist,
                organism = "hsapiens")

View(gostres$result)

p <- gostplot(gostres, capped = FALSE, interactive = FALSE)

p

