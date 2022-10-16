#GENIE3 for inex data 2020/11/20
library(GENIE3)

##input: count matrix  row: genes, column: cells
#scale data
scale.matrix = function(a){
  for(i in 1:nrow(a)){
    if(!identical(as.numeric(a[i,]),numeric(length = ncol(a)))){      #considering special case: zero vector
      a[i,] = scale(a[i,])
    }
  }
  return(a)
}



## intron and exon input
dir = './more_datasets_0301/hBCell'
load(file = paste(dir,'/intron_exon_data_processed.RData',sep = ''))


# all dataset
exon =  as.matrix(exon.sub)
intron = as.matrix(intron.sub)
scaleexon = scale.matrix(exon)
scaleexon = scaleexon[TFs,]
scaleintron = scale.matrix(intron)
#intron names:  add '-intron'
rownames(scaleintron) = paste(rownames(scaleintron),'-intron',sep = '')
#combine exon and intron
scaledata = rbind(scaleexon,scaleintron)
#run genie3
exprMatrix = scaledata
set.seed(123)
weightMatrix <- GENIE3(exprMatrix,regulators = TFs,nCores = 20) 
save(weightMatrix,file = "./genie3_inex_hBCell.RData")


