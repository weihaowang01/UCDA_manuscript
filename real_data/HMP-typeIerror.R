library(dplyr)
library(mia)
library(ggplot2)
library(ggsci)
library(tidyverse)
library(phyloseq)
library(doSNOW)
library(curatedMetagenomicData)
library(MASS)
# source("mbDecoda.R")


# nsim=100
# ncores=50
# 
# data.gen=function(count, Gamma){
#   group=rbinom(length(count[,1]),1,0.5)
#   test.data=list(count, group, Gamma)
#   names(test.data)=c("count", "grp", "Gamma")
#   return(test.data)
# }
# 

################## DATA #################
HMP_2012=curatedMetagenomicData("HMP_2012.relative_abundance", dryrun = FALSE, counts = TRUE)%>%
  mergeData()
altExps(HMP_2012) <-
  splitByRanks(HMP_2012)
# HMP_2012_class=altExp(HMP_2012,"class")
HMP_2012_class=altExp(HMP_2012,"genus")
otutable=HMP_2012_class@assays@data@listData[["relative_abundance"]]
otutable=t(as(otutable,"matrix"))


datasize <- nrow(otutable)
prevalence = apply(as(otutable, "matrix"), 2, function(x) {
  return(sum(x > 0))
})/(datasize)

keepOTUs = which(prevalence>  0.1)
otutable = otutable[,keepOTUs] # 748  89

# as.numeric(colSums(otutable==0)/748)
source("~/weihaowang/project2/helper_function.R")
Error1=matrix(0,ncol=5,nrow=ncol(otutable))
for(i in 1:100){
  count=otutable
  set.seed(i)
  group=rbinom(dim(otutable)[1],1,0.5)

  Y <- t(count)
  u=as.matrix(group)
  Z=cbind(as.matrix(u))
  colnames(Z)=c("u")
  formula <- 'u'
  formulaa <- '~u'
  
  ######## 比较 ########
  lib_cut=1
  res <- LinDA::linda(Y, Z, paste('~',formula),type = "count", alpha=0.05,lib.cut = lib_cut)
  rej.ld2 <- which(res$output[[1]]$reject)
  rej.ld2
  if(length(rej.ld2)!=0){
    Error1[rej.ld2,1]=1+Error1[rej.ld2,1]
  }
  

  source("~/weihaowang/project2/PODA2confun.R")
  rej2 = PODA2(as.matrix(Z),Y,type="hard-modi",formulaa)
  rej2
  print(rej2)
  if(length(rej2)!=0){
    Error1[rej2,2]=1+Error1[rej2,2]
  }


  y_implement = .imputation(Y,Z,"~u")
  ws_cate=(t(t(log(y_implement))-rowMeans(t(log(y_implement)))))
  # rr-clr
  factor.num=est.confounder.num(formulaa,data.frame(Z),t(ws_cate),method="ed",rmax=30,nRepeat=20)
  cate.results=cate(formulaa,data.frame(Z),t(ws_cate), r = factor.num,fa.method="ml",adj.method="rr")
  rej.cate2=which(p.adjust(cate.results$beta.p.value[,2], "BH") < 0.05)
  rej.cate2
  if(length(rej.cate2)!=0){
    Error1[rej.cate2,3]=1+Error1[rej.cate2,3]
  }

  factor.num=est.confounder.num(formulaa,data.frame(Z),t(log(y_implement)),method="ed",rmax=30,nRepeat=20)
  cate.results=cate(~u,data.frame(Z),t(log(y_implement)), r = factor.num,fa.method="ml",adj.method="rr")
  rej.cate.ct2=which(p.adjust(cate.results$beta.p.value[,2], "BH") < 0.05)
  rej.cate.ct2
  if(length(rej.cate.ct2)!=0){
    Error1[rej.cate.ct2,4]=1+Error1[rej.cate.ct2,4]
  }

  ##################################
  ############ ANCOMBC2 ############
  ##################################
  source("~/weihaowang/project2/ancombc2fun.R")
  rej.ancombc2 <- ancombc2.fun.ss(nrow(Y),Y, Z, formula, alpha=0.05,p_adj_method="BH")
  rej.ancombc2
  if(length(rej.ancombc2)!=0){
    Error1[rej.ancombc2,5]=1+Error1[rej.ancombc2,5]
  }

  
}
Error1/100