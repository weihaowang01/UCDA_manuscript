##### 喜报 我要开始做真实数据了 #####

source("~/weihaowang/project2/tool-function/cate/cate.R")
source("~/weihaowang/project2/tool-function/cate/adjust_functions.R")
source("~/weihaowang/project2/tool-function/cate/factor_functions.R")
source("~/weihaowang/project2/tool-function/cate/simulate_functions.R")
source("~/weihaowang/project2/tool-function/cate/ESA_BCV.R")
source("~/weihaowang/project2/tool-function/cate/plot_esaBcv.R")
source("~/weihaowang/project2/tool-function/vst_nc.R")
source("~/weihaowang/project2/helper_function.R")
source("~/weihaowang/project2/PODA2confun.R")
library(curatedMetagenomicData)
library(dplyr)
library(mia)
library(ggplot2)
library(ggsci)
library(tidyverse)
library(phyloseq)
# library(ANCOMBC)
library(MicrobiomeStat)
# library(glmmTMB)
# library(UpSetR)
library(aplot)
library(patchwork)
library(doSNOW)
# library(ggVennDiagram)

metadata=read.csv(file = "~/weihaowang/project2/realdata/T2D_meta.csv",header=T)
QinJ_2012=curatedMetagenomicData("QinJ_2012.relative_abundance", dryrun = FALSE, counts = TRUE)%>%
  mergeData()

altExps(QinJ_2012) <-
  splitByRanks(QinJ_2012)
QinJ_2012_genus=altExp(QinJ_2012,"genus")
otutable=QinJ_2012_genus@assays@data@listData[["relative_abundance"]]

case.id <- metadata[which(metadata$Diabetic..Y.or.N.=="Y"),1]
control.id <- metadata[which(metadata$Diabetic..Y.or.N.=="N"),1]
data.id <- c(case.id,control.id)
group = c(rep(1,length(case.id)),rep(0,length(control.id)))
group <- group[which(data.id%in%colnames(otutable))]
otutable <- otutable[,match(data.id,colnames(otutable))[which(data.id%in%colnames(otutable))]]
Gamma1 <- metadata[match(colnames(otutable),metadata$Sample.ID),2]
Gamma1 = as.numeric(factor(Gamma1))
Gamma1=ifelse(Gamma1 %in% c(1, 2), 1, 2)
Gamma1 <- data.frame(gender=Gamma1)
Gamma2 <- metadata[match(colnames(otutable),metadata$Sample.ID),7]
Gamma2 = as.numeric(factor(Gamma2))
Gamma2 <- data.frame(Diabetic..Y.or.N.=Gamma2)
Gamma3 = as.data.frame(metadata[match(colnames(otutable),metadata$Sample.ID),c(3:6,8:11,13:18)])
# Gamma1 <- as.factor(Gamma) #肠型
Gamma = cbind(Gamma1,Gamma2,Gamma3)
otutable=t(as(otutable,"matrix"))

# all_missing_idx <- sort(unique(unlist(
#   apply(Gamma[,8:15], 2, function(col) which(col %in% c("", "-")))
# )))

all_missing_idx <- sort(unique(unlist(
  apply(Gamma[,-c(6,9,10)], 2, function(col) which(col %in% c("", "-")))
)))


# missing_list <- apply(Gamma[,8:15], 2, function(col) {
#   col_num <- as.numeric(col)
#   which(is.na(col_num))
# })
# 
# all_missing_idx <- sort(unique(unlist(missing_list)))
# 
# all_missing_idx

group=group[-all_missing_idx]
otutable=otutable[-all_missing_idx,]
Gamma=Gamma[-all_missing_idx,]
# Gamma <- apply(Gamma[c(1,3,4,6,9,10)], 2, as.numeric)
Gamma <- apply(Gamma[,-c(2,5,6,9,10)], 2, as.numeric)

datasize <- nrow(otutable)
prevalence = apply(as(otutable, "matrix"), 2, function(x) {
  return(sum(x > 0))
})/(datasize)

keepOTUs = which(prevalence>  0.1)
otutable = otutable[,keepOTUs]
#Gamma=NULL

T2D_QinJ2=list(otutable,group,Gamma)
# saveRDS(T2D_QinJ,"~/weihaowang/project2/realdata/T2Dnew.rds")
# T2D_QinJ=readRDS("~/weihaowang/project2/realdata/T2D.rds")

count=T2D_QinJ[[1]]
group=T2D_QinJ[[2]]
Gamma=T2D_QinJ[[3]]
Y <- t(count)
# Z = cbind(as.matrix(group),as.matrix(Gamma[,c(1,3,4,6,9,10)]))
# Z = cbind(as.matrix(group),as.matrix(Gamma))
# colnames(Z)=c("u",paste0("z",1:ncol(Z[,-1])))
# formula=paste(colnames(Z), collapse = "+")
# formulaa=paste("~",formula)

#
z1=as.matrix(Gamma[,16])
#
#
u=as.numeric(u)
z1=as.numeric(z1)
Z=cbind(as.matrix(u),as.matrix(z1))
colnames(Z)=c("u","z1")
formula <- 'u+z1'
formulaa <- '~u+z1'

######## 比较 ########
lib_cut=1000
res <- LinDA::linda(Y, as.matrix(Z), paste('~',formula),type = "count", alpha=0.05,lib.cut = lib_cut)
# res <-lindaa(Y, Z, paste('~',formula),type = "count", alpha=0.05,lib.cut = lib_cut)
rej.ld1 <- which(res$output[[1]]$reject)
rej.ld1


source("~/weihaowang/project2/PODA2confun.R")
rej1 = PODA2(as.matrix(Z),Y,type="hard-modi")
rej1


y_implement = .imputation(Y,Z,"~u")
ws_cate=(t(t(log(y_implement))-rowMeans(t(log(y_implement)))))
# rr-clr
factor.num=est.confounder.num(formulaa,data.frame(Z),t(ws_cate),method="ed",rmax=30,nRepeat=20)
cate.results=cate(formulaa,data.frame(Z),t(ws_cate), r = factor.num,fa.method="ml",adj.method="rr")
rej.cate1=which(p.adjust(cate.results$beta.p.value[,2], "BH") < 0.05)
rej.cate1

factor.num=est.confounder.num(formulaa,data.frame(Z),t(log(y_implement)),method="ed",rmax=30,nRepeat=20)
cate.results=cate(~u,data.frame(Z),t(log(y_implement)), r = factor.num,fa.method="ml",adj.method="rr")
rej.cate.ct1=which(p.adjust(cate.results$beta.p.value[,2], "BH") < 0.05)
rej.cate.ct1


count=T2D_QinJ2[[1]]
group=T2D_QinJ2[[2]]
Gamma=T2D_QinJ2[[3]]

Y <- t(count)
# Z <- res$Z
u=as.matrix(group)
z1=as.matrix(Gamma)


u=as.numeric(u)
z1=as.numeric(z1)
Z=cbind(as.matrix(u))
colnames(Z)=c("u")
formula <- 'u'
formulaa <- '~u'

######## 比较 ########
lib_cut=1000
res <- LinDA::linda(Y, Z, paste('~',formula),type = "count", alpha=0.05,lib.cut = lib_cut)
rej.ld2 <- which(res$output[[1]]$reject)
rej.ld2


source("~/weihaowang/project2/PODA2.R")
rej2 = PODA2(as.matrix(Z),Y,type="hard-modi")
rej2


y_implement = .imputation(Y,Z,"~u")
ws_cate=(t(t(log(y_implement))-rowMeans(t(log(y_implement)))))
# rr-clr
factor.num=est.confounder.num(formulaa,data.frame(Z),t(ws_cate),method="ed",rmax=30,nRepeat=20)
cate.results=cate(formulaa,data.frame(Z),t(ws_cate), r = factor.num,fa.method="ml",adj.method="rr")
rej.cate2=which(p.adjust(cate.results$beta.p.value[,2], "BH") < 0.05)
rej.cate2

factor.num=est.confounder.num(formulaa,data.frame(Z),t(log(y_implement)),method="ed",rmax=30,nRepeat=20)
cate.results=cate(~u,data.frame(Z),t(log(y_implement)), r = factor.num,fa.method="ml",adj.method="rr")
rej.cate.ct2=which(p.adjust(cate.results$beta.p.value[,2], "BH") < 0.05)
rej.cate.ct2
#


length(intersect(rej1, rej2))/length(union(rej1, rej2))
length(intersect(rej.ld1, rej.ld2))/length(union(rej.ld1, rej.ld2))
length(intersect(rej.cate1, rej.cate2))/length(union(rej.cate1, rej.cate2))
length(intersect(rej.cate.ct1, rej.cate.ct2))/length(union(rej.cate.ct1, rej.cate.ct2))
