source("~/weihaowang/project2/ancombc2fun.R")
source("~/ANCOMBC-bugfix/R/ancom.R")
source("~/ANCOMBC-bugfix/R/ancombc_bias_correct.R")
source("~/ANCOMBC-bugfix/R/ancombc_mult.R")
source("~/ANCOMBC-bugfix/R/ancombc_prep.R")
source("~/ANCOMBC-bugfix/R/ancombc.R")
source("~/ANCOMBC-bugfix/R/ancombc2.R")
source("~/ANCOMBC-bugfix/R/data.R")
source("~/ANCOMBC-bugfix/R/data_sanity_check.R")
source("~/ANCOMBC-bugfix/R/secom_dist.R")
source("~/ANCOMBC-bugfix/R/secom_linear.R")
source("~/ANCOMBC-bugfix/R/secom_prep.R")
source("~/ANCOMBC-bugfix/R/sim_data.R")
source("~/ANCOMBC-bugfix/R/sparse_dist.R")
source("~/ANCOMBC-bugfix/R/sparse_linear.R")
source("~/ANCOMBC-bugfix/R/utils.R")

source("~/weihaowang/project2/tool-function/cate/cate.R")
source("~/weihaowang/project2/tool-function/cate/adjust_functions.R")
source("~/weihaowang/project2/tool-function/cate/factor_functions.R")
source("~/weihaowang/project2/tool-function/cate/simulate_functions.R")
source("~/weihaowang/project2/tool-function/cate/ESA_BCV.R")
source("~/weihaowang/project2/tool-function/cate/plot_esaBcv.R")
source("~/weihaowang/project2/tool-function/vst_nc.R")
source("~/weihaowang/project2/helper_function.R")
source("~/weihaowang/project2/PODA2confun.R")


library(BiocManager)
library(remotes)
library(devtools)
library(mia)
library(microbiome)
library(magrittr)
library(TreeSummarizedExperiment)
library(dplyr)
library(tidytree)
library(doRNG)
library(rlang)
library(MASS)
library(ggrepel)
library(MASS)
library(S4Vectors)
library(mvtnorm)
library(foreach)
library(doParallel)
library(lmerTest)
library(parallel)
library(ggplot2)
# install.packages("modeest",dependencies = TRUE)
# library(modeest)

# if (!requireNamespace("modeest", quietly = TRUE)) {
#   install.packages("modeest")
# }
# BiocManager::install("modeest")
# devtools::install_github("paulponcet/modeest")
library(modeest)


# remotes::install_version(
#   "rbiom",
#   version = "2.2.1",
#   upgrade = "never"
# )
# BiocManager::install("microbiome")
# BiocManager::install("mia")

# # BiocManager::install(version = '3.16')
devtools::install_github("zhouhj1994/LinDA")
library(LinDA)
if (!requireNamespace("glmnet", quietly = TRUE)) {
  install.packages("glmnet")
}
library(glmnet)

pppp=NULL
ffff=NULL
ttime=100
fdr.all = matrix(0,nrow=ttime,ncol=7)
power.all = matrix(0,nrow=ttime,ncol=7)
aaaa=NULL
mumu=2
gammain=0.1
nnum=100
for(mumu in c(2)){
  for(gammain in c(0.05,0.1,0.15,0.2)){
    for(nnum in c(50,100,200)){
      for(qq in c(1:ttime)){
        q=qq
        source("~/weihaowang/project2/try/tmp_houtai1.R")
        # source("~/weihaowang/project2/try/生成初始数据ancombc2.R")
        H.keep=H
        
        
        ############################# 数据初步处理 ############################
        confun = NULL
        # confun=as.matrix(z1)
        # confun=cbind(as.matrix(z1),as.matrix(z2))
        pre_cut = 0.1
        sample_cut = 1000
        
        
        n=ncol(Y)
        m=nrow(Y)
        
        ############# sample filtration ###############
        if(is.null(confun)){
          sample_fil = which(colSums(Y)>=sample_cut)
          Y=Y[,sample_fil]
          X_true=X_true[,sample_fil]
          Z.ld=Z.ld[sample_fil,]
          u=u[sample_fil]
          Z=as.matrix(u)
          colnames(Z)=c("u")
          formula="u"
          formulaa="~u"
          ZZ=data.frame(Z)
          ZZ[,1]=factor(ZZ[,1])
          colnames(ZZ)=c("u")
        }else{
          confun=as.matrix(confun)
          sample_fil = which(colSums(Y)>=sample_cut)
          Y=Y[,sample_fil]
          confun=as.matrix(confun[sample_fil,])
          colnames(confun)=paste0("z",1:ncol(confun))
          u=u[sample_fil]
          Z=cbind(u,confun)
          colnames(Z)=c("u",paste0("z",1:ncol(confun)))
          formula=paste(colnames(Z), collapse = "+")
          formulaa=paste("~",formula)
          ZZ=data.frame(Z)
          ZZ[,1]=factor(ZZ[,1])
          colnames(ZZ)=c("u",paste0("z",1:ncol(confun)))
        }
        
        ############# taxa filtration ################
        datasize <- ncol(Y)
        prevalence = apply(as(Y, "matrix"), 1, function(x) {
          return(sum(x > 0))
        })/(datasize)
        
        keepOTUs = which(prevalence> pre_cut)
        Y.ori=Y
        Y = Y[keepOTUs,]
        Y.ld=Y
        H.keep= H[keepOTUs]
        # print(c("H",sum(H)))
        # print(c("Hkeep",sum(H.keep)))
        
        n=ncol(Y)
        m=nrow(Y)
        print(dim(Y))
        rej_hm = PODA2(Z,Y,type="hard-modi",formulaa)
        power.all[qq,2] = sum(which(H.keep!=0)%in%rej_hm)/length(which(H.keep!=0))
        if(length(rej_hm)==0){
          fdr.all[qq,2] = 0
        }else{
          fdr.all[qq,2] = (length(rej_hm)-sum(which(H.keep!=0)%in%rej_hm))/length(rej_hm)
        }
        
        print(c(power.all[qq,2],fdr.all[qq,2]))
        
        
        # rej_h = PODA2(Z,Y,type="soft-modi")
        # power.all[qq,3] = sum(which(H.keep!=0)%in%rej_h)/length(which(H.keep!=0))
        # if(length(rej_h)==0){
        #   fdr.all[qq,3] = 0
        #  }else{
        #   fdr.all[qq,3] = (length(rej_h)-sum(which(H.keep!=0)%in%rej_h))/length(rej_h)
        #  }
        
        
        
        ###################################################################
        ############################### true ##############################
        ###################################################################
        Y.true=X_true[keepOTUs,] #假设观测得到绝对丰度X_true
        
        # 假设数据已经准备好
        # Y: 矩阵，每行是一个响应变量的观测值
        # u: 向量，关注的协变量
        # z1, z2: 其他协变量
        
        
        
        # 获取所有响应变量的回归结果
        results <- test_u_coefficients(log(Y.true), Z.ld)
        results$p_adjusted_BH <- p.adjust(results$p_value, method = "BH")
        significant_vars <- results[results$p_adjusted_BH < 0.05, ]
        rej.true=significant_vars$response_var
        power.true=sum(H.keep[rej.true])/sum(H.keep)
        fdr.true=(length(H.keep[rej.true])-sum(H.keep[rej.true]))/length(H.keep[rej.true])
        # print(c(power.true,fdr.true))
        power.all[qq,1]=power.true
        if(length(H.keep[rej.true])==0){
          fdr.all[qq,1]=0
        }else{
          fdr.all[qq,1]=fdr.true
        }
        
        
        
        
        ##################################################################
        ########################## LinDA series ###########################
        ##################################################################
        # false
        # res=LinDA::linda(Y.ori, data.frame(u=factor(u)), paste('~',formula),type = "count", alpha=0.05,lib.cut = 0,prev.cut = 0,p.adj.method = "BH")
        # res=LinDA::linda(Y, data.frame(u=u), paste('~',formula),type = "count", alpha=0.05,lib.cut = 0,prev.cut = 0,p.adj.method = "BH")
        res=LinDA::linda(Y, data.frame(Z), paste('~',formula),type = "count", alpha=0.05,lib.cut = 0,prev.cut = 0,p.adj.method = "BH")
        rej.ld <- which(res$output[[1]]$reject)
        rejj.ld=as.numeric(gsub("taxon","",rownames(res$output[[1]])[which(res$output[[1]]$reject=="TRUE")]))
        power.ld=sum(H.keep[rej.ld])/sum(H.keep)
        fdr.ld=(length(H.keep[rej.ld])-sum(H.keep[rej.ld]))/length(rej.ld)
        
        
        power.all[qq,4]=power.ld
        if(length(H.keep[rej.ld])==0){
          fdr.all[qq,4]=0
        }else{
          fdr.all[qq,4]=fdr.ld
        }
        print(c(fdr.ld,power.ld))
        
        
        
        
        # ############################################################################
        # ############################## CATE-relative ###############################
        # ############################################################################
        #
        #
        y_implement = .imputation(Y,Z,formulaa)
        ws_cate=(t(t(log(y_implement))-rowMeans(t(log(y_implement)))))
        # rr-clr
        tt <- try( est.confounder.num(formulaa,data.frame(Z),t(ws_cate),method="ed",rmax=30,nRepeat=20)
                   , silent = TRUE)
        if(class(tt)=="try-error"){
          fdr.all[qq,5]=0
          power.all[qq,5]=0
        }else{
          factor.num=est.confounder.num(formulaa,data.frame(Z),t(ws_cate),method="ed",rmax=30,nRepeat=20)
          tt1 = try(   cate(formulaa,data.frame(Z),t(ws_cate), r = factor.num,fa.method="ml",adj.method="rr")
                       , silent = TRUE)
          if(class(tt1)=="try-error"){
            fdr.all[qq,5]=0
            power.all[qq,5]=0
          }else{
            factor.num=est.confounder.num(formulaa,data.frame(Z),t(ws_cate),method="ed",rmax=30,nRepeat=20)
            cate.results=cate(formulaa,data.frame(Z),t(ws_cate), r = factor.num,fa.method="ml",adj.method="rr")
            rej.cate=which(p.adjust(cate.results$beta.p.value[,2], "BH") < 0.05)
            power.cate=sum(H.keep[rej.cate])/sum(H.keep)
            fdr.cate=(length(H.keep[rej.cate])-sum(H.keep[rej.cate]))/length(H.keep[rej.cate])
            power.all[qq,5]=power.cate
            if(length(H.keep[rej.cate])==0){
              fdr.all[qq,5]=0
            }else{
              fdr.all[qq,5]=fdr.cate
            }
          }}
        #
        
        
        #########################################################################
        ############################## CATE-count ###############################
        #########################################################################
        # rr-clr
        tt <- try( est.confounder.num(formulaa,data.frame(Z),t(log(y_implement)),method="ed",rmax=30,nRepeat=20)
                   , silent = TRUE)
        if(class(tt)=="try-error"){
          fdr.all[qq,6]=0
          power.all[qq,6]=0
        }else{
          factor.num=est.confounder.num(formulaa,data.frame(Z),t(log(y_implement)),method="ed",rmax=30,nRepeat=20)
          tt1 = try(   cate(formulaa,data.frame(Z),t(log(y_implement)), r = factor.num,fa.method="ml",adj.method="rr")
                       , silent = TRUE)
          if(class(tt1)=="try-error"){
            fdr.all[qq,6]=0
            power.all[qq,6]=0
          }else{
            factor.num=est.confounder.num(formulaa,data.frame(Z),t(log(y_implement)),method="ed",rmax=30,nRepeat=20)
            cate.results=cate(~u,data.frame(Z),t(log(y_implement)), r = factor.num,fa.method="ml",adj.method="rr")
            rej.cate=which(p.adjust(cate.results$beta.p.value[,2], "BH") < 0.05)
            power.cate=sum(H.keep[rej.cate])/sum(H.keep)
            fdr.cate=(length(H.keep[rej.cate])-sum(H.keep[rej.cate]))/length(H.keep[rej.cate])
            power.all[qq,6]=power.cate
            if(length(H.keep[rej.cate])==0){
              fdr.all[qq,6]=0
            }else{
              fdr.all[qq,6]=fdr.cate
            }
          }}
        
        
        #########################################################################
        ############################## ANCOM-BC2 ################################
        #########################################################################
        tt <- try( ancombc2.fun.ss(nrow(Y),Y, ZZ, formula, alpha=0.05,p_adj_method="BH")
                   , silent = TRUE)
        if(class(tt)=="try-error"){
          fdr.all[qq,7]=0
          power.all[qq,7]=0
        }else{
          rej.ancombc <- tt
          fdr.ancombc=(length(H.keep[rej.ancombc])-sum(H.keep[rej.ancombc]))/length(H.keep[rej.ancombc])
          power.ancombc=sum(H.keep[rej.ancombc])/sum(H.keep)
          power.all[qq,7]=power.ancombc
          if(length(H.keep[rej.ancombc])==0){
            fdr.all[qq,7]=0
          }else{
            fdr.all[qq,7]=fdr.ancombc
          }
          print(c(fdr.ancombc,power.ancombc))
        }
        
        
        
        
        
        
        
        print(qq)
      }
      aaa=c(colMeans(fdr.all,na.rm=TRUE),colMeans(power.all,na.rm=TRUE))
      print(aaa)
      
      colnames(power.all)=c("true","hard-modi","hard","linda","cate-rl","cate-count","ancombc")
      colnames(fdr.all)=c("true","hard-modi","hard","linda","cate-rl","cate-count","ancombc")
      # colMeans(power.all)
      # colMeans(fdr.all)
      # pppp=colMeans(power.all)
      # ffff=colMeans(fdr.all)
      pppp=rbind(pppp,colMeans(power.all,na.rm=TRUE))
      ffff=rbind(ffff,colMeans(fdr.all,na.rm=TRUE))
      saveRDS(pppp,"~/weihaowang/project2/simufinal/power_sqrt2_depstrong.rds")
      saveRDS(ffff,"~/weihaowang/project2/simufinal/fdr_sqrt2_depstrong.rds")
      print(c(mu,nnum,gammain))
      print(colMeans(power.all))
      print(colMeans(fdr.all))
    }
  }
}
saveRDS(pppp,"~/weihaowang/project2/simufinal/power_sqrt2_depstrong.rds")
saveRDS(ffff,"~/weihaowang/project2/simufinal/fdr_sqrt2_depstrong.rds")