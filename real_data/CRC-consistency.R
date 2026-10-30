cons2 = NULL
cons.ld2 = NULL
cons.cate2 = NULL
cons.ancombc2 = NULL

lib_cut=0


Y=t(CRC_zackular[[1]])
u=CRC_zackular[[2]]

Z=cbind(as.matrix(u))
ZZ=data.frame(factor(u))
confoun=NULL
colnames(Z)=c("u")
colnames(ZZ)=c("u")
formula <- 'u'
formulaa <- '~u'
res <- LinDA::linda(Y, ZZ, paste('~',formula),type = "count", alpha=0.05,lib.cut = lib_cut)
rej.ld <- which(res$output[[1]]$reject)
rej.ld


source("~/weihaowang/project2/PODA2confun.R")
rej = PODA2(as.matrix(Z),Y,type="hard-modi",formulaa,seed=1234)
rej


y_implement = .imputation(Y,Z,"~u")
ws_cate=(t(t(log(y_implement))-rowMeans(t(log(y_implement)))))
# rr-clr
factor.num=est.confounder.num(formulaa,data.frame(Z),t(ws_cate),method="ed",rmax=30,nRepeat=20)
cate.results=cate(formulaa,data.frame(Z),t(ws_cate), r = factor.num,fa.method="ml",adj.method="rr")
rej.cate=which(p.adjust(cate.results$beta.p.value[,2], "BH") < 0.05)
rej.cate

factor.num=est.confounder.num(formulaa,data.frame(Z),t(log(y_implement)),method="ed",rmax=30,nRepeat=20)
cate.results=cate(~u,data.frame(Z),t(log(y_implement)), r = factor.num,fa.method="ml",adj.method="rr")
rej.cate.ct=which(p.adjust(cate.results$beta.p.value[,2], "BH") < 0.05)
rej.cate.ct
#

rej.ancombc <- ancombc2.fun.ss(nrow(Y),Y, ZZ, formula, alpha=0.05,p_adj_method="BH")
rej.ancombc


################################################################
############## DA DETECTION AND CONSISTENCY 1 ##################
################################################################
lib_cut=0


Y=t(CRC_zackular[[1]])
u=CRC_zackular[[2]]
z1=CRC_zackular[[3]]
z2=CRC_zackular[[4]] # factor

Z=cbind(as.matrix(u),as.matrix(z1))
ZZ=data.frame(factor(u),as.matrix(z1))
colnames(Z)=c("u","z1")
colnames(ZZ)=c("u","z1")
formula <- 'u+z1'
formulaa <- '~u+z1'
res <- LinDA::linda(Y, ZZ, paste('~',formula),type = "count", alpha=0.05,lib.cut = lib_cut)
rej.ld2 <- which(res$output[[1]]$reject)
rej.ld2


source("~/weihaowang/project2/PODA2confun.R")
rej2 = PODA2(as.matrix(Z),Y,type="hard-modi",formulaa,seed=1234)
rej2


y_implement = .imputation(Y,Z,"~u+z1")
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

rej.ancombc2 <- ancombc2.fun.ss(nrow(Y),Y, ZZ, formula, alpha=0.05,p_adj_method="BH")
rej.ancombc2


cons2 = c(cons2, length(intersect(rej, rej2))/length(union(rej, rej2)))
cons.ld2 = c(cons.ld2 , length(intersect(rej.ld, rej.ld2))/length(union(rej.ld, rej.ld2)))
cons.cate2 = c(cons.cate2, length(intersect(rej.cate, rej.cate2))/length(union(rej.cate, rej.cate2)))
# cons.cate2 = c(cons.cate2, length(intersect(rej.cate.ct, rej.cate.ct2))/length(union(rej.cate.ct, rej.cate.ct2)))
cons.ancombc2 =  c(cons.ancombc2, length(intersect(rej.ancombc, rej.ancombc2))/length(union(rej.ancombc, rej.ancombc2)))






################################################################
############## DA DETECTION AND CONSISTENCY 2 ##################
################################################################
lib_cut=0


Y=t(CRC_zackular[[1]])
u=CRC_zackular[[2]]
z1=CRC_zackular[[3]]
z2=CRC_zackular[[4]] # factor

Z=cbind(as.matrix(u),as.matrix(z1))
ZZ=data.frame(factor(u),as.matrix(z1))
colnames(Z)=c("u","z2")
colnames(ZZ)=c("u","z2")
formula <- 'u+z2'
formulaa <- '~u+z2'
res <- LinDA::linda(Y, ZZ, paste('~',formula),type = "count", alpha=0.05,lib.cut = lib_cut)
rej.ld3 <- which(res$output[[1]]$reject)
rej.ld3


source("~/weihaowang/project2/PODA2confun.R")
rej3 = PODA2(as.matrix(Z),Y,type="hard-modi",formulaa,seed=1234)
rej3


y_implement = .imputation(Y,Z,"~u+z2")
ws_cate=(t(t(log(y_implement))-rowMeans(t(log(y_implement)))))
# rr-clr
factor.num=est.confounder.num(formulaa,data.frame(Z),t(ws_cate),method="ed",rmax=30,nRepeat=20)
cate.results=cate(formulaa,data.frame(Z),t(ws_cate), r = factor.num,fa.method="ml",adj.method="rr")
rej.cate3=which(p.adjust(cate.results$beta.p.value[,2], "BH") < 0.05)
rej.cate3

factor.num=est.confounder.num(formulaa,data.frame(Z),t(log(y_implement)),method="ed",rmax=30,nRepeat=20)
cate.results=cate(~u,data.frame(Z),t(log(y_implement)), r = factor.num,fa.method="ml",adj.method="rr")
rej.cate.ct3=which(p.adjust(cate.results$beta.p.value[,2], "BH") < 0.05)
rej.cate.ct3
#

rej.ancombc3 <- ancombc2.fun.ss(nrow(Y),Y, ZZ, formula, alpha=0.05,p_adj_method="BH")
rej.ancombc3


# length(intersect(rej, rej3))/length(union(rej, rej3))
# length(intersect(rej.ld, rej.ld3))/length(union(rej.ld, rej.ld3))
# length(intersect(rej.cate, rej.cate3))/length(union(rej.cate, rej.cate3))
# length(intersect(rej.cate.ct, rej.cate.ct3))/length(union(rej.cate.ct, rej.cate.ct3))


cons2 = c(cons2, length(intersect(rej, rej3))/length(union(rej, rej3)))
cons.ld2 = c(cons.ld2 , length(intersect(rej.ld, rej.ld3))/length(union(rej.ld, rej.ld3)))
cons.cate2 = c(cons.cate2, length(intersect(rej.cate, rej.cate3))/length(union(rej.cate, rej.cate3)))
# cons.cate2 = c(cons.cate2, length(intersect(rej.cate.ct, rej.cate.ct3))/length(union(rej.cate.ct, rej.cate.ct3)))
cons.ancombc2 =  c(cons.ancombc2, length(intersect(rej.ancombc, rej.ancombc3))/length(union(rej.ancombc, rej.ancombc3)))

################################################################
############## DA DETECTION AND CONSISTENCY 2 ##################
################################################################
lib_cut=0


Y=t(CRC_zackular[[1]])
u=CRC_zackular[[2]]
z1=CRC_zackular[[3]]
z2=CRC_zackular[[4]] # factor

Z=cbind(as.matrix(u),as.matrix(z1),as.matrix(z2))
ZZ=data.frame(factor(u),as.matrix(z1),factor(z2))
colnames(Z)=c("u","z1","z2")
colnames(ZZ)=c("u","z1","z2")
formula <- 'u+z1+z2'
formulaa <- '~u+z1+z2'
res <- LinDA::linda(Y, ZZ, paste('~',formula),type = "count", alpha=0.05,lib.cut = lib_cut)
rej.ld4 <- which(res$output[[1]]$reject)
rej.ld4


source("~/weihaowang/project2/PODA2confun.R")
rej4 = PODA2(as.matrix(Z),Y,type="hard-modi",formulaa,seed=1234)
rej4


y_implement = .imputation(Y,Z,"~u+z1+z2")
ws_cate=(t(t(log(y_implement))-rowMeans(t(log(y_implement)))))
# rr-clr
factor.num=est.confounder.num(formulaa,data.frame(Z),t(ws_cate),method="ed",rmax=30,nRepeat=20)
cate.results=cate(formulaa,data.frame(Z),t(ws_cate), r = factor.num,fa.method="ml",adj.method="rr")
rej.cate4=which(p.adjust(cate.results$beta.p.value[,2], "BH") < 0.05)
rej.cate4

factor.num=est.confounder.num(formulaa,data.frame(Z),t(log(y_implement)),method="ed",rmax=30,nRepeat=20)
cate.results=cate(~u,data.frame(Z),t(log(y_implement)), r = factor.num,fa.method="ml",adj.method="rr")
rej.cate.ct4=which(p.adjust(cate.results$beta.p.value[,2], "BH") < 0.05)
rej.cate.ct4
#

rej.ancombc4 <- ancombc2.fun.ss(nrow(Y),Y, ZZ, formula, alpha=0.05,p_adj_method="BH")
rej.ancombc4


# length(intersect(rej, rej4))/length(union(rej, rej4))
# length(intersect(rej.ld, rej.ld4))/length(union(rej.ld, rej.ld4))
# length(intersect(rej.cate, rej.cate4))/length(union(rej.cate, rej.cate4))
# length(intersect(rej.cate.ct, rej.cate.ct4))/length(union(rej.cate.ct, rej.cate.ct4))

cons2 = c(cons2, length(intersect(rej, rej4))/length(union(rej, rej4)))
cons.ld2 = c(cons.ld2 , length(intersect(rej.ld, rej.ld4))/length(union(rej.ld, rej.ld4)))
cons.cate2 = c(cons.cate2, length(intersect(rej.cate, rej.cate4))/length(union(rej.cate, rej.cate4)))
# cons.cate2 = c(cons.cate2, length(intersect(rej.cate.ct, rej.cate.ct4))/length(union(rej.cate.ct, rej.cate.ct4)))
cons.ancombc2 =  c(cons.ancombc2, length(intersect(rej.ancombc, rej.ancombc4))/length(union(rej.ancombc, rej.ancombc4)))


value=c(mean(cons2),mean(cons.ld2),mean(cons.cate2),mean(cons.ancombc2))

method=c("PoDA", "LinDA", "CATE","ANCOM-BC2")
method=factor(method,levels=c("PoDA", "LinDA", "CATE","ANCOM-BC2"))
out = data.frame(method=method,value=value)

p=ggplot(out, aes(x=method, y=value, fill=method))+ theme_bw()
p=p+geom_bar(stat="identity",position=position_dodge(0.8),width = 0.5)+
  scale_fill_manual(values=c('red','#8074C8','#00A08799','#E64B3599'))+
  # geom_bar(aes(col=method),stat="identity",position=position_dodge(0.75),width = 0.5)+
  # facet_grid(vars(class), labeller = labeller(.cols = label_both))
  # ylim(c(0,0.25))+ 
  labs( y = 'Jaccard')+theme(legend.position = "right",axis.text.x = element_text(hjust=1,angle=45))

p


