
count=T2D_QinJ[[1]]
group=T2D_QinJ[[2]]
Gamma=T2D_QinJ[[3]]

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


source("~/weihaowang/project2/PODA2confun.R")
rej2 = PODA2(as.matrix(Z),Y,type="hard-modi",formulaa,seed=1234)
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

rej.ancombc2 <- ancombc2.fun.ss(nrow(Y),Y, Z, formula, alpha=0.05,p_adj_method="BH")
rej.ancombc2

# rej.ancombc2 <- ancombc2.fun(nrow(Y),Y, Z, formula, alpha=0.05,p_adj_method="BH")
# rej.ancombc2

all_combinations <- unlist(
  lapply(1:7, function(k) combn(1:7, k, simplify = FALSE)),
  recursive = FALSE
)


jcd.all = NULL
jcd.all.ld = NULL
jcd.all.cate = NULL
jcd.all.cate.ct = NULL
jcd.all.ancombc = NULL
for(i in 1:length(all_combinations)){
  # for(i in 1:5){
index=all_combinations[[i]]
count=T2D_QinJ[[1]]
group=T2D_QinJ[[2]]
Gamma=T2D_QinJ[[3]]
Y <- t(count)
# Z = cbind(as.matrix(group),as.matrix(Gamma[,-c(2,5,8:15)]))
Gamma1=as.matrix(Gamma[,-c(2,8:15)])
Z = cbind(as.matrix(group),Gamma1[,index])
if(ncol(Z)==2){
  colnames(Z)=c("u","z1")  
}else{
  colnames(Z)=c("u",paste0("z",1:ncol(Z[,-1])))  
}
formula=paste(colnames(Z), collapse = "+")
formulaa=paste("~",formula)



######## 比较 ########
lib_cut=1000
res <- LinDA::linda(Y, as.matrix(Z), paste('~',formula),type = "count", alpha=0.05,lib.cut = lib_cut)
# res <-lindaa(Y, Z, paste('~',formula),type = "count", alpha=0.05,lib.cut = lib_cut)
rej.ld1 <- which(res$output[[1]]$reject)
rej.ld1


source("~/weihaowang/project2/PODA2confun.R")
rej1 = PODA2(as.matrix(Z),Y,type="hard-modi",formulaa)
rej1


y_implement = .imputation(Y,Z,formulaa)
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

rej.ancombc1 <- ancombc2.fun.ss(nrow(Y),Y, Z, formula, alpha=0.05,p_adj_method="BH")
rej.ancombc1

ancombc2.fun(nrow(Y),Y, Z, formula, alpha=0.05,p_adj_method="BH")

jcd = length(intersect(rej1, rej2))/length(union(rej1, rej2))
jcd.ld = length(intersect(rej.ld1, rej.ld2))/length(union(rej.ld1, rej.ld2))

if(length(union(rej.cate1, rej.cate2))==0){
  jcd.cate=0
}else{
jcd.cate = length(intersect(rej.cate1, rej.cate2))/length(union(rej.cate1, rej.cate2))
}

if(length(union(rej.cate.ct1, rej.cate.ct2))==0){
  jcd.cate.ct=0
}else{
jcd.cate.ct = length(intersect(rej.cate.ct1, rej.cate.ct2))/length(union(rej.cate.ct1, rej.cate.ct2))
}


if(length(union(rej.ancombc1, rej.ancombc2))==0){
  jcd.ancombc=0
}else{
  jcd.ancombc = length(intersect(rej.ancombc1, rej.ancombc2))/length(union(rej.ancombc1, rej.ancombc2))
}


jcd.all = c(jcd.all,jcd)
jcd.all.ld = c(jcd.all.ld,jcd.ld)
jcd.all.cate = c(jcd.all.cate,jcd.cate)
jcd.all.cate.ct = c(jcd.all.cate.ct,jcd.cate.ct)
jcd.all.ancombc = c(jcd.all.ancombc,jcd.ancombc)
  }
mean(jcd.all)
mean(jcd.all.ld)
mean(jcd.all.cate)
mean(jcd.all.cate.ct)
mean(jcd.ancombc)

t2d_jaccard=c(mean(jcd.all),mean(jcd.all.ld),mean(jcd.all.cate),mean(jcd.all.cate.ct),mean(jcd.ancombc))
saveRDS(t2d_jaccard,"~/weihaowang/project2/realdata/t2d_jaccard.rds")

t2d_noconfun_result = list(rej1,rej.ld1,rej.cate1,rej.cate.ct1,rej.ancombc1)
saveRDS(t2d_noconfun_result,"~/weihaowang/project2/realdata/t2d_noconfun_result.rds")