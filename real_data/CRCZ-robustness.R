################################################################
########################## CRC DATA ############################
################################################################
otutable=read.table(file = "~/weihaowang/project2/realdata/crc/crc_baxter.otu_table.txt",sep="\t",header=T)
metadata=read.table(file = "~/weihaowang/project2/realdata/crc/crc_baxter.metadata.txt",sep="\t",header=T)
otutable1=read.table(file = "~/weihaowang/project2/realdata/crc/crc_zackular.otu_table.txt",sep="\t",header=T)
metadata1=read.table(file = "~/weihaowang/project2/realdata/crc/crc_zackular.metadata.txt",sep="\t",header=T, fileEncoding="ISO-8859-1")



tax <- as.data.frame(strsplit(otutable[,1],";"))
tax1 <- as.data.frame(strsplit(otutable1[,1],";"))
tax <- t(tax)
tax1 <- t(tax1)
rownames(tax) <- otutable[,1]
rownames(tax1) <- otutable1[,1]
rownames(otutable) <- otutable[,1]
rownames(otutable1) <- otutable1[,1]


otutable <- otutable[,-1]
colnames(otutable) <- as.data.frame(strsplit(colnames(otutable),"X"))[2,]
otutable1 <- otutable1[,-1]


case.id <- metadata[which(metadata$DiseaseState=="CRC"),1]
control.id <- metadata[which(metadata$DiseaseState=="H"),1]
data.id <- c(case.id,control.id)
group = c(rep(1,length(case.id)),rep(0,length(control.id)))
group <- group[which(data.id%in%colnames(otutable))]
sex = metadata[match(data.id,colnames(otutable))[which(data.id%in%colnames(otutable))],12]
sex = as.numeric(factor(gsub("[^a-zA-Z]", "", sex)))
age = metadata[match(data.id,colnames(otutable))[which(data.id%in%colnames(otutable))],3]
age = as.numeric(gsub("[^0-9]", "", age))
bmi = metadata[match(data.id,colnames(otutable))[which(data.id%in%colnames(otutable))],5] #有missing
otutable <- otutable[,match(data.id,colnames(otutable))[which(data.id%in%colnames(otutable))]]




case.id1 <- metadata1[which(metadata1$DiseaseState=="CRC"),1]
control.id1 <- metadata1[which(metadata1$DiseaseState=="H"),1]
data.id1 <- c(case.id1,control.id1)
data.id1 <- gsub("-", ".", data.id1)
group1 = c(rep(1,length(case.id1)),rep(0,length(control.id1)))
group1 <- group1[which(data.id1%in%colnames(otutable1))]
sex1 = metadata1[match(data.id1,colnames(otutable1))[which(data.id1%in%colnames(otutable1))],34]
sex1 = as.numeric(factor(gsub("[^a-zA-Z]", "", sex1)))
age1 = metadata1[match(data.id1,colnames(otutable1))[which(data.id1%in%colnames(otutable1))],33]
age1 = as.numeric(gsub("[^0-9]", "", age1))
# substr(colnames(otutable),2,8)
otutable1 <- otutable1[,match(data.id1,colnames(otutable1))[which(data.id1%in%colnames(otutable1))]]


tax <- tax[,-(7:8)]
taxname <- phyloseq(tax_table(tax))
taxa_names(taxname) <- rownames(otutable)
otutable=t(as(otutable,"matrix"))


phy_filter <- phyloseq(tax_table(taxname), otu_table((otutable),taxa_are_rows=F))
phy.glom <- tax_glom(phy_filter, taxrank="ta6")#, NArm=F
otu.data.f <- phy.glom@otu_table
datasize <- nrow(otu.data.f)


prevalence = apply(as(otu.data.f, "matrix"), 2, function(x) {
  return(sum(x > 0))
})/(datasize)

keepOTUs = prevalence>  0.1

otu.data.f1=prune_taxa(keepOTUs, otu.data.f)
count=otu.data.f1@.Data





case.id1 <- metadata1[which(metadata1$DiseaseState=="CRC"),1]
control.id1 <- metadata1[which(metadata1$DiseaseState=="H"),1]
data.id1 <- c(case.id1,control.id1)
data.id1 <- gsub("-", ".", data.id1)
group1 = c(rep(1,length(case.id1)),rep(0,length(control.id1)))
group1 <- group1[which(data.id1%in%colnames(otutable1))]
sex1 = metadata1[match(data.id1,colnames(otutable1))[which(data.id1%in%colnames(otutable1))],34]
sex1 = as.numeric(factor(gsub("[^a-zA-Z]", "", sex1)))
age1 = metadata1[match(data.id1,colnames(otutable1))[which(data.id1%in%colnames(otutable1))],33]
age1 = as.numeric(gsub("[^0-9]", "", age1))
# substr(colnames(otutable),2,8)
height1 = metadata1[match(data.id1,colnames(otutable1))[which(data.id1%in%colnames(otutable1))],37] #也是有很多确实
mass1 = metadata1[match(data.id1,colnames(otutable1))[which(data.id1%in%colnames(otutable1))],38] #有缺失！
otutable1 <- otutable1[,match(data.id1,colnames(otutable1))[which(data.id1%in%colnames(otutable1))]]

tax1 <- tax1[,-(7:8)]
taxname1 <- phyloseq(tax_table(tax1))
taxa_names(taxname1) <- rownames(otutable1)

otutable1=t(as(otutable1,"matrix"))



phy_filter <- phyloseq(tax_table(taxname1), otu_table((otutable1),taxa_are_rows=F))
phy.glom <- tax_glom(phy_filter, taxrank="ta6")#, NArm=F
otu.data.f <- phy.glom@otu_table
datasize <- nrow(otu.data.f)

prevalence = apply(as(otu.data.f, "matrix"), 2, function(x) {
  return(sum(x > 0))
})/(datasize)

keepOTUs = prevalence>  0.1

otu.data.f1=prune_taxa(keepOTUs, otu.data.f)
count1=otu.data.f1@.Data



tax_g <- as.data.frame(strsplit(colnames(count),";"))
tax1_g <- as.data.frame(strsplit(colnames(count1),";"))
tax_g <- t(tax_g)[,-(7:8)]
tax1_g <- t(tax1_g)[,-(7:8)]
colnames(count)=apply(tax_g, 1, paste, collapse = ";")
colnames(count1)=apply(tax1_g, 1, paste, collapse = ";")


# common <- intersect(colnames(count), colnames(count1))
# count <- count[,match(common,colnames(count))]
# count1 <- count1[,match(common,colnames(count1))]
CRC_baxter=list(count,group,age,sex)
CRC_zackular=list(count1,group1,age1,sex1)

taxs=union(colnames(count),colnames(count1))

CRC_baxter <- readRDS("~/weihaowang/project2/realdata/crc/CRC_baxter.rds")
CRC_zackular <- readRDS("~/weihaowang/project2/realdata/crc/CRC_zackular.rds")



################################################################
####################### NMDS PLOT ##############################
################################################################
set.seed(1111)
otu=rbind(CRC_baxter[[1]],CRC_zackular[[1]])
otu=otu/rowSums(otu)
otu %>%data.frame() ->otu
group=c(CRC_baxter[[2]],CRC_zackular[[2]]+2)
group[group==0]="Baxter:Healthy"
group[group==1]="Baxter:CRC"
group[group==2]="Zackular:Healthy"
group[group==3]="Zackular:CRC"
distance <- vegdist(otu, method = 'bray')
dfNmds<-metaMDS(otu,distance="bray",k = 2)
plot_data <- data.frame(dfNmds$points)
plot_data$group=as.factor(group)
names(plot_data)[1:2] <- c('NMDS1', 'NMDS2')
centroid1 <- aggregate(cbind(NMDS1,NMDS2) ~ group,
                       data = plot_data,
                       FUN = mean)
centroid1$reference=rep(c("Baxter","Zackular"),each=2)
centroid1$X=rep(c("case","control"),2)
dist(centroid1[,2:3],p=2)
p1=ggplot(data = centroid1, aes(x=NMDS1, y=NMDS2))+
  geom_point(aes(x=NMDS1, y=NMDS2,color=as.factor(reference)),size=5) +
  geom_text_repel(aes(label = group))+
  xlim(c(-0.08,0.18))+ylim(c(-0.08,0.18))+
  geom_line(aes(group = X),color= "#5A9599FF")+
  geom_line(aes(group = reference),color= c("#FF6348FF"))+
  scale_color_manual(values = c("#D55E0080","#91D1C2FF"))+
  labs(title = expression("CRC"))+
  theme_bw()+theme(plot.title = element_text(hjust=0.5),legend.position = "none")
p1




################################################################
############### DA DETECTION AND CONSISTENCY ###################
################################################################
lib_cut=0


Y=t(CRC_baxter[[1]])
u=CRC_baxter[[2]]

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

rej2=match(rownames(Y)[rej],taxs)
rej2.ld=match(rownames(Y)[rej.ld],taxs)
rej2.cate=match(rownames(Y)[rej.cate],taxs)
rej2.cate.ct=match(rownames(Y)[rej.cate.ct],taxs)
rej2.ancombc=match(rownames(Y)[rej.ancombc],taxs)

################################################################
############### DA DETECTION AND CONSISTENCY ###################
################################################################
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

rej3=match(rownames(Y)[rej],taxs)
rej3.ld=match(rownames(Y)[rej.ld],taxs)
rej3.cate=match(rownames(Y)[rej.cate],taxs)
rej3.cate.ct=match(rownames(Y)[rej.cate.ct],taxs)
rej3.ancombc=match(rownames(Y)[rej.ancombc],taxs)


inter=length(intersect(rej3,rej2))
consist1=inter/length(union(rej3,rej2))
consist1

inter=length(intersect(rej3.ld,rej2.ld))
consist2=inter/length(union(rej3.ld,rej2.ld))
consist2

inter=length(intersect(rej3.cate,rej2.cate))
consist3=inter/length(union(rej3.cate,rej2.cate))
consist3

inter=length(intersect(rej3.cate.ct,rej2.cate.ct))
consist4=inter/length(union(rej3.cate.ct,rej2.cate.ct))
consist4

inter=length(intersect(rej3.ancombc,rej2.ancombc))
consist5=inter/length(union(rej3.ancombc,rej2.ancombc))
consist5



value=c(consist1,consist2,consist3,consist5)



# sig=rep(0,5)
# sig[match(rej,Nind)]=1
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



##################################################################################################
##################################################################################################
##################################################################################################
library(patchwork)
library(tidyr)
library(dplyr)
library(tibble)
library(ggplot2)
library(ggpubr)
library(stringr)
library(ComplexUpset)
rm(count)
# rej.list=readRDS("~/weihaowang/project2/realdata/t2d_noconfun_result.rds")
# rej.list=rej.list[[-4]]
findings=matrix(0,nrow=200,ncol=5)
rej.list=list()
rej.list[[1]]=rej2
rej.list[[2]]=rej2.ld
rej.list[[3]]=rej2.cate
rej.list[[4]]=rej2.cate.ct
rej.list[[5]]=rej2.ancombc
for(i in 1:5){
  findings[rej.list[[i]],i]=1
}
findings=findings[,-4]


findings = as.data.frame(findings)
colnames(findings)=c("PoDA", "LinDA", "CATE","ANCOM-BC2")

findings=findings[which(rowSums(findings)!=0),]


bar_colors <- c("#3D5c6f", rep("blue", 10))

upset_plot1 <-ComplexUpset::upset(
  findings,
  intersect = c(c("PoDA", "LinDA", "CATE","ANCOM-BC2")),
  name = NULL,
  queries=list(
    # ComplexUpset::upset_query(
    #   intersect=c('PoDA'),
    #   color="#E47159",
    #   fill="#E47159",
    #   only_components=c('intersections_matrix', 'Intersection size')
    # ),
    ComplexUpset::upset_query("LinDA", fill="#F6BCA9"),
    ComplexUpset::upset_query("CATE", fill="#F6BCA9"),
    ComplexUpset::upset_query("ANCOM-BC2", fill="#F6BCA9"),
    ComplexUpset::upset_query("PoDA", fill="#E47159")
  ),
  base_annotations = list(
    "Intersection size" = ComplexUpset::intersection_size(
      counts = TRUE,
      aes(fill = as.factor(seq_along(after_stat(count))))
    ) + scale_fill_manual(values = bar_colors, guide = 'none')
  )
)


upset_plot1
upset_plot1_wrp <- wrap_elements(full = upset_plot1)




findings=matrix(0,nrow=200,ncol=5)
rej.list=list()
rej.list[[1]]=rej3
rej.list[[2]]=rej3.ld
rej.list[[3]]=rej3.cate
rej.list[[4]]=rej3.cate.ct
rej.list[[5]]=rej3.ancombc
for(i in 1:5){
  findings[rej.list[[i]],i]=1
}
findings=findings[,-4]
# findings[200,]=findings[200,]+1

findings = as.data.frame(findings)
colnames(findings)=c("PoDA", "LinDA", "CATE","ANCOM-BC2")

findings=findings[which(rowSums(findings)!=0),]


bar_colors <- c("#3D5c6f", rep("blue", 10))

upset_plot2 <-ComplexUpset::upset(
  findings,
  intersect = c(c("PoDA", "LinDA", "CATE","ANCOM-BC2")),
  name = NULL,
  queries=list(
    ComplexUpset::upset_query(
      intersect=c('PoDA'),
      color="#E47159",
      fill="#E47159",
      only_components=c('intersections_matrix', 'Intersection size')
    ),
    ComplexUpset::upset_query("LinDA", fill="#F6BCA9"),
    ComplexUpset::upset_query("CATE", fill="#F6BCA9"),
    ComplexUpset::upset_query("ANCOM-BC2", fill="#F6BCA9"),
    ComplexUpset::upset_query("PoDA", fill="#E47159")
  ),
  base_annotations = list(
    "Intersection size" = ComplexUpset::intersection_size(
      counts = TRUE,
      aes(fill = as.factor(seq_along(after_stat(count))))
    ) + scale_fill_manual(values = bar_colors, guide = 'none')
  )
)


upset_plot2
upset_plot2_wrp <- wrap_elements(full = upset_plot2)




