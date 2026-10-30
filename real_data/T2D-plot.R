#jaccard
# taxs=union(colnames(otutable1),colnames(otutable))
# 
# repro = match(intersect(colnames(otutable1),colnames(otutable)),taxs)
# repro1=length(intersect(intersect(rej4,rej5),repro))/length(intersect(union(rej4,rej5),repro))
# repro2=length(intersect(intersect(rej4.ld,rej5.ld),repro))/length(intersect(union(rej4.ld,rej5.ld),repro))
# repro3=length(intersect(intersect(rej4.ancombc,rej5.ancombc),repro))/length(intersect(union(rej4.ancombc,rej5.ancombc),repro))
# repro4=length(intersect(intersect(rej4.ancombc2,rej5.ancombc2),repro))/length(intersect(union(rej4.ancombc2,rej5.ancombc2),repro))
# repro5=length(intersect(intersect(rej4.locom,rej5.locom),repro))/length(intersect(union(rej4.locom,rej5.locom),repro))
# 



methods <- c("PoDA", "LinDA", "CATE","ANCOM-BC2")
# jaccard_index <- c(repro1, repro2, repro5, repro3, repro4)
jaccard_index = readRDS("~/weihaowang/project2/realdata/t2d_jaccard.rds")
jaccard_index=jaccard_index[-4]
jaccard_index[4]=jaccard_index[4]+0.001

data <- data.frame(Method = methods, Jaccard_Index = jaccard_index)


data$Method <- factor(data$Method, levels = c("PoDA", "LinDA", "CATE","ANCOM-BC2"))


p.sig2=ggplot(data, aes(x = Method, y = Jaccard_Index)) +
  geom_bar(stat = "identity", aes(fill = Method), position=position_dodge(0.1),width = 0.5) +
  scale_fill_manual(values = c("PoDA" = "#E47159", "LinDA" = "#3D5c6f", "CATE" = "#3D5c6f", "ANCOM-BC2" = "#3D5c6f")) +
  theme(
    legend.position = "none",
    axis.text.x = element_text(angle = 0, hjust = 0.5, vjust = 1,size = 8),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_blank(),
    axis.line = element_line(),
    axis.title.y = element_text(size = 10)
  )+
  ylab("Jaccard index")+xlab("")
p.sig2





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
rej.list=readRDS("~/weihaowang/project2/realdata/t2d_noconfun_result.rds")
# rej.list=rej.list[[-4]]
findings=matrix(0,nrow=200,ncol=5)
# rej.list=list()
# rej.list[[1]]=rej4.ld
# rej.list[[2]]=rej4
# rej.list[[3]]=rej4.locom
# rej.list[[4]]=rej4.ancombc
# rej.list[[5]]=rej4.ancombc2
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


upset_plot1
upset_plot1_wrp <- wrap_elements(full = upset_plot1)





##################################################################################################
##################################################################################################
##################################################################################################
Error1=readRDS("~/weihaowang/project2/realdata/t2dtypeIerror.rds")
modified_string <- sub(".*_", "", rownames(Y))
rejtaxon = Error1[,1]/100

sig=rep(0,length(rejtaxon))

out = data.frame(
  sig = factor(sig),
  taxon = factor(modified_string),
  levels = modified_string,
  value = rejtaxon,
  facet = "LinDA"
)

# Plot
p=ggplot(out, aes(x=taxon, y=value, fill=sig))+ theme_bw()
p.sig1 <-p +
  theme_bw() +
  geom_bar(stat = "identity", position = position_dodge(0.75), width = 0.5) +
  scale_fill_manual(values = c('#E47159', '#3D5C6F')) +
  theme(
    legend.position = "none", axis.text.y = element_text(
      angle = 90,
      hjust = 0.5,
      vjust = 0.5   ,size=7      ),
    axis.text.x = element_blank(), axis.ticks.x = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank() ,plot.margin = margin(0, 0, 0, 0),
    strip.text = element_text(size = 16)
  ) +
  ylim(c(0, 0.05)) +
  labs(y = '', x = "") +
  geom_hline(yintercept = 0.05, linetype = "dashed", color = "black") +
  facet_grid(facet ~ ., scales = "free", labeller = labeller(facet = label_value))



rejtaxon = Error1[,2]/100

sig=rep(0,length(rejtaxon))

out = data.frame(
  sig = factor(sig),
  taxon = factor(modified_string),
  levels = modified_string,
  value = rejtaxon,
  facet = "Wallen1",
  facet2 = "PoDA"
)

# Plot
p=ggplot(out, aes(x=taxon, y=value, fill=sig))+ theme_bw()
p.sig2 <-p +
  theme_bw() +
  geom_bar(stat = "identity", position = position_dodge(0.75), width = 0.5) +
  scale_fill_manual(values = c('#E47159', '#3D5C6F')) +
  theme(
    legend.position = "none",  axis.text.y = element_text(
      angle = 90,
      hjust = 0.5,
      vjust = 0.5   ,size=7      ),
    axis.text.x = element_blank(), axis.ticks.x = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank() ,plot.margin = margin(0, 0, 0, 0),
    strip.text.x = element_text(size = 22),  
    strip.text.y = element_text(size = 16)  
  ) +
  ylim(c(0, 0.05)) +
  labs(y = '', x = "") +
  geom_hline(yintercept = 0.05, linetype = "dashed", color = "black") +
  # facet_grid(facet2 ~ facet, scales = "free", labeller = labeller(facet = label_value))
  facet_grid(facet2 ~ ., scales = "free", labeller = labeller(facet = label_value))



rejtaxon = Error1[,3]/100

sig=rep(0,length(rejtaxon))

out = data.frame(
  sig = factor(sig),
  taxon = factor(modified_string),
  levels = modified_string,
  value = rejtaxon,
  facet = "CATE"
)

# Plot
p=ggplot(out, aes(x=taxon, y=value, fill=sig))+ theme_bw()
p.sig3 <-p +
  theme_bw() +
  geom_bar(stat = "identity", position = position_dodge(0.75), width = 0.5) +
  scale_fill_manual(values = c('#E47159', '#3D5C6F')) +
  theme(
    legend.position = "none",  axis.text.y = element_text(
      angle = 90,
      hjust = 0.5,
      vjust = 0.5   ,size=7      ),
    axis.text.x = element_blank(),axis.ticks.x = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank() ,plot.margin = margin(0, 0, 0, 0),
    strip.text = element_text(size = 16)
  ) +
  ylim(c(0, 0.05)) +
  labs(y = '', x = "") +
  geom_hline(yintercept = 0.05, linetype = "dashed", color = "black") +
  facet_grid(facet ~ ., scales = "free", labeller = labeller(facet = label_value))







rejtaxon = Error1[,4]/100

sig=rep(0,length(rejtaxon))

out = data.frame(
  sig = factor(sig),
  taxon = factor(modified_string),
  levels = modified_string,
  value = rejtaxon,
  facet = "ANCOM-BC2"
)

# Plot
p=ggplot(out, aes(x=taxon, y=value, fill=sig))+ theme_bw()
p.sig4 <-p +
  theme_bw() +
  geom_bar(stat = "identity", position = position_dodge(0.75), width = 0.5) +
  scale_fill_manual(values = c('#E47159', '#3D5C6F')) +
  theme(
    legend.position = "none",
    axis.ticks.y = element_blank(),axis.text.y = element_text(
      angle = 90,
      hjust = 0.5,
      vjust = 0.5   ,size=7      ),
    axis.text.x = element_text(hjust=1,angle=90,vjust=0.5,size=12,face = "bold"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank() ,plot.margin = margin(0, 0, 0, 0),
    strip.text = element_text(size = 16)
  ) +
  ylim(c(0, 0.7)) +
  labs(y = '', x = "") +
  geom_hline(yintercept = 0.05, linetype = "dashed", color = "black") +
  facet_grid(facet ~ ., scales = "free", labeller = labeller(facet = label_value))

# Display the plot
p.sig4




arranged_plots=ggarrange(p.sig2, p.sig1, p.sig3, p.sig4,
                         ncol = 1,
                         heights = c(1, 1, 1,2))



# ptype2=arranged_plots2

ptype <- annotate_figure(
  arranged_plots,
  left = text_grob("Type I error",
                   rot = 90,
                   vjust = 0.5,
                   hjust = 0,
                   size = 22)
)
ptype




