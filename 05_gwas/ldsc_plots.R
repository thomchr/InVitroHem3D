source(file.path(Sys.getenv("EHT_REPO", "."), "config", "paths.R"))  # DATA_DIR, OUT_DIR, REPO_DIR
setwd(file.path(DATA_DIR, "ldsc"))
library(ggplot2)
library(ggsci)
library(egg)
library(scales)
library(reshape2)
library(patchwork)
library(paletteer)

options(scipen = 100)

traits=c("PLT", "RBC", "WBC" ,"Stroke", "VTE", "CAD", "Height" ,"BMI" ,"WHR" ,
         "SBP", "DBP", "BIP", "SCZ", "ALZ" ,"SmkInit" ,"DrnkWk", "T2D", "T1D" ,"Asthma")
N=c(1186350,1186390,1186574,1167058,1163612,1177354,1073504,1175763,1175784,1158394,
  1160153,1166522,1168073,1167687,1181366,1181159,1186298,1185028,967201)
# ldsc panel ###########
allstat=read.delim(file.path(DATA_DIR, "Celltypes_repertoire.txt"),header = T)
ipsblod = read.delim(file.path(DATA_DIR, "ldsc/Hematopoietic_ChrisThom_list.txt"))
ipsblod$file = factor(ipsblod$file,levels = paste(ipsblod$file))
g1=ggplot(ipsblod)+geom_bar(aes(atac_peaks_N,file),fill="#99c3ff",stat = "identity",width = .8)+
  geom_bar(aes(ocr_N,file),fill="#f3776a",stat = "identity",width = .7)+
  theme_bw()+
  scale_y_discrete(name="Cell types")+
  scale_x_continuous(expand = c(0,0),name = "Number of open chromatin regions",
                     breaks = seq(0,210000,by=100000),
                     labels = label_number(suffix = " K", scale = 1e-3),
                     limits = c(0,210000))+
  theme(panel.grid.minor = element_blank(),
        legend.position = c(0.75,0.9),legend.background = element_blank())

ALL=NULL
for(i in 1:length(traits)){
  tmp = read.delim(paste0(file.path(DATA_DIR, "ldsc/narrow/summary_"),traits[i],".txt"),header=T,sep="\t",fill = TRUE)
  tmp$trait = paste(traits[i])
  tmp$cate = "cREs"
  tmp$N = tmp$Prop._SNPs * N[i]
  ALL=rbind(ALL,tmp)
  
  tmp = read.delim(paste0(file.path(DATA_DIR, "ldsc/atac/summary_"),traits[i],".txt"),header=T,sep="\t",fill = TRUE)
  tmp$trait = paste(traits[i])
  tmp$cate = "Total OCRs"
  tmp$N = tmp$Prop._SNPs * N[i]
  ALL=rbind(ALL,tmp)
  
  tmp = read.delim(paste0(file.path(DATA_DIR, "ldsc/Prom/summary_"),traits[i],".txt"),header=T,sep="\t",fill = TRUE)
  tmp$trait = paste(traits[i])
  tmp$cate = "Promoter OCRs"
  tmp$N = tmp$Prop._SNPs * N[i]
  ALL=rbind(ALL,tmp)
  
  tmp = read.delim(paste0(file.path(DATA_DIR, "ldsc/expand/summary_"),traits[i],".txt"),header=T,sep="\t",fill = TRUE)
  tmp$trait = paste(traits[i])
  tmp$cate = "cREs ±500bps"
  tmp$N = tmp$Prop._SNPs * N[i]
  ALL=rbind(ALL,tmp)
  
  tmp = read.delim(paste0(file.path(DATA_DIR, "ldsc/notcRE/summary_"),traits[i],".txt"),header=T,sep="\t",fill = TRUE)
  tmp$trait = paste(traits[i])
  tmp$cate = "not_cRE\nOCRs"
  tmp$N = tmp$Prop._SNPs * N[i]
  ALL=rbind(ALL,tmp)
  
  tmp = read.delim(paste0(file.path(DATA_DIR, "ldsc/notcRE_notProm/summary_"),traits[i],".txt"),header=T,sep="\t",fill = TRUE)
  tmp$trait = paste(traits[i])
  tmp$cate = "not_cRE/Prom\nOCRs"
  tmp$N = tmp$Prop._SNPs * N[i]
  ALL=rbind(ALL,tmp)
}
rm(i,tmp)

ALL$sig = round(ALL$Enrichment_p,digits = 2) <= 0.05

allbl=NULL
for(i in 1:nrow(ipsblod)){
  a = ALL[which(ALL$Celltype==ipsblod$file[i]),]
  
  for(e in 1:length(traits)){
    b = a[which(a$trait==traits[e]),]
    cre=b[which(b$cate=="cREs"),]
    b = b[-which(b$cate=="cREs"),]
    b$pval=NA
    for(x in 1:5){
      b$pval[x] = BSDA::tsum.test(mean.x = cre$Enrichment,s.x = cre$Enrichment_std_error,n.x = round(cre$N),
                                  mean.y = b$Enrichment[x],s.y = b$Enrichment_std_error[x],n.y = round(b$N[x]),
                                  alternative = "greater")$p.value
    }
    cre$pval=0
    allbl=rbind(allbl,cre,b)
  }
} 
allbl$sig = round(allbl$Enrichment_p,digits = 2) <= 0.05
allbl$tsig = allbl$Enrichment_p < 0.05 & allbl$cate!="cREs"
allbl$cate = factor(allbl$cate,levels = c("Total OCRs",
                                      "Promoter OCRs",
                                      "cREs",
                                      "cREs ±500bps",
                                      "not_cRE\nOCRs",
                                      "not_cRE/Prom\nOCRs"))
allbl$Celltype=factor(allbl$Celltype,levels = levels(ipsblod$file))
minbin = round(min(allbl$Prop._SNPs),2)
maxbin = round(max(allbl$Prop._SNPs),2)
minp = floor(min(-log10(allbl$Enrichment_p)))
maxp = ceiling(max(-log10(allbl$Enrichment_p)))
maxenr = ceiling(max(allbl$Enrichment[allbl$cate=="cREs"]))
minenr = floor(min(allbl$Enrichment[allbl$cate=="cREs"]))
allbl$trait=factor(allbl$trait,levels = traits)
g2=ggplot(data=allbl[allbl$cate=="cREs",])+geom_vline(xintercept = 1, linetype='dashed', linewidth=.5)+
  geom_point(aes(Enrichment,factor(Celltype),size=`Prop._SNPs`,color=-log10(Enrichment_p)),alpha=0.7)+
  geom_errorbarh(aes(y=Celltype,xmax=Enrichment+Enrichment_std_error,xmin=Enrichment-Enrichment_std_error),color="#FF7F0EFF",height=.5,linewidth=.4,position = "dodge")+
  theme_bw()+
  geom_point(data=allbl[which(allbl$Enrichment_p < 0.05 & allbl$cate=="cREs"),],aes(Enrichment,Celltype),shape=8,color="white",size=1)+
  scale_size_binned(name = "Proportion of SNP",breaks = seq(0,maxbin,length.out=4),
                    labels = function(x) paste0(as.character(round(x,3)*100),"%"),limits = c(0,maxbin))+
  scale_x_continuous(name = "Fold enrichment of heritability")+
  scale_color_gradient(name="-log10 p-value",high = "#B2182B",low = "blue",limits = c(minp,maxp))+
  facet_grid(cols = vars(trait),scales = "free",space = "free_y",drop = TRUE)+
  theme(axis.title.y = element_blank(),axis.text.y = element_blank(),plot.caption = element_text(hjust = 0.5,vjust = 0),
        axis.title.x.top = element_blank(),axis.text.x.top = element_blank(),axis.ticks.x.top = element_blank(),
        panel.grid.minor = element_blank(),plot.title = element_text(hjust = 0.5),
        strip.text.y = element_blank(),legend.position = "right")+
  guides(size = guide_bins(show.limits = TRUE))
pdf(file.path(DATA_DIR, "ldsc/cREs_ldsc.pdf"),height = 4,width = 20)
g1+g2+plot_layout(ncol = 2,widths = c(1.5,19),guides = "collect")
dev.off()

