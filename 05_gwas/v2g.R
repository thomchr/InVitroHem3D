source(file.path(Sys.getenv("EHT_REPO", "."), "config", "paths.R"))  # DATA_DIR, OUT_DIR, REPO_DIR
setwd(file.path(DATA_DIR, "V2G"))
library(tibble)
library(reshape2)
library(ggplot2)
library(patchwork)
library(nVennR)
library(UpSetR)
library(eulerr)
library(mixOmics)
library(pathfindR)
library(clusterProfiler)
library(DESeq2)
library(ggrepel)
library(ggpubr)
library(ggnewscale)
library(clustifyr)

# re V2g #########
setwd(file.path(DATA_DIR, "V2G"))
files=list.files("blood_traits2_v2g")
allv2g=NULL
for(file in files){
  tmp=read.table(paste0("blood_traits2_v2g/",file),header=T)
  colnames(tmp)[which(colnames(tmp)=="source")]="Celltype"
  allv2g = rbind(allv2g,tmp)
}
rm(tmp,files,file)

allv2g = allv2g[-grep("close",allv2g$type),]
allv2g$ensembl_gene_id = allv2g$gene_id

tmp=unlist(strsplit(paste(allv2g$gene_id[grep("ENSG",allv2g$gene_id)]),split="\\."))
allv2g$ensembl_gene_id[grep("ENSG",allv2g$gene_id)] = paste(tmp[seq(1,by=2,length(tmp))])


allv2g = unique(allv2g[,-which(colnames(allv2g) %in% c("tx_id","gene_id"))])

for(i in which(allv2g$Celltype=="Hematopoietic_ChrisThom_multiome")){
  gen=setdiff(unique(allv2g$ensembl_gene_id[which(allv2g$gene_name == allv2g$gene_name[i])]),"NA")
  allv2g$ensembl_gene_id[i] = ifelse(length(gen)==1,gen,NA)
}

gene=unique(allv2g$gene_name[which(is.na(allv2g$ensembl_gene_id))])
library(biomaRt)
ensembl <- useEnsembl(biomart = "genes",dataset = "hsapiens_gene_ensembl")
gen <- getBM(
  attributes = c("ensembl_gene_id","hgnc_symbol",'chromosome_name','start_position','end_position'), 
  values = gene,filters = "hgnc_symbol",
  mart = ensembl
)
for(i in 1:length(gene)){
  tmp=gen[which(gen$hgnc_symbol == gene[i]),]
  if(nrow(tmp)==1){
    allv2g$ensembl_gene_id[which(allv2g$gene_name == gene[i])] = tmp$ensembl_gene_id
  }
}

# match proxy with traits #######
ld=read.table(file.path(DATA_DIR, "V2G/leads/eachStudy.leads.ld"),header=T)
traits=c("PLT","RBC","WBC","Stroke","VTE","CAD","Height","BMI","WHR","SBP","DBP","BIP","SCZ","ALZ","SmkInit","DrnkWk","T2D","T1D","Asthma")
ld = unique(ld[,c('rsID1','rsID2')])
add=NULL
for(trait in traits){
  file=list.files(file.path(DATA_DIR, "V2G/leads"),pattern = paste0(trait,".leads"))
  if(length(file)>0){
    tmp=read.table(paste0(file.path(DATA_DIR, "V2G/leads"),file))
    tmp$trait=trait
    add=rbind(add,tmp)
  }
}
tmp=data.frame(rsID1 = unique(add$V1) ,rsID2=unique(add$V1) )
ld=rbind(ld,tmp)
ld=merge(ld,add,by.x="rsID1",by.y="V1")
add=NULL
tmp=read.table(file.path(DATA_DIR, "V2G/leads/Height.bed"))
t=unlist(strsplit(tmp$V4,split="rsid:"))
tmp$V1 = paste(t[seq(2,by=2,length(t))])
tmp$trait="Height"
add=unique(rbind(add,tmp[,c("V1","trait")]))
tmp=read.table(file.path(DATA_DIR, "V2G/leads/T1D_credibles.bed"))
t=unlist(strsplit(tmp$V4,split="proxy_rsid:"))
tmp$V1 = paste(t[seq(2,by=2,length(t))])
tmp$trait="T1D"
add=unique(rbind(add,tmp[,c("V1","trait")]))
colnames(add)[1] = colnames(ld)[2] = "proxy_rsid"
add$rsID1=add$proxy_rsid
ld=unique(rbind(ld,add[,colnames(ld)])  )
add=unique(ld[,c('rsID1','trait')])
add$proxy_rsid = add$rsID1
ld=unique(rbind(ld,add[,colnames(ld)])  )

rm(add,t,tmp)
proxy=reshape2::dcast(proxy_rsid~trait,data=unique(ld[,c('proxy_rsid','trait')]))
proxy$numtrait = apply(proxy[,traits],1,function(x) length(na.exclude(x)))
table(proxy$numtrait)
proxy$trait = apply(proxy[,traits],1,function(x) paste(na.exclude(x),collapse = ","))
table(proxy$trait)

allv2g = unique(merge(allv2g,ld[,c('proxy_rsid','trait')],by.x="variant_id",by.y="proxy_rsid",all.x=TRUE))
allv2g$Celltype[which(allv2g$Celltype == "Hematopoietic_ChrisThom")] = "poolHiC_consensus"
allv2g$Celltype[which(allv2g$Celltype == "Hematopoietic_ChrisThom_multiome")] = "eRegulon_Scenicplus"

saveRDS(allv2g,file=file.path(DATA_DIR, "V2G/allv2g.rds"))
saveRDS(allv2g[-which(allv2g$Celltype %in% c("poolHiC_consensus","eRegulon_Scenicplus")),],file=file.path(DATA_DIR, "V2G/allv2g_10celltypes.rds"))
saveRDS(allv2g[which(allv2g$Celltype =="poolHiC_consensus"),],file=file.path(DATA_DIR, "V2G/allv2g_poolHiC_consensus.rds"))
saveRDS(allv2g[which(allv2g$Celltype =="eRegulon_Scenicplus"),],file=file.path(DATA_DIR, "V2G/allv2g_eRegulon_Scenicplus.rds"))

# output for paper ###########
library(nVennR)
allv2g=readRDS("allv2g.rds")
scenicv2g = allv2g[which(allv2g$Celltype =="eRegulon_Scenicplus"),]
tradgene = readLines(file.path(DATA_DIR, "hiC/blood.genes.tmp"))
eReg = readRDS(file.path(DATA_DIR, "multiome/region_gene_nth.rds"))

l=list(HiC_gene=tradgene,
       HiC_V2G_gene = unique(allv2g$gene_name[-which(allv2g$Celltype %in% c("poolHiC_consensus","eRegulon_Scenicplus"))]),
       scenic_gene=unique(eReg$Gene),
       scenic_V2G_gene=unique(scenicv2g$gene_name))
a=plotVenn(l,
           labelRegions = T,opacity = 0.3,fontScale=2,
           outFile="olap_genes_4ways.svg")

core=a$reg[[length(a$reg  )]]
saveRDS(core,file="core1908.rds")
allv2g=allv2g[which(allv2g$Celltype %in% c("poolHiC_consensus","eRegulon_Scenicplus")),]
allv2g=allv2g[which(allv2g$gene_name %in% core),]

novel=read.csv("Novel_SNP_to_Gene_Pairs.core1908.csv")
novel$FoundInOpenTarget=TRUE
allv2g = merge(allv2g,novel,by=colnames(novel)[1:4],all.x=TRUE)
ab=readRDS(file.path(DATA_DIR, "hiC/ABcompartment/gene_compartment_df.rds"))
colnames(ab)[5]="Expression_Log2FC"
allv2g = merge(allv2g,ab[,c('Gene_Name','Compartment_Change','Expression_Log2FC')],by.x = "gene_name",by.y="Gene_Name",all.x=TRUE)

ab=readRDS(file.path(DATA_DIR, "hiC/ABcompartment/region_compartment_df.rds"))
allv2g$RegionName=paste0(allv2g$chr,"-",allv2g$ocr_start,"-",allv2g$ocr_end)
colnames(ab)[5]="Accessibility_Log2FC"
allv2g = merge(allv2g,ab[,c('RegionName','Accessibility_Log2FC')],by = "RegionName",all.x=TRUE)

ab=readRDS(file.path(DATA_DIR, "hiC/compareTAD/RNA_TAD_changes.rds"))
ab = unique(ab[,1:2])
colnames(ab)[2]="TADchanges"
allv2g = merge(allv2g,ab,by.x = "gene_name",by.y="Gene_Name",all.x=TRUE)
ab=readRDS(file.path(DATA_DIR, "hiC/compareTAD/RNA_subTAD_changes.rds"))
ab = unique(ab[,1:2])
colnames(ab)[2]="subTADchanges"
allv2g = merge(allv2g,ab,by.x = "gene_name",by.y="Gene_Name",all.x=TRUE)

pair=reshape2::dcast(gene_name+ensembl_gene_id+RegionName+chr+variant_id~'traits',data=allv2g,value.var = "trait",fun.aggregate = function(x) paste(unique(x),collapse = ", "))

allv2g=unique(allv2g[,c("gene_name","ensembl_gene_id","RegionName","chr",
                        "variant_id","variant_pos","Expression_Log2FC","Accessibility_Log2FC",
                        "FoundInOpenTarget","Compartment_Change","TADchanges","subTADchanges" )])
allv2g = merge(pair,allv2g,by=intersect(colnames(pair),colnames(allv2g)))

write.table(allv2g,file="V2G_1908.txt",col.names = T,row.names = F,sep="\t",quote=F)


