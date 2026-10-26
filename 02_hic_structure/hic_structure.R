# res 40K HiCRep #####
source(file.path(Sys.getenv("EHT_REPO", "."), "config", "paths.R"))  # DATA_DIR, OUT_DIR, REPO_DIR
scv=read.csv(file.path(DATA_DIR, "hiC/scc_heatmap.40K.csv"))
scv = rbind(scv,data.frame(sample1=c("HE_D9_plus7_rep1","HE_D9_plus7_rep2","HPC_Pro_D9_plus8_rep1","HPC_Pro_D9_plus8_rep2"),
                           sample2=c("HE_D9_plus7_rep1","HE_D9_plus7_rep2","HPC_Pro_D9_plus8_rep1","HPC_Pro_D9_plus8_rep2"),scc=1))

mat=reshape2::dcast(sample1~sample2,data=scv,value.var = "scc")
rownames(mat)=mat$sample1
mat=as.matrix(mat[,-1])
mat[upper.tri(mat)] = mat[lower.tri(mat)] 
library(gplots)
library(RColorBrewer)
pdf(file.path(DATA_DIR, "hiC/scc_heatmap.40K.pdf"))
heatmap.2(
  mat,
  distfun = function(x) as.dist((1-x)/2),
  hclust=function(x) hclust(x,method="complete"),
  Rowv = TRUE, Colv= T,
  col=brewer.pal(9,"Blues"),
  symm = T,
  margins = c(12,12),
  trace="none",
  density.info="none",
  key.title = NA,
  keysize=1,
  key.xlab="SCC"
)
dev.off()
rm(scv,mat)

# TAD compare #########
library(ggplot2)
library(scales)
library(smplot2)
library(forcats)
library(ggsci)
library(ggthemes)
library(paletteer)
library(dplyr)

HE=read.table(file.path(DATA_DIR, "hiC/compareTAD/HE_D9_plus7.hiTADs_10K.bed"))
HPC_Pro=read.table(file.path(DATA_DIR, "hiC/compareTAD/HPC_Pro_D9_plus8.hiTADs_10K.bed"))
colnames(HE) = colnames(HPC_Pro) = c("chr",'start','end','tadlevel')
HE$he_domain = paste(HE$chr,HE$start,HE$end,sep = ".")
HPC_Pro$hpc_domain = paste(HPC_Pro$chr,HPC_Pro$start,HPC_Pro$end,sep=".")


pdf(file.path(DATA_DIR, "hiC/compareTAD/numberTAD.pdf"),height = 3,width = 3)
ggplot()+
  geom_bar(data=HE,aes(x="HE",group=tadlevel,fill=factor(tadlevel)),stat = "count",position = "stack")+
  geom_bar(data=HPC_Pro,aes(x="HPC_Pro",group=tadlevel,fill=factor(tadlevel)),stat = "count",position = "stack")+
  scale_fill_colorblind(name="Hierarchical\nTAD level",breaks=seq(6,0))+
  scale_y_continuous("TADs and sub-domains number x1K",labels = label_number(scale = 1e-3) ,expand = c(0,500) )+
  annotate("text",x=c("HE","HPC_Pro"),y=c(nrow(HE)+10,nrow(HPC_Pro)+10),label=c(nrow(HE),nrow(HPC_Pro)))+
  labs(x="")+
theme_minimal()
dev.off()



conservedTAD=read.table(file.path(DATA_DIR, "hiC/compareTAD/domain_conserved.txt"))
conservedTAD$he_domain = paste(conservedTAD$V1,conservedTAD$V2,conservedTAD$V3,sep=".")
conservedTAD$hpc_domain = paste(conservedTAD$V4,conservedTAD$V5,conservedTAD$V6,sep=".")
HE$conservedTAD = HE$he_domain %in% conservedTAD$he_domain
HPC_Pro$conservedTAD = HPC_Pro$hpc_domain %in% conservedTAD$hpc_domain

semi=read.table(file.path(DATA_DIR, "hiC/compareTAD/domain_inner_changed.txt"))
semi$he_domain = paste(semi$V1,semi$V2,semi$V3,sep=".")
semi$hpc_domain = paste(semi$V4,semi$V5,semi$V6,sep=".")
HE$semiTAD = HE$he_domain %in% semi$he_domain
HPC_Pro$semiTAD = HPC_Pro$hpc_domain %in% semi$hpc_domain

conservedsub=read.table(file.path(DATA_DIR, "hiC/compareTAD/subTAD_conserved.txt"))
conservedsub$he_domain = paste(conservedsub$V1,conservedsub$V2,conservedsub$V3,sep=".")
conservedsub$hpc_domain = paste(conservedsub$V4,conservedsub$V5,conservedsub$V6,sep=".")
conservedsub = merge(conservedsub,HE[,c('he_domain','tadlevel')],by="he_domain")
conservedsub = merge(conservedsub,HPC_Pro[,c('hpc_domain','tadlevel')],by="hpc_domain")
colnames(conservedsub)[9:10] = c("HE_tad","HPC_tad")
HE$conservedsubTAD = HE$he_domain %in% conservedsub$he_domain[which(conservedsub$HE_tad==conservedsub$HPC_tad)]
HPC_Pro$conservedsubTAD = HPC_Pro$hpc_domain %in% conservedsub$hpc_domain[which(conservedsub$HE_tad==conservedsub$HPC_tad)]
HE$subTADlower = HE$he_domain %in% conservedsub$he_domain[which(conservedsub$HE_tad > conservedsub$HPC_tad)]
HPC_Pro$subTADlower = HPC_Pro$hpc_domain %in% conservedsub$hpc_domain[which(conservedsub$HE_tad > conservedsub$HPC_tad)]
HE$subTADhigher = HE$he_domain %in% conservedsub$he_domain[which(conservedsub$HE_tad < conservedsub$HPC_tad)]
HPC_Pro$subTADhigher = HPC_Pro$hpc_domain %in% conservedsub$hpc_domain[which(conservedsub$HE_tad < conservedsub$HPC_tad)]

splited=read.table(file.path(DATA_DIR, "hiC/compareTAD/domain_split.txt"))
splited$he_domain = paste(splited$V1,splited$V2,splited$V3,sep=".")
splited$hpc_domain = paste(splited$V4,splited$V5,splited$V6,sep=".")
HE$splitTAD = HE$he_domain %in% splited$he_domain
HPC_Pro$splitTAD = FALSE
for(i in 1:nrow(splited)){
  HPC_Pro$splitTAD[which(HPC_Pro$chr == splited$V4[i] & HPC_Pro$start >= splited$V5[i] & HPC_Pro$end <= splited$V6[i] & HPC_Pro$tadlevel==0 )] = TRUE
}

subsplit = read.table(file.path(DATA_DIR, "hiC/compareTAD/subTAD_split.txt"))
subsplit$he_domain = paste(subsplit$V1,subsplit$V2,subsplit$V3,sep=".")
subsplit$hpc_domain = paste(subsplit$V4,subsplit$V5,subsplit$V6,sep=".")
HE$splitsubTAD = HE$he_domain %in% subsplit$he_domain
HPC_Pro$splitsubTAD = FALSE
for(i in 1:nrow(subsplit)){
  lev = as.numeric(paste(HE$tadlevel[which(HE$he_domain == subsplit$he_domain[i])]))
  HPC_Pro$splitsubTAD[which(HPC_Pro$chr == subsplit$V4[i] & HPC_Pro$start >= subsplit$V5[i] & HPC_Pro$end <= subsplit$V6[i] & HPC_Pro$tadlevel >= lev)] = TRUE
}

merg=read.table(file.path(DATA_DIR, "hiC/compareTAD/domain_merged.txt"))
merg$he_domain = paste(merg$V1,merg$V2,merg$V3,sep=".")
merg$hpc_domain = paste(merg$V4,merg$V5,merg$V6,sep=".")
HE$mergeTAD = FALSE
for(i in 1:nrow(merg)){
  HE$mergeTAD[which(HE$chr == merg$V1[i] & HE$start >= merg$V2[i] & HE$end <= merg$V3[i] & HE$tadlevel==0)] = TRUE
}
HPC_Pro$mergeTAD = HPC_Pro$hpc_domain %in% merg$hpc_domain

submerge=read.table(file.path(DATA_DIR, "hiC/compareTAD/subTAD_merged.txt"))
submerge$he_domain = paste(submerge$V1,submerge$V2,submerge$V3,sep=".")
submerge$hpc_domain = paste(submerge$V4,submerge$V5,submerge$V6,sep=".")
HE$mergesubTAD = FALSE
HPC_Pro$mergesubTAD = HPC_Pro$hpc_domain %in% submerge$hpc_domain
for(i in 1:nrow(submerge)){
  lev = HPC_Pro$tadlevel[which(HPC_Pro$hpc_domain == submerge$hpc_domain[i])]
  HE$mergesubTAD [which(HE$chr == submerge$V1[i] & HE$start >= submerge$V2[i] & HE$end <= submerge$V3[i] & HE$tadlevel >= lev)] = TRUE
}


unique(HE[,6:14])
unique(HPC_Pro[,6:14])

HE$type = apply(HE[,6:14],1,function(x) paste(colnames(HE)[which(x==TRUE)+5],collapse = ","))
HPC_Pro$type = apply(HPC_Pro[,6:14],1,function(x) paste(colnames(HPC_Pro)[which(x==TRUE)+5],collapse = ","))

HPC_Pro$splitsubTAD[which(HPC_Pro$type=="conservedsubTAD,splitsubTAD")]=FALSE
HPC_Pro$splitsubTAD[which(HPC_Pro$type=="subTADhigher,splitsubTAD")]=FALSE
HPC_Pro$splitsubTAD[which(HPC_Pro$type == "splitsubTAD,mergesubTAD")]=FALSE
HPC_Pro$type = apply(HPC_Pro[,6:14],1,function(x) paste(colnames(HPC_Pro)[which(x==TRUE)+5],collapse = ","))
HPC_Pro$type[which(HPC_Pro$type=="")] = "Undefined"
HE$type[which(HE$type=="")] = "Undefined"

he=GRanges(HE[which(HE$type=="Undefined"),1:5])
hpc=GRanges(HPC_Pro[which(HPC_Pro$type=="Undefined"),1:5])
a=findOverlaps(hpc,he)
int <- pintersect(hpc[queryHits(a)], he[subjectHits(a)])
int_width <- width(int)
hpc_width <- width(hpc[queryHits(a)])
he_width <- width(he[subjectHits(a)])
overlap_fraction <- pmin(int_width / hpc_width, int_width / he_width)
keep <- overlap_fraction >= 0.9
filtered_hits <- a[keep]
hpc_overlaps <- hpc[queryHits(filtered_hits)]
he_overlaps <- he[subjectHits(filtered_hits)]
filtered = cbind(data.frame(hpc_overlaps),data.frame(he_overlaps))
colnames(filtered)[13]="he_level"
filtered=filtered[,c('hpc_domain','tadlevel','he_domain','he_level')]
a=which(filtered$tadlevel==filtered$he_level & filtered$hpc_domain==filtered$he_domain)
HE$conservedsubTAD[which(HE$he_domain %in% filtered$he_domain[a])] = TRUE
HPC_Pro$conservedsubTAD[which(HPC_Pro$hpc_domain %in% filtered$hpc_domain[a])] = TRUE

HE$type = apply(HE[,6:14],1,function(x) paste(colnames(HE)[which(x==TRUE)+5],collapse = ","))
HPC_Pro$type = apply(HPC_Pro[,6:14],1,function(x) paste(colnames(HPC_Pro)[which(x==TRUE)+5],collapse = ","))
he=GRanges(HE[which(HE$type==""),1:5])
hpc=GRanges(HPC_Pro[which(HPC_Pro$type==""),1:5])
a=findOverlaps(hpc,he)
int <- pintersect(hpc[queryHits(a)], he[subjectHits(a)])
int_width <- width(int)
hpc_width <- width(hpc[queryHits(a)])
he_width <- width(he[subjectHits(a)])
overlap_fraction <- pmin(int_width / hpc_width, int_width / he_width)
keep <- overlap_fraction >= 0.9
filtered_hits <- a[keep]
hpc_overlaps <- hpc[queryHits(filtered_hits)]
he_overlaps <- he[subjectHits(filtered_hits)]
filtered = cbind(data.frame(hpc_overlaps),data.frame(he_overlaps))
colnames(filtered)[13]="he_level"
filtered=filtered[,c('hpc_domain','tadlevel','he_domain','he_level')]
a=which(filtered$tadlevel==filtered$he_level & filtered$tadlevel==0)
HE$conservedTAD_shifted = HE$he_domain %in% filtered$he_domain[a]
HPC_Pro$conservedTAD_shifted = HPC_Pro$hpc_domain %in% filtered$hpc_domain[a]
a=which(filtered$tadlevel==filtered$he_level & filtered$tadlevel!=0)
HE$conservedsubTAD_shifted = HE$he_domain %in% filtered$he_domain[a]
HPC_Pro$conservedsubTAD_shifted = HPC_Pro$hpc_domain %in% filtered$hpc_domain[a]

HE=HE[,-15]
HPC_Pro=HPC_Pro[,-15]
HE$type = apply(HE[,c(6:16)],1,function(x) paste(colnames(HE)[which(x==TRUE)+5],collapse = ","))
HPC_Pro$type = apply(HPC_Pro[,6:16],1,function(x) paste(colnames(HPC_Pro)[which(x==TRUE)+5],collapse = ","))

he=GRanges(HE[which(HE$type==""),1:5])
hpc=GRanges(HPC_Pro[which(HPC_Pro$type==""),1:5])
a=findOverlaps(hpc,he)
int <- pintersect(hpc[queryHits(a)], he[subjectHits(a)])
int_width <- width(int)
hpc_width <- width(hpc[queryHits(a)])
he_width <- width(he[subjectHits(a)])
overlap_fraction <- pmin(int_width / hpc_width, int_width / he_width)
keep <- overlap_fraction >= 0.9
filtered_hits <- a[keep]
hpc_overlaps <- hpc[queryHits(filtered_hits)]
he_overlaps <- he[subjectHits(filtered_hits)]
filtered = cbind(data.frame(hpc_overlaps),data.frame(he_overlaps))
colnames(filtered)[13]="he_level"
filtered=filtered[,c('hpc_domain','tadlevel','he_domain','he_level')]
a=which(filtered$he_level==0 & filtered$tadlevel>filtered$he_level)
HE$TADbecomesubTAD = HE$he_domain %in% filtered$he_domain[a]
HPC_Pro$TADbecomesubTAD = HPC_Pro$hpc_domain %in% filtered$hpc_domain[a]
filtered=filtered[-a,]
a=which(filtered$tadlevel==0)
HE$subTADbecomeTAD = HE$he_domain %in% filtered$he_domain[a]
HPC_Pro$subTADbecomeTAD = HPC_Pro$hpc_domain %in% filtered$hpc_domain[a]
filtered=filtered[-a,]
a=which(filtered$tadlevel>filtered$he_level)
HE$subTADhigher = HE$he_domain %in% filtered$he_domain[a]
HPC_Pro$subTADhigher = HPC_Pro$hpc_domain %in% filtered$hpc_domain[a]
filtered=filtered[-a,]
a=which(filtered$tadlevel<filtered$he_level)
HE$subTADlower = HE$he_domain %in% filtered$he_domain[a]
HPC_Pro$subTADlower = HPC_Pro$hpc_domain %in% filtered$hpc_domain[a]

HE=HE[,-17]
HPC_Pro=HPC_Pro[,-17]
HE$type = apply(HE[,c(6:18)],1,function(x) paste(colnames(HE)[which(x==TRUE)+5],collapse = ","))
HPC_Pro$type = apply(HPC_Pro[,6:18],1,function(x) paste(colnames(HPC_Pro)[which(x==TRUE)+5],collapse = ","))

he=GRanges(HE[which(HE$type==""),1:5])
hpc=GRanges(HPC_Pro[which(HPC_Pro$type==""),1:5])
a=findOverlaps(hpc,he)
hpc_overlaps <- hpc[queryHits(a)]
he_overlaps <- he[subjectHits(a)]
filtered = cbind(data.frame(hpc_overlaps),data.frame(he_overlaps))
colnames(filtered)[13]="he_level"
filtered=filtered[,c('hpc_domain','tadlevel','he_domain','he_level')]
HE$complex = HE$he_domain %in% filtered$he_domain
HPC_Pro$complex = HPC_Pro$hpc_domain %in% filtered$hpc_domain
HE$disappeared = HE$he_domain %notin% filtered$he_domain & HE$type==""
HPC_Pro$newlyform = HPC_Pro$hpc_domain %notin% filtered$hpc_domain & HPC_Pro$type==""
HE=HE[,-19]
HPC_Pro=HPC_Pro[,-19]
HE$type = apply(HE[,c(6:20)],1,function(x) paste(colnames(HE)[which(x==TRUE)+5],collapse = ","))
HPC_Pro$type = apply(HPC_Pro[,6:20],1,function(x) paste(colnames(HPC_Pro)[which(x==TRUE)+5],collapse = ","))



HPC_Pro$type = factor(HPC_Pro$type , levels = c("Undefined","splitTAD","mergeTAD","semiTAD","conservedTAD",
                                                "splitsubTAD","mergesubTAD","subTADlower","subTADhigher","conservedsubTAD"))

chang_summary <- rbind(HE[which(HE$disappeared & HE$tadlevel==0),c("chr","start","end","type")],
                       HPC_Pro[which(HPC_Pro$tadlevel==0),c("chr","start","end","type")]) %>%
  count(type) %>% # Count the number of observations in each 'trans' category
  mutate(
    percentage = n / sum(n) * 100, # Calculate the percentage
    # Create a label string with percentage, formatted to 1 decimal place
    label = paste0(round(percentage, 1), "%")
  ) %>%
  # Arrange in descending order of count so that the labels can be
  # positioned correctly on the stacked bar before polar transformation.
  # Then calculate the cumulative percentage to help with label positioning.
  arrange(desc(type)) %>% # Arrange by 'trans' if needed for specific stacking order
  mutate(
    csum = rev(cumsum(rev(percentage))), 
    text_position = percentage/2 + lead(csum, 1),
    text_position = if_else(is.na(text_position), percentage/2, text_position)
  )
g1=ggplot(chang_summary, aes(x = "", y = percentage, fill = fct_inorder(type))) + # Note the mapping: x="", y=percentage
  geom_bar(stat = "identity", width = 1, color = "white") + # Use stat="identity" as y is already the percentage
  coord_polar(theta = "y") +
  scale_fill_observable(name="TAD changes")+
  theme_bw() +
  theme(
    axis.text = element_blank(),
    panel.grid = element_blank(),
    panel.border = element_blank(),
    axis.ticks = element_blank(),
    axis.title = element_blank(),
    legend.title = element_text(hjust = 0)
  ) +
  # Add the text layer for percentages
  geom_text(aes(1.1,y = text_position, label = label), # Use the calculated text_position and label
            color = "black", size = 4) # Customize text color and size

chang_summary2 <- rbind(HE[which(HE$disappeared & HE$tadlevel!=0),c("chr","start","end","type")],
                        HPC_Pro[which(HPC_Pro$tadlevel!=0),c("chr","start","end","type")]) %>%
  count(type) %>% # Count the number of observations in each 'trans' category
  mutate(
    percentage = n / sum(n) * 100, # Calculate the percentage
    label = paste0(round(percentage, 1), "%")
  ) %>%
  arrange(desc(type)) %>% # Arrange by 'trans' if needed for specific stacking order
  mutate(
    csum = rev(cumsum(rev(percentage))), 
    text_position = percentage/2 + lead(csum, 1),
    text_position = if_else(is.na(text_position), percentage/2, text_position)
  )
g2=ggplot(chang_summary2, aes(x = "", y = percentage, fill = fct_inorder(type))) + # Note the mapping: x="", y=percentage
  geom_bar(stat = "identity", width = 1, color = "white") + # Use stat="identity" as y is already the percentage
  coord_polar(theta = "y") +
  scale_fill_observable(name="subTAD changes")+
  theme_bw() +
  theme(
    axis.text = element_blank(),
    panel.grid = element_blank(),
    panel.border = element_blank(),
    axis.ticks = element_blank(),
    axis.title = element_blank(),
    legend.title = element_text(hjust = 0)
  ) +
  # Add the text layer for percentages
  geom_text(aes(1.1,y = text_position, label = label), # Use the calculated text_position and label
            color = "black", size = 4) # Customize text color and size


pdf(file.path(DATA_DIR, "hiC/compareTAD/HPC_domain_changes.pdf"),height = 4,width = 12)
g1+g2
dev.off()

# semi

# split

# subsplit submerge

#merge

a = rbind(HE[which(HE$disappeared & HE$tadlevel==0),c("chr","start","end","type")],
          HPC_Pro[which(HPC_Pro$tadlevel==0),c("chr","start","end","type")])
rownames(a)=NULL
a=GenomicRanges::GRanges(a)
saveRDS(a,file=file.path(DATA_DIR, "hiC/compareTAD/TADcompare_diff.granges.rds"))

a = rbind(HE[which(HE$disappeared & HE$tadlevel!=0),c("chr","start","end","type")],
          HPC_Pro[which(HPC_Pro$tadlevel!=0),c("chr","start","end","type")])
rownames(a)=NULL
a=GenomicRanges::GRanges(a)
saveRDS(a,file=file.path(DATA_DIR, "hiC/compareTAD/subTADcompare_diff.granges.rds"))

consensus=readRDS("compareTAD/consensus_range.rds")
overlaps <- findOverlaps(consensus, a)
saveRDS(overlaps,file=file.path(DATA_DIR, "hiC/compareTAD/overlaps.rds"))

promoter_ranges=readRDS("compareTAD/promoter_ranges.rds")
overlaps <- findOverlaps(consensus, a)
saveRDS(overlaps,file=file.path(DATA_DIR, "hiC/compareTAD/sub.overlaps.rds"))



library(ggsignif)
library(ggpubr)
# Boundaries aligner analysis ####
sub2tad=read.table(file.path(DATA_DIR, "hiC/compareTAD/sub2tad.txt"))
sub2tad$type="subTAD to TAD boundaries"
tad2sub=read.table(file.path(DATA_DIR, "hiC/compareTAD/tad2sub.txt"))
tad2sub$type="TAD to subTAD boundaries"
distad=read.table(file.path(DATA_DIR, "hiC/compareTAD/bound_disappeared_tad.txt"))
distad$type="Disappeared TAD boundaries"
dissub=read.table(file.path(DATA_DIR, "hiC/compareTAD/bound_disappeared_sub.txt"))
dissub$type="Disappeared subTAD boundaries"
distad2=read.table(file.path(DATA_DIR, "hiC/compareTAD/bound_disappeared_tad.HPC.txt"))
distad2$type="New TAD boundaries"
dissub2=read.table(file.path(DATA_DIR, "hiC/compareTAD/bound_disappeared_sub.HPC.txt"))
dissub2$type="New subTAD boundaries"
colnames(distad2)=colnames(dissub2)=colnames(distad)=colnames(dissub)=c("V3","V4","type")
sub=read.table(file.path(DATA_DIR, "hiC/compareTAD/bound_conserved_sub_bounds.txt"))
sub$type="Conserved subTAD boundaries"
tad=read.table(file.path(DATA_DIR, "hiC/compareTAD/bound_conserved_tad_bounds.txt"))
tad$type="Conserved TAD boundaries"

hebound=unique(rbind(sub2tad[,c(1,2,5)],
              tad2sub[,c(1,2,5)],
              distad,
              dissub,
              sub[,c(1,2,5)],
              tad[,c(1,2,5)]))
hpcbound=unique(rbind(sub2tad[,c(3,4,5)],
                     tad2sub[,c(3,4,5)],
                     distad,
                     dissub,
                     distad2,
                     dissub2,
                     sub[,c(3:5)],
                     tad[,c(3:5)]))

chang_summary <- hpcbound %>%
  count(type) %>% # Count the number of observations in each 'trans' category
  mutate(
    percentage = n / sum(n) * 100, # Calculate the percentage
    # Create a label string with percentage, formatted to 1 decimal place
    label = paste0(round(percentage, 1), "%")
  ) %>%
  # Arrange in descending order of count so that the labels can be
  # positioned correctly on the stacked bar before polar transformation.
  # Then calculate the cumulative percentage to help with label positioning.
  arrange(desc(type)) %>% # Arrange by 'trans' if needed for specific stacking order
  mutate(
    csum = rev(cumsum(rev(percentage))), 
    text_position = percentage/2 + lead(csum, 1),
    text_position = if_else(is.na(text_position), percentage/2, text_position)
  )
chang_summary$type = factor(chang_summary$type,levels=chang_summary$type)
pdf(file.path(DATA_DIR, "hiC/compareTAD/HPC_boundaries_changes.pdf"),height = 4,width = 6)
ggplot(chang_summary, aes(x = "", y = percentage, fill = fct_inorder(type))) + # Note the mapping: x="", y=percentage
  geom_bar(stat = "identity", width = 1, color = "white") + # Use stat="identity" as y is already the percentage
  coord_polar(theta = "y") +
  theme_bw() +
  theme(
    axis.text = element_blank(),
    panel.grid = element_blank(),
    panel.border = element_blank(),
    axis.ticks = element_blank(),
    axis.title = element_blank(),
    legend.title = element_text(hjust = 0)
  ) +
  guides(fill = guide_legend(title = "TAD-hierarchical\nboundary changes")) +
  # Add the text layer for percentages
  geom_text(aes(1.1,y = text_position, label = label), # Use the calculated text_position and label
            color = "black", size = 4) # Customize text color and size
dev.off()

chang_summary <- hebound %>%
  count(type) %>% # Count the number of observations in each 'trans' category
  mutate(
    percentage = n / sum(n) * 100, # Calculate the percentage
    # Create a label string with percentage, formatted to 1 decimal place
    label = paste0(round(percentage, 1), "%")
  ) %>%
  # Arrange in descending order of count so that the labels can be
  # positioned correctly on the stacked bar before polar transformation.
  # Then calculate the cumulative percentage to help with label positioning.
  arrange(desc(type)) %>% # Arrange by 'trans' if needed for specific stacking order
  mutate(
    csum = rev(cumsum(rev(percentage))), 
    text_position = percentage/2 + lead(csum, 1),
    text_position = if_else(is.na(text_position), percentage/2, text_position)
  )
chang_summary$type = factor(chang_summary$type,levels=chang_summary$type)
pdf(file.path(DATA_DIR, "hiC/compareTAD/HE_boundaries_changes.pdf"),height = 5,width = 8)
ggplot(chang_summary, aes(x = "", y = percentage, fill = fct_inorder(type))) + # Note the mapping: x="", y=percentage
  geom_bar(stat = "identity", width = 1, color = "white") + # Use stat="identity" as y is already the percentage
  coord_polar(theta = "y") +
  theme_bw() +
  theme(
    axis.text = element_blank(),
    panel.grid = element_blank(),
    panel.border = element_blank(),
    axis.ticks = element_blank(),
    axis.title = element_blank(),
    legend.title = element_text(hjust = 0)
  ) +
  guides(fill = guide_legend(title = "TAD-hierarchical\nboundary changes")) +
  # Add the text layer for percentages
  geom_text(aes(1.1,y = text_position, label = label), # Use the calculated text_position and label
            color = "black", size = 4) # Customize text color and size
dev.off()



# AB compartment #######
setwd(file.path(DATA_DIR, "hiC/ABcompartments"))
HE1=read.table("HE_D9_plus7_rep1.eigen.40K.cis.vecs.tsv",sep = "\t",header=T)
HE2=read.table("HE_D9_plus7_rep2.eigen.40K.cis.vecs.tsv",sep = "\t",header=T)
HPC1=read.table("HPC_Pro_D9_plus8_rep1.eigen.40K.cis.vecs.tsv",sep = "\t",header = T)
HPC2=read.table("HPC_Pro_D9_plus8_rep2.eigen.40K.cis.vecs.tsv",sep = "\t",header = T)
all(HE1$end==HPC2$end)
ma=cbind(HE_D9_plus7_rep1 = HE1$E1>0,
         HE_D9_plus7_rep2 = HE2$E1>0,
         HPC_Pro_D9_plus8_rep1 = HPC1$E1>0,
         HPC_Pro_D9_plus8_rep2 = HPC2$E1>0)
rownames(ma)=seq(1,nrow(ma))
ma=ma[-which(is.na(rowSums(ma))),]
ma=as.data.frame(ma)

jaccard_similarity <- function(A, B) { 
  intersection = length(which(A==B) )
  union = length(A) + length(B) - intersection 
  return (intersection/union) 
}

pairwise_jaccard <- function(matrix) {
  n <- ncol(matrix)  # Number of columns
  jaccard_matrix <- matrix(0, nrow = n, ncol = n)  # Initialize a matrix to store Jaccard coefficients
  
  # Calculate pairwise Jaccard coefficients
  for (i in 1:(n - 1)) {
    for (j in (i + 1):n) {
      jaccard_matrix[i, j] <- jaccard_similarity(matrix[, i], matrix[, j])
      jaccard_matrix[j, i] <- jaccard_matrix[i, j]  # Since the Jaccard matrix is symmetric
    }
  }
  diag(jaccard_matrix) <- 1
  # Return the pairwise Jaccard coefficient matrix
  return(jaccard_matrix)
}
jac=pairwise_jaccard(ma)
colnames(jac)=rownames(jac)=colnames(ma)


library(corrplot)
pdf(file.path(DATA_DIR, "hiC/ABcompartments/jaccard_pairwise.pdf"))
corrplot(jac, type = "upper", is.corr = F,col.lim = c(0.75,1),
        col= COL1('OrRd', 200),addCoef.col = "black",
         tl.col = "black", tl.srt = 45)
dev.off()
rm(HE1,HE2,HPC1,HPC2)


HE = read.table(file.path(DATA_DIR, "hiC/ABcompartments/HE_D9_plus7.eigen.40K.cis.vecs.tsv"),sep = "\t",header=T)
HPC = read.table(file.path(DATA_DIR, "hiC/ABcompartments/HPC_Pro_D9_plus8.eigen.40K.cis.vecs.tsv"),sep = "\t",header=T)
HE$HE = ifelse(HE$E1>0,"A","B")
HPC$HPC = ifelse(HPC$E1>0,"A","B")
chang = cbind(HE[,c(1:3,8)],HPC=HPC$HPC)
chang=chang[-which(is.na(chang$HE)),]
chang=chang[-which(is.na(chang$HPC)),]
chang$trans = paste0(chang$HE," -> ",chang$HPC)
chang$trans = factor(chang$trans,levels = c("A -> B", "B -> A", "A -> A", "B -> B"))
library(dplyr)
chang_summary <- chang %>%
  count(trans) %>% # Count the number of observations in each 'trans' category
  mutate(
    percentage = n / sum(n) * 100, # Calculate the percentage
    # Create a label string with percentage, formatted to 1 decimal place
    label = paste0(round(percentage, 1), "%")
  ) %>%
  # Arrange in descending order of count so that the labels can be
  # positioned correctly on the stacked bar before polar transformation.
  # Then calculate the cumulative percentage to help with label positioning.
  arrange(desc(trans)) %>% # Arrange by 'trans' if needed for specific stacking order
  mutate(
    cum_percent = cumsum(percentage),
    # Calculate the position for the text labels.
    # In a stacked bar (before polar), the label position is typically
    # in the middle of each segment.
    # The y-position for text in a stacked bar is the cumulative sum minus half of the current segment's size.
    text_position = cum_percent - (percentage / 2)
  )


pdf(file.path(DATA_DIR, "hiC/ABcompartments/compartment_witch.pdf"),height = 3,width = 4)
ggplot(chang_summary, aes(x = "", y = percentage, fill = trans)) + # Note the mapping: x="", y=percentage
  geom_bar(stat = "identity", width = 1, color = "white") + # Use stat="identity" as y is already the percentage
  coord_polar(theta = "y",start = 90) +
  theme_bw() +
  theme(
    axis.text = element_blank(),
    panel.grid = element_blank(),
    panel.border = element_blank(),
    axis.ticks = element_blank(),
    axis.title = element_blank(),
    legend.title = element_text(hjust = 0.5)
  ) +
  guides(fill = guide_legend(title = "Compartment\nswitch")) +
  # Add the text layer for percentages
  geom_text(aes(y = text_position, label = label), # Use the calculated text_position and label
            color = "black", size = 4) # Customize text color and size
dev.off()

write.table(chang[which(chang$trans=="A -> B"),1:3],file=file.path(DATA_DIR, "hiC/AtoB.bed"),sep="\t",col.names = F,row.names = F,quote=F)
write.table(chang[which(chang$trans=="B -> A"),1:3],file=file.path(DATA_DIR, "hiC/BtoA.bed"),sep="\t",col.names = F,row.names = F,quote=F)


library(ChIPseeker)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
 
AtoB=annotatePeak(peak = GRanges(chang[which(chang$trans=="A -> B"),1:3]),TxDb = TxDb.Hsapiens.UCSC.hg38.knownGene,tssRegion = c(-2000, 2000))
annotations_edb <- AnnotationDbi::select(EnsDb.Hsapiens.v86,
                                         keys = unique(AtoB$geneId),
                                         columns = c("GENENAME"),
                                         keytype = "ENTREZID")
AtoB = merge(AtoB,annotations_edb,by.x="geneId",by.y="ENTREZID")
write.table(unique(AtoB[grep("Promoter",AtoB$annotation),c("GENENAME","seqnames","start","end")]),file=file.path(DATA_DIR, "hiC/AtoB.promoter.anno.txt"),col.names = T,row.names = F,sep="\t",quote=F)

BtoA=annotatePeak(peak = GRanges(chang[which(chang$trans=="B -> A"),1:3]),TxDb = TxDb.Hsapiens.UCSC.hg38.knownGene,tssRegion = c(-2000, 2000))
annotations_edb <- AnnotationDbi::select(EnsDb.Hsapiens.v86,
                                         keys = unique(BtoA@anno$geneId),
                                         columns = c("GENENAME"),
                                         keytype = "ENTREZID")
BtoA=as.data.frame(BtoA)
BtoA = merge(BtoA,annotations_edb,by.x="geneId",by.y="ENTREZID",all.x=T)
write.table(unique(BtoA[grep("Promoter",BtoA$annotation),c("GENENAME","seqnames","start","end")]),file=file.path(DATA_DIR, "hiC/BtoA.promoter.anno.txt"),col.names = T,row.names = F,sep="\t",quote=F)

