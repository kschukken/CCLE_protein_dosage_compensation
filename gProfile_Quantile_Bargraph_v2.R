###### g:profile Gene enrichemnt analysis and plots ####
## 210303
## g:profile analysis 
## Categorical (Buffer/scaling/antiscaling, down /up-regulated) Protein expression 
## upon chromosome gain or loss. 

## Klaske Schukken

library(ggplot2)
library(tidyverse) 
library(readxl)
library(ggpubr)
library("cowplot")
library(plyr)
library('gprofiler2')


# Notes: RNA/Protein expression difference upon chrm gain/loss (difference, not corr!)
# Difference in expression can be + or -. 
# with filtered cell lines only (only cells with RNA & Protein data)

# To do: 
# per category I am looking at, run samples through g:profiler 
# either as categorical groups or as ordered list. 
# Sub-Step 1: Get g:profiler data for RNA/Protein Chrm Gain/Loss 
# Sub-Step 2: filter for top 5-10 GO terms 
# SubStep 3: bargraph -log2() of p-value for top terms for all conditions. 

###### Step 1: Get data #####

setwd() #set location 
#make in Protein_RNA_expression.PerCell_v2.R
# also available from supplementary data ** sheet 2.
CN.Diff.xRNA.yProt.ThreeGroups<- read.csv("RNA.Protein_Mono.Di.Tri_Difference_Pvalue_min10points_3category.csv")
CN.Diff.xRNA.yProt.ThreeGroups<- CN.Diff.xRNA.yProt.ThreeGroups[,-c(1)]


### These datasets are from g:profiler, the online interface, where I submitted RNA and proteins 
# as ordered lists (highest to lowest, and lowest to highest) based on their correlation with cellular aneuploidy
# then I uploaded these lists to this script. 

## available from Supplementary Data **
RNA.negSig<-read.delim2("gProfiler_RNA.negCor.AnScore.expression.BHSig.csv", 
                        dec=".", header = TRUE, sep=",") 
RNA.posSig<-read.delim2("gProfiler_RNA.posCor.AnScore.expression.BHSig.csv", 
                        dec=".", header = TRUE, sep=",") 

Protein.negSig<-read.delim2("gProfiler_Protein.NegCor.AneuploidyScore.expression.BHSig.csv", 
                            dec=".", header = TRUE, sep=",") 
Protein.posSig<-read.delim2("gProfiler_Protein.PosCorr.AneuploidyScore.BHSig.csv", 
                            dec=".", header = TRUE, sep=",") 



###### Figure 4: gene expression correlates with aneuploidy score #####
## Figure 4: gene expression correlates with aneuploidy score
## Look at RNA and Proteins positively or negatively correlated with aneuploidy score

#make bargraph of key terms/hits RNA
# RNA positively and negatively correlated with aneuploidy
R.KeyNegSig<-RNA.negSig[order(RNA.negSig$adjusted_p_value),]
R.KeyNegSig1<- subset(R.KeyNegSig, source %in%  c("CORUM" ,"GO:BP", "GO:MF", "WP", "KEGG", "REAC"))
R.KeyNegSig2<- R.KeyNegSig1[c(1,9,10,12,38, 13,25,26,32,37),]#ribosome (2,4,7,), translation, protease
R.KeyNegSig2$term_name<- factor(R.KeyNegSig2$term_name, levels= R.KeyNegSig2$term_name)
R.KeyNegSig2$Termtype<-c("Ribosome","Ribosome","Ribosome","Ribosome","Ribosome",
                       "RNA processing","RNA processing","RNA processing", "RNA processing","RNA processing")

R.KeyPosSig<-RNA.posSig[order(RNA.posSig$adjusted_p_value),]
R.KeyPosSig1<- subset(R.KeyPosSig, source %in% c("CORUM" ,"GO:BP", "GO:MF", "WP", "KEGG", "REAC"))
R.KeyPosSig2<-R.KeyPosSig1[c(2,8,19,91,174,201,  412,413),]#protein complexes, DNA damage response, heat shock
R.KeyPosSig2$term_name<- factor(R.KeyPosSig2$term_name, levels= R.KeyPosSig2$term_name)
R.KeyPosSig2$Termtype<-c("Adhesion", "Adhesion", "Adhesion", "Adhesion", "Adhesion", "Adhesion", 
                         "Heat acclimation", "Heat acclimation")

#subset(R.KeyPosSig1, grepl( "heat", R.KeyPosSig1$term_name, fixed = TRUE))
# look for heat shock terms.

plot.KeyTermsRNA<- ggplot(R.KeyPosSig2, 
                       aes(x=term_name, y=-log2(adjusted_p_value), fill=Termtype))+
  geom_bar(stat="identity")+
  ylab("-log2 (P-Value)") +
  xlab("Key RNA terms")+
  ggtitle ("Positively correlated with aneuploidy")+
  theme_classic()+
  ylim(0,140)+
  coord_flip()
plot.KeyTermsRNA#7x5
# plot.RNA.Pos.AneuScoreCorr.KeyTerms_3


#make bargraph of key terms/hits PROTEIN positively and negatively regulated with aneuploidy score
KeyNegSig<-Protein.negSig[order(Protein.negSig$adjusted_p_value),]
KeyNegSig2<-KeyNegSig[c(6,8,9,15, 2,24,27),]#ribosome, translation, protease
KeyNegSig2$term_name<- factor(KeyNegSig2$term_name, levels= KeyNegSig2$term_name)
KeyNegSig2$Termtype<-c("Ribosome","Ribosome","Ribosome","Ribosome",
                       "RNA processing", "RNA processing", "RNA processing")

KeyPosSig<-Protein.posSig[order(Protein.posSig$adjusted_p_value),]
#KeyPosSig2<-KeyPosSig[c(1,6,9,12,13,15,16,21,23,30, 14,18,20,26,28,53, 45,49),]#protein complexes, DNA damage response, heat shock
#KeyPosSig2<-KeyPosSig[c(1,20,26, 45,49),]#protein complexes, DNA damage response, heat shock. (1,6,9,10,14,18,20,26,28,53, 45,49)
KeyPosSig2<-KeyPosSig[c(22,29,50, 45,49),]#membranes, heat shock.
KeyPosSig2$Termtype<-c("Membrane","Membrane", "Membrane",
                       "Heat shock response", "Heat shock response")
KeyPosSig2$term_name<- factor(KeyPosSig2$term_name, levels= KeyPosSig2$term_name)


plot.KeyTerms<- ggplot(KeyPosSig2, 
                          aes(x=term_name, y=-log2(adjusted_p_value), fill=Termtype))+
  geom_bar(stat="identity")+
  ylab("-log2 (P-Value)") +
  xlab("Biological terms")+
  ggtitle ("Positively correlated with aneuploidy")+
  theme_classic()+
  coord_flip()+
  ylim(0,13)
plot.KeyTerms
# 6x4
# plot.PosCorr.AneuScore.terms



##
##origin of terms pos and negatively correlated with aneuploidy score
# this did not end up in the paper, but interesting
Protein.AnScore.Source<- rbind.fill.matrix(t(table(Protein.posSig$source)), t(table(Protein.negSig$source)))
row.names(Protein.AnScore.Source)<- c("Positive", "Negative") 
Protein.AnScore.Source2<-t(Protein.AnScore.Source)
GOMF<- c(0,0)
MIRNA<-c(0,0)
TF<-c(0,0)
Protein.AnScore.Source3 <- rbind( Protein.AnScore.Source2[1:3,], GOMF, Protein.AnScore.Source2[ 4:6,], MIRNA, Protein.AnScore.Source2[7,], TF ,Protein.AnScore.Source2[8,] )
row.names(Protein.AnScore.Source3)<- c("CORUM", "GO:BP", "GO:CC", "GO:MF", "HP", "HPA", "KEGG", "MIRNA", "REAC", "TF", "WP") 
Protein.AnScore.Source4<-Protein.AnScore.Source3[c(1,9),]
Protein.AnScore.Source_melt <- melt(Protein.AnScore.Source2, id = colnames)

RNA.AnScore.Source<- rbind(table(RNA.posSig$source), table(RNA.negSig$source))
row.names(RNA.AnScore.Source)<- c("Positive", "Negative") #Weird error: 
RNA.AnScore.Source2<-t(RNA.AnScore.Source)
RNA.AnScore.Source3<-RNA.AnScore.Source2[c(1,9),]
RNA.AnScore.Source_melt <- melt(RNA.AnScore.Source3, id = colnames)


plot.CVterms<- ggplot(Protein.AnScore.Source_melt, 
                      aes(x=Var1, y= value, fill=Var2))+
  geom_col(position = "dodge")+
  ylab("Count") +
  xlab("Enriched terms types")+
  scale_fill_discrete(name = "Correlated with \naneuploidy score")+
  ggtitle ("Origin of protein terms correlated with aneuploidy")+
  theme_classic()

plot.CVterms
# Lots of CORUM terms (protein complexes) pos correlated with aneuploidy
# lots of REAC terms (reactosome genes, signalling pathways) negatively correlated with aneuploidy. 



###### PROTEIN AntiScale/buffer/scale key terms ####
# Figure C/D: Buffering/scaling terms
## Three categories based on difference-- PROTEIN
# categorized by scaling, anti-scaling and buffering
# data generated in Protein_RNA_Corr_min10.Filtered.R

### Run through g:profiler
## !! Run g:profiler below, or upload results from supplementary data **

### Find proteins buffered upon gain: 
Prot.Gain.Buffering<- gost(
  subset(CN.Diff.xRNA.yProt.ThreeGroups, Three.Protein.Gain=="Buffering")$RNA_Name,
  organism = "hsapiens",
  ordered_query = FALSE,
  multi_query = FALSE,
  significant = TRUE,
  exclude_iea = FALSE,
  measure_underrepresentation = FALSE,
  evcodes = FALSE,
  user_threshold = 0.05,
  correction_method = c("g_SCS"),
  domain_scope = "custom",
  custom_bg = as.character(CN.Diff.xRNA.yProt.ThreeGroups$RNA_Name),
  numeric_ns = "",
  sources = NULL,
  as_short_link = FALSE
)
Prot.Gain.Buffering.result<-Prot.Gain.Buffering$result

## Proteins buffered upon loss:
Prot.Loss.Buffering<- gost(
  subset(CN.Diff.xRNA.yProt.ThreeGroups, Three.Protein.Loss=="Buffering")$RNA_Name,
  organism = "hsapiens",
  ordered_query = FALSE,
  multi_query = FALSE,
  significant = TRUE,
  exclude_iea = FALSE,
  measure_underrepresentation = FALSE,
  evcodes = FALSE,
  user_threshold = 0.05,
  correction_method = c("g_SCS"),
  domain_scope = "custom",
  custom_bg = as.character(CN.Diff.xRNA.yProt.ThreeGroups$RNA_Name),
  numeric_ns = "",
  sources = NULL,
  as_short_link = FALSE
)
Prot.Loss.Buffering.result<-Prot.Loss.Buffering$result

### Find proteins scaling upon gain: 
Prot.Gain.Scale<- gost(
  subset(CN.Diff.xRNA.yProt.ThreeGroups, Three.Protein.Gain=="Scaling")$RNA_Name,
  organism = "hsapiens",
  ordered_query = FALSE,
  multi_query = FALSE,
  significant = TRUE,
  exclude_iea = FALSE,
  measure_underrepresentation = FALSE,
  evcodes = FALSE,
  user_threshold = 0.05,
  correction_method = c("g_SCS"),
  domain_scope = "custom",
  custom_bg = as.character(CN.Diff.xRNA.yProt.ThreeGroups$RNA_Name),
  numeric_ns = "",
  sources = NULL,
  as_short_link = FALSE
)
Prot.Gain.Scale.result<-Prot.Gain.Scale$result

## Proteins scaling upon loss:
Prot.Loss.Scale<- gost(
  subset(CN.Diff.xRNA.yProt.ThreeGroups, Three.Protein.Loss=="Scaling")$RNA_Name,
  organism = "hsapiens",
  ordered_query = FALSE,
  multi_query = FALSE,
  significant = TRUE,
  exclude_iea = FALSE,
  measure_underrepresentation = FALSE,
  evcodes = FALSE,
  user_threshold = 0.05,
  correction_method = c("g_SCS"),
  domain_scope = "custom",
  custom_bg = as.character(CN.Diff.xRNA.yProt.ThreeGroups$RNA_Name),
  numeric_ns = "",
  sources = NULL,
  as_short_link = FALSE
)
Prot.Loss.Scale.result<-Prot.Loss.Scale$result


### Find proteins anti-scaling upon gain: 
Prot.Gain.AntiScale<- gost(
  subset(CN.Diff.xRNA.yProt.ThreeGroups, Three.Protein.Gain=="Anti-Scaling")$RNA_Name,
  organism = "hsapiens",
  ordered_query = FALSE,
  multi_query = FALSE,
  significant = TRUE,
  exclude_iea = FALSE,
  measure_underrepresentation = FALSE,
  evcodes = FALSE,
  user_threshold = 0.05,
  correction_method = c("g_SCS"),
  domain_scope = "custom",
  custom_bg = as.character(CN.Diff.xRNA.yProt.ThreeGroups$RNA_Name),
  numeric_ns = "",
  sources = NULL,
  as_short_link = FALSE
)
Prot.Gain.AntiScale.result<-Prot.Gain.AntiScale$result

## Proteins anti-scaling upon loss:
Prot.Loss.AntiScale<- gost(
  subset(CN.Diff.xRNA.yProt.ThreeGroups, Three.Protein.Loss=="Anti-Scaling")$RNA_Name,
  organism = "hsapiens",
  ordered_query = FALSE,
  multi_query = FALSE,
  significant = TRUE,
  exclude_iea = FALSE,
  measure_underrepresentation = FALSE,
  evcodes = FALSE,
  user_threshold = 0.05,
  correction_method = c("g_SCS"),
  domain_scope = "custom",
  custom_bg = as.character(CN.Diff.xRNA.yProt.ThreeGroups$RNA_Name),
  numeric_ns = "",
  sources = NULL,
  as_short_link = FALSE
)
Prot.Loss.AntiScale.result<-Prot.Loss.AntiScale$result

##Buffering terms
P.3Cat.Gain.Buff<-Prot.Gain.Buffering.result[order(Prot.Gain.Buffering.result$p_value),]
P.3Cat.Gain.Buff1<- subset(P.3Cat.Gain.Buff, source %in% c("CORUM" ,"GO:BP", "GO:MF", "WP", "KEGG", "REAC"))
P.3Cat.Gain.Buff2<-P.3Cat.Gain.Buff1[c(1, 9,13,60, 28,35,72),]#168, CORUM root, Ribosomal 32, cell cycle, splisosome
P.3Cat.Gain.Buff2$term_name<- factor(P.3Cat.Gain.Buff2$term_name, levels= P.3Cat.Gain.Buff2$term_name)
P.3Cat.Gain.Buff2$Termtype<-c("CORUM root", "Ribosomal", "Ribosomal", "Ribosomal", "Cell cycle", "Cell cycle", "Cell cycle")

P.3Cat.Loss.Buff<-Prot.Loss.Buffering.result[order(Prot.Loss.Buffering.result$p_value),]
P.3Cat.Loss.Buff1<- subset(P.3Cat.Loss.Buff, source %in% c("CORUM" ,"GO:BP", "GO:MF", "WP", "KEGG", "REAC"))
P.3Cat.Loss.Buff2<-P.3Cat.Loss.Buff1[c(1, 7,9,75, 91,105,113),]#180, CORUM root,  Ribosomal 24, cell cycle, splicing,
P.3Cat.Loss.Buff2$term_name<- factor(P.3Cat.Loss.Buff2$term_name, levels= P.3Cat.Loss.Buff2$term_name)
P.3Cat.Loss.Buff2$Termtype<-c("CORUM root", "Ribosomal", "Ribosomal", "Ribosomal", "Cell cycle", "Cell cycle", "Cell cycle")

## Scaling terms
P.3Cat.Gain.Scale<-Prot.Gain.Scale.result[order(Prot.Gain.Scale.result$p_value),]
#P.3Cat.Gain.Scale1<- subset(P.3Cat.Gain.Scale, source %in% c("CORUM" ,"GO:BP", "GO:MF", "WP", "KEGG", "REAC"))
P.3Cat.Gain.Scale2<-P.3Cat.Gain.Scale[c(2,3,4,5,14),]# 
P.3Cat.Gain.Scale2$term_name<- factor(P.3Cat.Gain.Scale2$term_name, levels= P.3Cat.Gain.Scale2$term_name)
P.3Cat.Gain.Scale2$Termtype<-c("Membrane", "Membrane", "Membrane", "Membrane", "Authophagy")

P.3Cat.Loss.Scale<-Prot.Loss.Scale.result[order(Prot.Loss.Scale.result$p_value),]
P.3Cat.Loss.Scale1<- subset(P.3Cat.Loss.Scale, source %in% c("CORUM" ,"GO:BP", "GO:MF", "WP", "KEGG", "REAC"))
P.3Cat.Loss.Scale2<-P.3Cat.Loss.Scale1[c(1,3,6,10,13),]#36, metabolism, 
P.3Cat.Loss.Scale2$term_name<- factor(P.3Cat.Loss.Scale2$term_name, levels= P.3Cat.Loss.Scale2$term_name)
P.3Cat.Loss.Scale2$Termtype<-c("Metabolism","Metabolism","Metabolism","Metabolism","Metabolism")

## anti-scaling terms
P.3Cat.Gain.AntiScale<-Prot.Gain.AntiScale.result[order(Prot.Gain.AntiScale.result$p_value),]
P.3Cat.Gain.AntiScale1<- subset(P.3Cat.Gain.AntiScale, source %in% c("CORUM" ,"GO:BP", "GO:MF", "WP", "KEGG", "REAC"))
P.3Cat.Gain.AntiScale2<-P.3Cat.Gain.AntiScale1[c(2,4,26,35, 3,5,16),]#58 total,  extracellular(adhesion, matrix, cell movement), ribosome, 
P.3Cat.Gain.AntiScale2$term_name<- factor(P.3Cat.Gain.AntiScale2$term_name, levels= P.3Cat.Gain.AntiScale2$term_name)
P.3Cat.Gain.AntiScale2$Termtype<-c("Extracellular","Extracellular","Extracellular","Extracellular", "Ribosomal","Ribosomal","Ribosomal")

P.3Cat.Loss.AntiScale<-Prot.Loss.AntiScale.result[order(Prot.Loss.AntiScale.result$p_value),]
P.3Cat.Loss.AntiScale1<- subset(P.3Cat.Loss.AntiScale, source %in% c("CORUM" ,"GO:BP", "GO:MF", "WP", "KEGG", "REAC"))
P.3Cat.Loss.AntiScale2<-P.3Cat.Loss.AntiScale1[c(2,4,11,26),]#135 total, extracellular (adhesion, matrix, cell movement), 
P.3Cat.Loss.AntiScale2$term_name<- factor(P.3Cat.Loss.AntiScale2$term_name, levels= P.3Cat.Loss.AntiScale2$term_name)
P.3Cat.Loss.AntiScale2$Termtype<-c("Extracellular","Extracellular","Extracellular","Extracellular")


ggplot(P.3Cat.Gain.Scale2, 
                          aes(x=term_name, y=-log2(p_value), fill=Termtype))+
  geom_bar(stat="identity")+
  ylab("-log2 (P-Value)") +
  xlab("Key RNA terms")+
  ggtitle ("Biological terms enriched: Protein gain Scale")+
  theme_classic()+
  ylim(0,160)+
  coord_flip()
#7x5
# plot.gprofiler.3Term.Protein.Gain.Scale_v3

# Write those terms as csv files
#setwd()
#write.csv(Prot.Loss.AntiScale.result[,1:13], 
#          file =paste("gprofiler.Protein.Loss.3cat.AntiScale.csv",
#                      sep=','), row.names = TRUE)

#P.3Cat.Gain.Scale<- read.csv("gprofiler.Protein.Gain.3cat.Scale.csv")
#P.3Cat.Loss.Scale<- read.csv("gprofiler.Protein.Loss.3cat.Scale.csv")



###### RNA     AntiScale/buffer/scale key terms ####
## Figure C/D: Buffering/scaling terms
## Three categories based on difference-- RNA
# categorized by scaling, anti-scaling and buffering
# data generated in Protein_RNA_Corr_min10.Filtered.R

### Run through g:profiler
## !! Run g:profiler below, or upload results from supplementary data **

### Find RNAs buffered upon gain: 
RNA.Gain.Buffering<- gost(
  subset(CN.Diff.xRNA.yProt.ThreeGroups, Three.RNA.Gain=="Buffering")$RNA_Name,
  organism = "hsapiens",
  ordered_query = FALSE,
  multi_query = FALSE,
  significant = TRUE,
  exclude_iea = FALSE,
  measure_underrepresentation = FALSE,
  evcodes = FALSE,
  user_threshold = 0.05,
  correction_method = c("g_SCS"),
  domain_scope = "custom",
  custom_bg = as.character(CN.Diff.xRNA.yProt.ThreeGroups$RNA_Name),
  numeric_ns = "",
  sources = NULL,
  as_short_link = FALSE
)
RNA.Gain.Buffering.result<-RNA.Gain.Buffering$result

## RNAs buffered upon loss:
RNA.Loss.Buffering<- gost(
  subset(CN.Diff.xRNA.yProt.ThreeGroups, Three.RNA.Loss=="Buffering")$RNA_Name,
  organism = "hsapiens",
  ordered_query = FALSE,
  multi_query = FALSE,
  significant = TRUE,
  exclude_iea = FALSE,
  measure_underrepresentation = FALSE,
  evcodes = FALSE,
  user_threshold = 0.05,
  correction_method = c("g_SCS"),
  domain_scope = "custom",
  custom_bg = as.character(CN.Diff.xRNA.yProt.ThreeGroups$RNA_Name),
  numeric_ns = "",
  sources = NULL,
  as_short_link = FALSE
)
RNA.Loss.Buffering.result<-RNA.Loss.Buffering$result

### Find RNAs scaling upon gain: 
RNA.Gain.Scale<- gost(
  subset(CN.Diff.xRNA.yProt.ThreeGroups, Three.RNA.Gain=="Scaling")$RNA_Name,
  organism = "hsapiens",
  ordered_query = FALSE,
  multi_query = FALSE,
  significant = TRUE,
  exclude_iea = FALSE,
  measure_underrepresentation = FALSE,
  evcodes = FALSE,
  user_threshold = 0.05,
  correction_method = c("g_SCS"),
  domain_scope = "custom",
  custom_bg = as.character(CN.Diff.xRNA.yProt.ThreeGroups$RNA_Name),
  numeric_ns = "",
  sources = NULL,
  as_short_link = FALSE
)
RNA.Gain.Scale.result<-RNA.Gain.Scale$result

## RNAs scaling upon loss:
RNA.Loss.Scale<- gost(
  subset(CN.Diff.xRNA.yProt.ThreeGroups, Three.RNA.Loss=="Scaling")$RNA_Name,
  organism = "hsapiens",
  ordered_query = FALSE,
  multi_query = FALSE,
  significant = TRUE,
  exclude_iea = FALSE,
  measure_underrepresentation = FALSE,
  evcodes = FALSE,
  user_threshold = 0.05,
  correction_method = c("g_SCS"),
  domain_scope = "custom",
  custom_bg = as.character(CN.Diff.xRNA.yProt.ThreeGroups$RNA_Name),
  numeric_ns = "",
  sources = NULL,
  as_short_link = FALSE
)
RNA.Loss.Scale.result<-RNA.Loss.Scale$result


### Find RNAs anti-scaling upon gain: 
RNA.Gain.AntiScale<- gost(
  subset(CN.Diff.xRNA.yProt.ThreeGroups, Three.RNA.Gain=="Anti-Scaling")$RNA_Name,
  organism = "hsapiens",
  ordered_query = FALSE,
  multi_query = FALSE,
  significant = TRUE,
  exclude_iea = FALSE,
  measure_underrepresentation = FALSE,
  evcodes = FALSE,
  user_threshold = 0.05,
  correction_method = c("g_SCS"),
  domain_scope = "custom",
  custom_bg = as.character(CN.Diff.xRNA.yProt.ThreeGroups$RNA_Name),
  numeric_ns = "",
  sources = NULL,
  as_short_link = FALSE
)
RNA.Gain.AntiScale.result<-RNA.Gain.AntiScale$result

## RNAs anti-scaling upon loss:
RNA.Loss.AntiScale<- gost(
  subset(CN.Diff.xRNA.yProt.ThreeGroups, Three.RNA.Loss=="Anti-Scaling")$RNA_Name,
  organism = "hsapiens",
  ordered_query = FALSE,
  multi_query = FALSE,
  significant = TRUE,
  exclude_iea = FALSE,
  measure_underrepresentation = FALSE,
  evcodes = FALSE,
  user_threshold = 0.05,
  correction_method = c("g_SCS"),
  domain_scope = "custom",
  custom_bg = as.character(CN.Diff.xRNA.yProt.ThreeGroups$RNA_Name),
  numeric_ns = "",
  sources = NULL,
  as_short_link = FALSE
)
RNA.Loss.AntiScale.result<-RNA.Loss.AntiScale$result

##Buffering terms RNA
R.3cat.Gain.Buff<-RNA.Gain.Buffering.result[order(RNA.Gain.Buffering.result$p_value),]
R.3cat.Gain.Buff1<- subset(R.3cat.Gain.Buff, source %in% c("CORUM" ,"GO:BP", "GO:MF", "WP", "KEGG", "REAC"))
R.3cat.Gain.Buff2<-R.3cat.Gain.Buff1[c(1,2),]#168, CORUM root, Ribosomal 32, cell cycle, splisosome
R.3cat.Gain.Buff2$term_name<- factor(R.3cat.Gain.Buff2$term_name, levels= R.3cat.Gain.Buff2$term_name)
R.3cat.Gain.Buff2$Termtype<-c("mRNA processing", "mRNA processing")

R.3cat.Loss.Buff<-RNA.Loss.Buffering.result[order(RNA.Loss.Buffering.result$p_value),]
R.3cat.Loss.Buff1<- subset(R.3cat.Loss.Buff, source %in% c("CORUM" ,"GO:BP", "GO:MF", "WP", "KEGG", "REAC"))
R.3cat.Loss.Buff2<-R.3cat.Loss.Buff1[c(1),]#180, CORUM root,  Ribosomal 24, cell cycle, splicing,
R.3cat.Loss.Buff2$term_name<- factor(R.3cat.Loss.Buff2$term_name, levels= R.3cat.Loss.Buff2$term_name)
R.3cat.Loss.Buff2$Termtype<-c("Extracellular")

## Scaling terms RNA
R.3cat.Gain.Scale<-RNA.Gain.Scale.result[order(RNA.Gain.Scale.result$p_value),]
R.3cat.Gain.Scale1<- subset(R.3cat.Gain.Scale, source %in% c("CORUM" ,"GO:BP", "GO:MF", "WP", "KEGG", "REAC"))
R.3cat.Gain.Scale2<-R.3cat.Gain.Scale1[c(1, 5),]#4 total. signaling pathways?
R.3cat.Gain.Scale2$term_name<- factor(R.3cat.Gain.Scale2$term_name, levels= R.3cat.Gain.Scale2$term_name)
R.3cat.Gain.Scale2$Termtype<-c("CORUM root", "Senescence")

R.3cat.Loss.Scale<-RNA.Loss.Scale.result[order(RNA.Loss.Scale.result$p_value),]
R.3cat.Loss.Scale1<- subset(R.3cat.Loss.Scale, source %in% c("CORUM" ,"GO:BP", "GO:MF", "WP", "KEGG", "REAC"))
R.3cat.Loss.Scale2<-R.3cat.Loss.Scale1[c(1,7,5, 9, 65),]#36, ribosome, rna process 11,12,58, 16,26,43,
R.3cat.Loss.Scale2$term_name<- factor(R.3cat.Loss.Scale2$term_name, levels= R.3cat.Loss.Scale2$term_name)
R.3cat.Loss.Scale2$Termtype<-c("metabolic process","metabolic process","metabolic process",
                               "CORUM root", "cell cycle")

## anti-scaling terms RNA
R.3cat.Gain.AntiScale<-RNA.Gain.AntiScale.result[order(RNA.Gain.AntiScale.result$p_value),]
R.3cat.Gain.AntiScale1<- subset(R.3cat.Gain.AntiScale, source %in% c("CORUM" ,"GO:BP", "GO:MF", "WP", "KEGG", "REAC"))
R.3cat.Gain.AntiScale2<-R.3cat.Gain.AntiScale1[c(1,3,5,6,17),]#58 total,  extracellular(adhesion, matrix, cell movement), ribosome, 
R.3cat.Gain.AntiScale2$term_name<- factor(R.3cat.Gain.AntiScale2$term_name, levels= R.3cat.Gain.AntiScale2$term_name)
R.3cat.Gain.AntiScale2$Termtype<-c("Extracellular","Extracellular","Extracellular","Extracellular","Extracellular")

R.3cat.Loss.AntiScale<-RNA.Loss.AntiScale.result[order(RNA.Loss.AntiScale.result$p_value),]
R.3cat.Loss.AntiScale1<- subset(R.3cat.Loss.AntiScale, source %in% c("CORUM" ,"GO:BP", "GO:MF", "WP", "KEGG", "REAC"))
R.3cat.Loss.AntiScale2<-R.3cat.Loss.AntiScale1[c(1,4,7,14,21),]#135 total, extracellular (adhesion, matrix, cell movement), 
R.3cat.Loss.AntiScale2$term_name<- factor(R.3cat.Loss.AntiScale2$term_name, levels= R.3cat.Loss.AntiScale2$term_name)
R.3cat.Loss.AntiScale2$Termtype<-c("Extracellular","Extracellular","Extracellular","Extracellular","Extracellular")


ggplot(R.3cat.Loss.Scale2, 
       aes(x=term_name, y=-log2(p_value), fill=Termtype))+
  geom_bar(stat="identity")+
  ylab("-log2 (P-Value)") +
  xlab("Key RNA terms")+
  ggtitle ("Biological terms enriched: RNA loss scale")+
  theme_classic()+
  ylim(0,160)+
  coord_flip()
#7x5
# plot.gprofiler.3Term.RNA.Loss.Scale

# Write those terms as csv files
#setwd()
#write.csv(RNA.Gain.Scale.result[,1:13], 
#          file =paste("gprofiler.RNA.Gain.3cat.Scale.csv",
#                      sep=','), row.names = TRUE)



###### High aneuploid cells vs low aneuploid cells, PROTEIN buffering key terms ####
# Sup Fig 5
#  PROTEIN Buffering 
# data generated in Protein_RNA_expression.PerCell.R and Protein_expression_data_GeneScore.R

### Run through g:profiler
## !! Run g:profiler below, or upload results from supplementary data **


### Find proteins buffered upon gain: 
HighAneu.Prot.Gain.Buffering.result<- gost(
  subset(CN.Diff.RNA.Prot_HighPloidy, Three.Protein.Gain=="Buffering")$RNA_Name,
  organism = "hsapiens",
  ordered_query = FALSE,
  multi_query = FALSE,
  significant = TRUE,
  exclude_iea = FALSE,
  measure_underrepresentation = FALSE,
  evcodes = FALSE,
  user_threshold = 0.05,
  correction_method = c("g_SCS"),
  domain_scope = "custom",
  custom_bg = as.character(CN.Diff.RNA.Prot_HighPloidy$RNA_Name),
  numeric_ns = "",
  sources = NULL,
  as_short_link = FALSE
)
HighAneu.Prot.Gain.Buffering.result<-HighAneu.Prot.Gain.Buffering.result$result


HighAneu.Prot.Loss.Buffering.result<- gost(
  subset(CN.Diff.RNA.Prot_HighPloidy, Three.Protein.Loss=="Buffering")$RNA_Name,
  organism = "hsapiens",
  ordered_query = FALSE,
  multi_query = FALSE,
  significant = TRUE,
  exclude_iea = FALSE,
  measure_underrepresentation = FALSE,
  evcodes = FALSE,
  user_threshold = 0.05,
  correction_method = c("g_SCS"),
  domain_scope = "custom",
  custom_bg = as.character(CN.Diff.RNA.Prot_HighPloidy$RNA_Name),
  numeric_ns = "",
  sources = NULL,
  as_short_link = FALSE
)
HighAneu.Prot.Loss.Buffering.result<-HighAneu.Prot.Loss.Buffering.result$result


## Proteins buffered upon gain, low aneuploid cells:
LowAneu.Prot.Gain.Buffering.result<- gost(
  subset(CN.Diff.RNA.Prot_LowPloidy, Three.Protein.Gain=="Buffering")$RNA_Name,
  organism = "hsapiens",
  ordered_query = FALSE,
  multi_query = FALSE,
  significant = TRUE,
  exclude_iea = FALSE,
  measure_underrepresentation = FALSE,
  evcodes = FALSE,
  user_threshold = 0.05,
  correction_method = c("g_SCS"),
  domain_scope = "custom",
  custom_bg = as.character(CN.Diff.RNA.Prot_LowPloidy$RNA_Name),
  numeric_ns = "",
  sources = NULL,
  as_short_link = FALSE
)
LowAneu.Prot.Gain.Buffering.result<-LowAneu.Prot.Gain.Buffering.result$result


## Proteins buffered upon loss:
LowAneu.Prot.Loss.Buffering.result<- gost(
  subset(CN.Diff.RNA.Prot_LowPloidy, Three.Protein.Loss=="Buffering")$RNA_Name,
  organism = "hsapiens",
  ordered_query = FALSE,
  multi_query = FALSE,
  significant = TRUE,
  exclude_iea = FALSE,
  measure_underrepresentation = FALSE,
  evcodes = FALSE,
  user_threshold = 0.05,
  correction_method = c("g_SCS"),
  domain_scope = "custom",
  custom_bg = as.character(CN.Diff.RNA.Prot_LowPloidy$RNA_Name),
  numeric_ns = "",
  sources = NULL,
  as_short_link = FALSE
)
LowAneu.Prot.Loss.Buffering.result<-LowAneu.Prot.Loss.Buffering.result$result

##Buffering terms HIGH ANEUPLOID CELL GROUP: 
P.Gain.Buff.HighAn<-HighAneu.Prot.Gain.Buffering.result[order(HighAneu.Prot.Gain.Buffering.result$p_value),]
P.Gain.Buff.HighAn1<- subset(P.Gain.Buff.HighAn, source %in% c("CORUM" ,"GO:BP", "GO:MF", "WP", "KEGG", "REAC"))
P.Gain.Buff.HighAn2<-P.Gain.Buff.HighAn1[c(1, 2,4,6, 14,17,54),]#168, CORUM root, Ribosomal 32, cell cycle, splisosome
P.Gain.Buff.HighAn2$term_name<- factor(P.Gain.Buff.HighAn2$term_name, levels= P.Gain.Buff.HighAn2$term_name)
P.Gain.Buff.HighAn2$Termtype<-c("CORUM root", "RNA processing","RNA processing","RNA processing","Ribosomal", "Ribosomal", "Ribosomal")

P.Loss.Buff.HighAn<-HighAneu.Prot.Loss.Buffering.result[order(HighAneu.Prot.Loss.Buffering.result$p_value),]
P.Loss.Buff.HighAn1<- subset(P.Loss.Buff.HighAn, source %in% c("CORUM" ,"GO:BP", "GO:MF", "WP", "KEGG", "REAC"))
P.Loss.Buff.HighAn2<-P.Loss.Buff.HighAn1[c(1, 2,4,10, 7,12,17),]#168, CORUM root, Ribosomal 32, cell cycle, splisosome
P.Loss.Buff.HighAn2$term_name<- factor(P.Loss.Buff.HighAn2$term_name, levels= P.Loss.Buff.HighAn2$term_name)
P.Loss.Buff.HighAn2$Termtype<-c("CORUM root", "RNA processing","RNA processing","RNA processing","Ribosomal", "Ribosomal", "Ribosomal" )

##Buffering terms LOW ANEUPLOID CELL GROUP: 
P.Gain.Buff.LowAn<-LowAneu.Prot.Gain.Buffering.result[order(LowAneu.Prot.Gain.Buffering.result$p_value),]
#P.Gain.Buff.LowAn1<- subset(P.Gain.Buff.LowAn, source %in% c("CORUM" ,"GO:BP", "GO:MF", "WP", "KEGG", "REAC"))
P.Gain.Buff.LowAn2<-P.Gain.Buff.LowAn[c(1,3, 2,6,23, 9),]#168, CORUM root, Ribosomal 32, cell cycle, splisosome
P.Gain.Buff.LowAn2$term_name<- factor(P.Gain.Buff.LowAn2$term_name, levels= P.Gain.Buff.LowAn2$term_name)
P.Gain.Buff.LowAn2$Termtype<-c("CORUM root", "CORUM root", "RNA processing","RNA processing","RNA processing", "Ribosomal")


P.Loss.Buff.LowAn<-LowAneu.Prot.Loss.Buffering.result[order(LowAneu.Prot.Loss.Buffering.result$p_value),]
#P.Loss.Buff.LowAn1<- subset(P.Loss.Buff.LowAn, source %in% c("CORUM" ,"GO:BP", "GO:MF", "WP", "KEGG", "REAC"))
P.Loss.Buff.LowAn2<-P.Loss.Buff.LowAn[c(1, 2, 3,4),]#168, CORUM root, Ribosomal 32, cell cycle, splisosome
P.Loss.Buff.LowAn2$term_name<- factor(P.Loss.Buff.LowAn2$term_name, levels= P.Loss.Buff.LowAn2$term_name)
P.Loss.Buff.LowAn2$Termtype<-c("RNA processing", "Ribosomal", "CORUM root", "CORUM root")

ggplot(P.Gain.Buff.HighAn2, 
       aes(x=term_name, y=-log2(p_value), fill=Termtype))+
  geom_bar(stat="identity")+
  ylab("-log2 (P-Value)") +
  xlab("Key RNA terms")+
  ggtitle ("Biological terms enriched: Protein gain Buffered")+
  theme_classic()+
  ylim(0,170)+
  coord_flip()
#7x5
# plot.gprofiler.HighAneuploid.Gain.Protein.Buffered

# Write those terms as csv files
setwd()
write.csv(P.Loss.Buff.LowAn[,1:13], 
          file =paste("gprofiler.LowAneu.Protein.Loss.Buffer.csv",
                      sep=','), row.names = TRUE)
