##### Factors affecting for protein buffering  ####
## Factors that influence gene attenuation at Protein/RNA level: 3 cat BUFFERING
## look up ROC Area under the curve
## 210429: Buffering factors for 3-category distribution: anti-scale, scale and buffering
## 210330
## By Klaske M. Schukken


## Protein expression data: 
##  https://www.cell.com/cell/fulltext/S0092-8674(19)31385-6#secsectitle0190
## Nusinow et al. Cell, 2020, Quantitative proteomics of the cancer cell line encyclopedia
## Chromosome arm data: Uri Ben-David paper 
## RNA expression data:  CCLE depmap data.


# Protein expression for genes with minimum 10 cells per category (gain & loss) (RNA & Protein)
# and only cells with RNA & DNA expression data. 
# from Protein_RNA_Corr_min10.Filtered.R

library('ggplot2')
library('tidyverse')
library('readxl')
library('reshape2')
library('BBmisc')
library('dplyr')
library("ROCR")
library("pROC") #for auc() calculations, and roc() curves
library("readODS")
library("ggVennDiagram")

# Steps: 
# 1) Get filtered data (RNA & Protein & aneuploidy)
#       only cells that have aneuploudy, RNA and Protein data
#       only genes with minimum of 10 cells per category. 
#       get data about location of each gene
#       get data about difference upon gain or loss, per gene
# 2) Get factor data from various sources 
#       extract MobiBD Disorder score per gene, etc.
# 3) Merge factor data with difference per gene
#       for gain and loss, for RNA and Protein
#       get ROC AUC and correlations (?) 
# 4) Plot difference (gain and loss) vs. factor, and ROC AUC

##### Step 1: Get filtered RNA and Protein data. ####
## Step 1: Get filtered RNA and Protein data. 

#### Protein expression data 
# filtered to only get data from cells that have RNA and Protein expression
## Protein info: 12755 genes
setwd()
Protein_Info3<-read_csv2("Protein_location_info.csv")
Protein_ProID<-read_csv2("Protein_ID_info.csv")
#Protein_ProID<-Protein_ProID[,-c(1)]# Remove "X" collumn

#Get updated Protein info with only 371 cell lines (cell lines also have RNA Data)
# protein data 12755 Proteins
Protein.Expression.filtered<-read.delim2("Protein_Expression_filtered.csv", 
                                         dec=",", header = TRUE, sep=";")


### RNA expression data 
# filtered to only get data from cells that have RNA and Protein expression
# Import RNA expression data (no aneuploidy data). CCLE depmap data.
# 371 cell lines
# 19 144 RNAs

# RNA Info. names, chrm location, etc. 
# 40 674 RNAs info 
RNA_Info<-read_excel("HGNC_Protein_Info_edited.xlsx", sheet ="HGNC names")
RNA_Info2<-RNA_Info[,c(2,3,11,12,22,14)]

# Get updated RNA expression with only 371 cell lines (cell lines also have RNA Data)
# 19 144 RNAs expression
RNA.Expression.filtered<-read.delim2("RNA_Expression_filtered.csv", 
                                     dec=",", header = TRUE, sep=";")


###  Aneuploidy data:
aneuploid<- read.csv("arm_calls_data_bendavid.csv", header=TRUE)


# Protein / RNA expression difference upon chrm arm gain/loss, per gene 
CN.Diff.xRNA.yProt.ThreeGroups <- read.csv("RNA.Protein_loss.neutral.gain_Difference_Pvalue_min10points_3cat.csv", header=TRUE) #9,414 genes
##!! or get this dataframe from supplementary data *! 

# Note: loss of about 2k genes (protein ID)-- 
# due to needing both RNA and Protein expression/gene
# and due to getting 10+ cell line with chrm gain/loss for chrm arm of gene location.

#add info about protein ID, uniprot ID, etc. 
CN.Diff.xRNA.yProt.ThreeGroups2<- merge(x=CN.Diff.xRNA.yProt.ThreeGroups, y=Protein_ProID, by.x="Protein_ID", by.y="Protein_Id")


###Optional: 
#download the factor results & scores, so you don't need to calculcate them yourself
# If you upload these datasets to R, you can go down to STEP 4: Plot ROC bargraph, mean expression boxplot
All.Factors.Diff<- read.csv(file= "Protein.AllFactors.csv")
Corr.protein.ROCAUC<- read.csv(file="Factor.ROCAUC.correlationscore.csv")

##### Step 2: Get potential buffering factor info ####
## Step 2: Get potential buffering factor info
## Protein Intrinsic Disorder data from: MobiBD. Version: 4.0 - Release: 2020_09
#  “mobidb_result.tsv”
#  https://mobidb.bio.unipd.it
#get Disorder data from MobiBD

#setwd()#set working directory to where you downloaded MobiBD data
Intrinsic_Disorder_All<-read.delim2("mobidb_result.tsv", 
                                         dec=".", header = TRUE, sep="\t")

Intrinsic_Disorder_MobiBD<- filter(Intrinsic_Disorder_All, feature=="prediction-disorder-mobidb_lite")#NS difference
#Intrinsic_Disorder_2<- filter(Intrinsic_Disorder_All, feature=="curated-disorder-merge")# NS
#Intrinsic_Disorder_3<- filter(Intrinsic_Disorder_All, feature=="prediction-disorder-vsl")# NS
#Intrinsic_Disorder_4<- filter(Intrinsic_Disorder_All, feature=="prediction-disorder-glo")# NS (RNA loss correlates, but don't trust, low corr)
#Intrinsic_Disorder_5<- filter(Intrinsic_Disorder_All, feature=="prediction-disorder-espN")# NS
lip_anchor<- filter(Intrinsic_Disorder_All, feature=="prediction-lip-anchor")#NS difference
Intrinsic_Disorder_binding<- filter(Intrinsic_Disorder_All, feature=="derived-binding_mode_disorder_to_disorder-mobi")#NS difference
Intrinsic_Disorder_bindingdo<- filter(Intrinsic_Disorder_All, feature=="derived-binding_mode_disorder_to_order-mobi")# Yes, but don't trust it
#Intrinsic_Disorder_transmembrane<- filter(Intrinsic_Disorder_All, feature=="prediction-transmembrane-uniprot")#NS
#Intrinsic_Disorder_signalPep<- filter(Intrinsic_Disorder_All, feature=="prediction-signal_peptide-uniprot")#NS
Intrinsic_Disorder_CoilCoil<- filter(Intrinsic_Disorder_All, feature=="prediction-coiled_coil-uniprot")#
Intrinsic_Disorder_LowComplex<- filter(Intrinsic_Disorder_All, feature=="prediction-low_complexity-merge")#
Intrinsic_Disorder_Polarity<- filter(Intrinsic_Disorder_All, feature=="prediction-polar-mobidb_lite_sub")#
Intrinsic_Disorder_Polyampholyte<- filter(Intrinsic_Disorder_All, feature=="prediction-polyampholyte-mobidb_lite_sub")#
Intrinsic_Disorder_homology<- filter(Intrinsic_Disorder_All, feature=="homology-domain-merge")#technically RNA corrolates, but I don't buy it, very small corr


# Note: acc is Uniprot_Acc

####
### Get mRNA decay rates 
## data from: Yang et al. 2003, Genome Res. Decay Rates of Human mRNAs: Correlation With Functional Characteristics and Sequence Attributes
## Sup Table 9. 
#setwd()#set working directory to where you downloaded this data
mRNA_Decay<-read_xlsx("/Users/user/Documents/Depmap Aneuploidy/Protein.Quantile/RNA_Protein_comparison/Buffering_factors/Yang.2003.mRNADecay.Rates.xlsx")
mRNA_Decay$Rate_1<- gsub(".\\s", "", mRNA_Decay$Rate_1)
mRNA_Decay$Rate_2<- gsub(".\\s", "", mRNA_Decay$Rate_2)
mRNA_Decay$Rate_3<- mRNA_Decay$Rate_1 #combine data from both rate columns
for (i in 1:length(mRNA_Decay$Accession)) {
  if (is.na(mRNA_Decay$Rate_3[i])){
    mRNA_Decay$Rate_3[i]<-mRNA_Decay$Rate_2[i] 
  }
}

headerRows<-c()# find rows with header titles instead of data
for (i in 1:length(mRNA_Decay$Accession)) {
  if (mRNA_Decay$Rate_3[i]=="Rate"){
    headerRows<-append(headerRows, c(i))
  }
}
mRNA_Decay2<-mRNA_Decay[-c(headerRows), ]#remove header rows
mRNA_Decay2<-mRNA_Decay2[,-c(3,4) ]#remove rate_1 and _2, already have _3
mRNA_Decay2$StdDev<-as.numeric(mRNA_Decay2$StdDev)
mRNA_Decay2$Rate_3<-as.numeric(mRNA_Decay2$Rate_3)
mRNA_Decay2$Accession<-as.factor(mRNA_Decay2$Accession)

#write.table(mRNA_Decay2$Accession, "mRNA.Decay_RefSeq.csv")

mRNADecay_info<- read_delim("/Users/user/Documents/Depmap Aneuploidy/Protein.Quantile/RNA_Protein_comparison/Buffering_factors/uniprot_mRNA.Decay_info.tab", 
                                delim="\t") 

mRNA_Decay3<-merge(x=mRNA_Decay2, y=mRNADecay_info, 
                      by.x="Accession", by.y="yourlist:M202104065C475328CEF75220C360D524E9D456CE1A9041H") # genes

#write_csv(mRNA_Decay3, "mRNA_Decay_Mean_Yang.etal.csv")
#mRNA_Decay3<-read_csv(file="mRNA_Decay_Mean_Yang.etal.csv")



###
## NCBI genome browser data. Human RefSeq data. 
## NOTE:! multiple 5'UTR and/or 3'UTR per gene are possible. 
##      keep UTR data seperate from other factor lists, otherwise it skews data for other genes 
## Downloaded 210409
## 5'UTR and 3'UTR length data
## https://genome.ucsc.edu/cgi-bin/hgTables?hgsid=1081368649_vYvfytavAH0WudzGAUbcIpcRFbCA

#bin            Indexing field to speed chromosome range queries.
#name	          Name of gene (usually transcript_id from GTF)
#chrom	  	    Reference sequence chromosome or scaffold
#strand   	    + or - for strand
#txStart	    	Transcription start position (or end position for minus strand item)
#txEnd	        Transcription end position (or start position for minus strand item)
#cdsStart	      Coding region start (or end position for minus strand item)
#cdsEnd	        Coding region end (or start position for minus strand item)
#exonCount	    Number of exons
#exonStarts	    Exon start positions (or end positions for minus strand item)
#exonEnds	 	    Exon end positions (or start positions for minus strand item)
#score		      score
#name2		      Alternate name (e.g. gene_id from GTF)
#cdsStartStat	  Status of CDS start annotation (none, unknown, incomplete, or complete)
#cdsEndStat	    Status of CDS end annotation (none, unknown, incomplete, or complete)
#exonFrames		  Exon frame {0,1,2}, or -1 if no frame for exon

#172766 genes. from NCBI.
NCBI.genedata<- read_tsv("NCBI.Human.RefSeq.tsv")
#NCBI.genedata$UTR5<-NA
#NCBI.genedata$UTR3<-NA

for (i in 1:length(NCBI.genedata$name)){
  if (NCBI.genedata$strand[i]=="+"){ #if positive strand
    NCBI.genedata$UTR5[i]<-abs(NCBI.genedata$cdsStart[i]-NCBI.genedata$txStart[i])
    NCBI.genedata$UTR3[i]<-abs(NCBI.genedata$txEnd[i]-NCBI.genedata$cdsEnd[i])
  } else { #if negaitve strand
    NCBI.genedata$UTR5[i]<-abs(NCBI.genedata$txEnd[i]-NCBI.genedata$cdsEnd[i])
    NCBI.genedata$UTR3[i]<-abs(NCBI.genedata$cdsStart[i]-NCBI.genedata$txStart[i])
  }
}

#write_csv(NCBI.genedata, "NCBI.GeneData.3UTR.5UTR.csv")
#NCBI.genedata<- read_csv("NCBI.GeneData.3UTR.5UTR.csv")

NCBI.genedata.NR<-subset(NCBI.genedata, startsWith(NCBI.genedata$name, "NR_")) #only munually curated, non coding genes
NCBI.genedata.NM<-subset(NCBI.genedata, startsWith(NCBI.genedata$name, "NM_")) #only munually curated, mRNA genes
NCBI.genedata.NM.5UTR<-unique(NCBI.genedata.NM[,c(13,17)]) #5' UTR data, unique entries, 33k
NCBI.genedata.NM.3UTR<-unique(NCBI.genedata.NM[,c(13,18)])#3'UTR data, unique entries, 26k




### Protein half-life data from: Methieson et al. 2018 Systematic analysis of protein turnover in primary cells
# 41467_2018_3106_MOESM5_ESM.xlsx
#setwd()#set working directory to where you downloaded this data
ProteinHL<-read_xlsx("/Users/user/Documents/Depmap Aneuploidy/Protein.Quantile/RNA_Protein_comparison/Buffering_factors/41467_2018_3106_MOESM5_ESM.xlsx")
ProteinHL$Mean_HalfLife<-NA

for (i in 1:length(ProteinHL$Mean_HalfLife)) {
ProteinHL$Mean_HalfLife[i]<-mean(c(ProteinHL$`Bcells replicate 1 half_life`[i], 
                                 ProteinHL$`Bcells replicate 2 half_life`[i], 
                              ProteinHL$`NK cells replicate 1 half_life`[i], 
                              ProteinHL$`NK cells replicate 2 half_life`[i],
                              ProteinHL$`Hepatocytes replicate 1 half_life`[i], 
                              ProteinHL$`Hepatocytes replicate 2 half_life`[i], 
                              ProteinHL$`Monocytes replicate 1 half_life`[i], 
                              ProteinHL$`Monocytes replicate 2 half_life`[i], 
                              ProteinHL$`Mouse Neurons, replicate 3 half_life`[i], 
                              ProteinHL$`Mouse Neurons, replicate 4 half_life`[i]), 
                              na.rm=TRUE)
}

#write_csv(ProteinHL, "ProteinHalfLife.Methieson.csv")

ProteinHL<- read_csv("ProteinHalfLife.Methieson.csv")


#### Human transcription and translation rates, calculated. 
## from: Hausser, Mayo, Keren, &Alon. 2019. 
# Central dogma rates and the trade-off between precision and economy 
# in gene expression. Nature communications, 10 (68)
# trancription rate is equal to mRNA abundance* rate of mRNA decay, directly proportional.
# just do mRNA abundance.  mRNA abundance= transcription * decay

Central.Dogma.info<- read_ods("/Users/user/Documents/Depmap Aneuploidy/Protein.Quantile/RNA_Protein_comparison/Buffering_factors/human.Transcription.Translation.Rates.ods") 
# variable                                                description                unit
# RNA.RPKM  log10 RPKM in RNAseq experiment of Eichhorn et al. (HeLa) reads / nucleotides
# RPF.RPKM      log10 RPKM in RP experiment of Eichhorn et al. (HeLa) reads / nucleotides
#        m       log10 mRNA abundance (estimated from Eichorn et al.)     copies per cell
#       bm  log10 transcription rate (estimated from Eichhorn et al.)              mRNA/h
#       bp    log10 translation rate (estimated from Eichhorn et al.)    protein/(mRNA.h)
#       lp                        log10 protein length (from UniProt)          aminoacids
#       lm log10 premRNA length (from RefSeq chromosomal coordinates)         nucleotides
#        p           log10 protein abundance estimated from bm and bp     proteins / cell

Central.Dogma<- read_ods("/Users/user/Documents/Depmap Aneuploidy/Protein.Quantile/RNA_Protein_comparison/Buffering_factors/human.Transcription.Translation.Rates.ods", 
                         sheet=2) 
colnames(Central.Dogma)<- c("EnsembleID", "RNA.RPKM", "RPF.RPKM", "m", "bm", "bp", "lp", "lm", "p")
# added ensemble_ID column
# now to add Uniprot_accession ID, and/or gene symbol to the Ensemble ID data
#write.table(Central.Dogma$EnsembleID, "/Users/user/Documents/Depmap Aneuploidy/Protein.Quantile/RNA_Protein_comparison/Buffering_factors/CentralDogma_EnsembleID.csv")
#write.table(Central.Dogma$EnsembleID[5200:8449], "/Users/user/Documents/Depmap Aneuploidy/Protein.Quantile/RNA_Protein_comparison/Buffering_factors/CentralDogma_EnsembleID2.csv")

#Ran emsemble_ID's through BioDB.net https://biodbnet-abcc.ncifcrf.gov/db/db2dbRes.php 
# got Gene IDs and Gene symbols and uniprot IDs for each ensemble ID. 

CD_EnsembleID_info<- read_delim("/Users/user/Documents/Depmap Aneuploidy/Protein.Quantile/RNA_Protein_comparison/Buffering_factors/bioDBnet_EnsembleID_Info.txt", 
                                delim="\t") #only 5200 ENSEMBLE ID's found back. from 8k genes in Central Dogma
CD_EnsembleID_info2<- read_delim("/Users/user/Documents/Depmap Aneuploidy/Protein.Quantile/RNA_Protein_comparison/Buffering_factors/bioDBnet_EnsembleID_Info2.txt", 
                                delim="\t") #only 5200 ENSEMBLE ID's found back. from 8k genes in Central Dogma
CD_EnsembleID_info3<- rbind(CD_EnsembleID_info, CD_EnsembleID_info2)

Central.Dogma2<-merge(x=Central.Dogma, y=CD_EnsembleID_info3, 
                      by.x="EnsembleID", by.y="Ensembl Gene ID") #8449 genes

####
#Intrinsic Protein variance, and intrinsic RNA variance
# calculated in Variance_aneuploidy.R
# calculated from filtered Depmap data (cells with both RNA & Protein data only)
# from cells not aneuploid for chrm arm gene is located on
# from Variance_aneuploidy.R
#setwd()#set working directory to where you downloaded this data
Protein.var.sd<-read.delim2("Protein.Variance.SD.VarMean.csv", 
                            dec=".", header = TRUE, sep=",")
RNA.var.sd<-read.delim2("RNA.Variance.SD.VarMean.csv", 
                        dec=".", header = TRUE, sep=",")

Protein.var.sd$log2.CoeffVar<- -log2(Protein.var.sd$Variance.Mean.noAn)
RNA.var.sd$log2.CoeffVar<- -log2(RNA.var.sd$Variance.Mean.noAn)

####
#Phosphosite plus 
#https://www.phosphosite.org/staticDownloads
# downloaded April 27, 2021
#Datasets: Acetylation, methylation, phosphorylation, regulatory, ubiquitination, sumoylation
#setwd()#set working directory to where you downloaded this data
#Acetyl
Acetylation<-read.csv("Acetylation_site_dataset2.csv")
Acetylation2<-subset(Acetylation, Acetylation$ORGANISM=="human")
Acetylation3<-table(Acetylation2$ACC_ID) #9505 proteins, 14880 ACC_IDs

Acetylation4<-data.frame(ACC_ID=rownames(Acetylation3), 
                         Acetylation=Acetylation3) #14880 genes, 7629 of which are 0
#Methyl
Methylation<-read_tsv("Methylation_site_dataset2.tsv")
Methylation2<-subset(Methylation, Methylation$ORGANISM=="human")
Methylation3<-table(Methylation2$ACC_ID) #ACC_IDs

Methylation4<-data.frame(ACC_ID=rownames(Methylation3), 
                         Methylation=Methylation3)#5691 genes, none are 0. 

#Phosphorylation
Phosphorylation<-read_tsv("Phosphorylation_site_dataset2.tsv")
Phosphorylation2<-subset(Phosphorylation, Phosphorylation$ORGANISM=="human")
Phosphorylation3<-table(Phosphorylation2$ACC_ID) #ACC_IDs

Phosphorylation4<-data.frame(ACC_ID=rownames(Phosphorylation3), 
                            Phosphorylation=Phosphorylation3)#19833 genes, of which none 0

#ubiquitination
ubiquitination<-read_tsv("ubiquitination_site_dataset2.tsv")
ubiquitination2<-subset(ubiquitination, ubiquitination$ORGANISM=="human")
ubiquitination3<-table(ubiquitination2$ACC_ID) #ACC_IDs

ubiquitination4<-data.frame(ACC_ID=rownames(ubiquitination3), 
                            ubiquitination=ubiquitination3)#12435 genes, of which none 0

#Sumoylation
Sumoylation<-read_tsv("Sumoylation_site_dataset2.tsv")
Sumoylation2<-subset(Sumoylation, Sumoylation$ORGANISM=="human")
Sumoylation3<-table(Sumoylation2$ACC_ID) #ACC_IDs

Sumoylation4<-data.frame(ACC_ID=rownames(Sumoylation3), 
                         Sumoylation=Sumoylation3)#2669 genes, of which none 0

#regulatory
regulatory<-read_tsv("regulatory_sites2.tsv")
regulatory2<-subset(regulatory, regulatory$ORGANISM=="human")
regulatory3<-table(regulatory2$ACC_ID) #ACC_IDs

regulatory4<-data.frame(ACC_ID=rownames(regulatory3), 
                        regulatory=regulatory3)#3110 genes, of which none 0

# did not end up using phosphosite plus data on galactose acytylation because not enough datapoints
# O-Gal-N-Ac

#O-GlcNAc


### Gene amplification percent ###
## From TSG.OG_CNV_Difference
## Percent of cells with gene amplifications(>2 DNA copy number)
## file from TSG.OG.CNV.Difference_v2.R
GeneAmplification<-read.csv(file="AmplificationRatio_perGene1.75.csv")

GeneAmplification$AmplificationRatio

### number of complexes gene is in: frequency in CORUM (all complexes dataset)
## Downloaded 210512
CORUM.Complex.all<-read_tsv(file="/Users/user/Documents/Depmap Aneuploidy/Protein.Quantile/RNA_Protein_comparison/Buffering_factors/CORUMallComplexes.tsv")

CORUM.Complex.all<-subset(CORUM.Complex.all, Organism=="Human")
CORUM.Complex.subunits<-as.vector(CORUM.Complex.all$`subunits(UniProt IDs)`)
CORUM.Complex.subunits<-strsplit(CORUM.Complex.subunits, ";")
CORUM.complex.subunits2 = c()
for (i in 1:length(CORUM.Complex.subunits)){
  CORUM.complex.subunits2<-append(CORUM.complex.subunits2, CORUM.Complex.subunits[[i]])
}
CORUM.complex.subunits3<-table(CORUM.complex.subunits2) #number of times genes occur in CORUM dataset
CORUM.complex.subunits4<-data.frame(ACC_ID=rownames(CORUM.complex.subunits3), 
                     gene=CORUM.complex.subunits3)#3110 genes, of which none 0


### number of Protein-Protein interactions
### HIPPIE version 2.2 updated 2/14/2009. downloaded 210513
## http://cbdm-01.zdv.uni-mainz.de/~mschaefer/hippie/download.php
## Get interactions with >0.6 confidence (0=none, 1=complete confidence)
HIPPIE.PPI<-read_tsv(file="HIPPIE-current.mitab.tsv")
HIPPIE.PPI$"Confidence Value"<-as.numeric(HIPPIE.PPI$"Confidence Value")
HIPPIE.PPI.2<- subset(HIPPIE.PPI, HIPPIE.PPI$`Confidence Value`>=0.6 ) #only get PPI with >0.6 confidence

HIPPIE.PPI.3<-table(HIPPIE.PPI.2$`Gene Name Interactor A`) #number of times genes occur in CORUM dataset
HIPPIE.PPI.3.2<-data.frame(PPI=HIPPIE.PPI.3)#3110 genes, of which none 0
HIPPIE.PPI.4<-table(HIPPIE.PPI.2$`Gene Name Interactor B`) #number of times genes occur in CORUM dataset
HIPPIE.PPI.4.2<-data.frame(PPI=HIPPIE.PPI.4)#3110 genes, of which none 0
HIPPIE.PPI.all<-merge(HIPPIE.PPI.3.2, HIPPIE.PPI.4.2, 
                      by.x="PPI.Var1", by.y="PPI.Var1", all=TRUE)
HIPPIE.PPI.all$PPI.Freq.x[is.na(HIPPIE.PPI.all$PPI.Freq.x)]<-0
HIPPIE.PPI.all$PPI.Freq.y[is.na(HIPPIE.PPI.all$PPI.Freq.y)]<-0
HIPPIE.PPI.all$PPI.Freq<-(as.numeric(HIPPIE.PPI.all$PPI.Freq.x)+ as.numeric(HIPPIE.PPI.all$PPI.Freq.y))

###
#Non-exponential decay
#https://www.sciencedirect.com/science/article/pii/S009286741631248X
#Kinetic Analysis of Protein Stability Reveals Age-Dependent Degradation
#Science direct
NonExponentialDecay<-read_xlsx("Non-exponential decay.xlsx")
NonExponentialDecay2<-NonExponentialDecay[,-c(5,6)] # Δ-score (Non-exponential decay delta), higher= non exponential decay
NonExponentialDecay2$`Δ-score (Non-exponential decay delta)`

###
#Dependency score, depmap
Dependency<-read_xlsx("Dependency scores.xlsx")
Dependency2<-Dependency[,-c(3)] # Dependency Score
Dependency2$`Dependency.Score`



###
# Aggregation score
# https://www.sciencedirect.com/science/article/pii/S2211124713005664
# Widespread Aggregation and Neurodegenerative Diseases Are Associated with Supersaturated Proteins

Aggregation<-read_xlsx("Aggregation score.xlsx")
Aggregation2<-Aggregation[,-c(4,5)] # Aggregation score
Aggregation2$`Aggregation score`

###
# Mutational score
# get all mutations in all cell lines
# subset only cells in our analysis
# cell line list from Protein_RNA_filtered_CellLine.R
# count table for gene mutations per gene
mutation<-read_csv("CCLE_mutations.csv")
mutations1<-mutation[,c(1,2,16)] # 
mutations1<- subset(mutations1, DepMap_ID %in% RNA_Protein_Cell_Info$DepMap_ID)
mutations2<- table(mutations1$Hugo_Symbol)
mutations2<-data.frame(Hugo_Symbol=rownames(mutations2), 
                        Mutations=mutations2)#18586 genes

table(mutation$'Variant_Classification')
mutations.subset<- subset(mutation, Variant_Classification %in% c("Nonsense_Mutation", "Missense_Mutation"))
mutations.subset1<-mutations.subset[,c(1,2,16)] # 
mutations.subset1<- subset(mutations.subset1, DepMap_ID %in% RNA_Protein_Cell_Info$DepMap_ID)
mutations.subset2<- table(mutations.subset1$Hugo_Symbol)
mutations.nonsense.missense<-data.frame(Hugo_Symbol=rownames(mutations.subset2), 
                       mutations.nonsense.missense=mutations.subset2)#18586 genes

mutations.subset<- subset(mutation, Variant_Classification %in% c("Frame_Shift_Ins", "Frame_Shift_Del"))
mutations.subset1<-mutations.subset[,c(1,2,16)] # 
mutations.subset1<- subset(mutations.subset1, DepMap_ID %in% RNA_Protein_Cell_Info$DepMap_ID)
mutations.subset2<- table(mutations.subset1$Hugo_Symbol)
mutations.FrameShift<-data.frame(Hugo_Symbol=rownames(mutations.subset2), 
                                        mutations.FrameShift=mutations.subset2)#18586 genes

mutations.nonsense.missense
mutations.FrameShift
### Nonsense_Mutation
## Protein gain ROC= 0.500
## Protein loss ROC= 0.504

### Silent
## Protein gain ROC= 0.509
## Protein loss ROC= 0.503

### Missence_mutation
## Protein gain ROC= 0.502
## Protein loss ROC= 0.510

## Frame_Shift_Del
## Protein gain ROC= 0.506
## Protein loss ROC= 0.517

## Frame_Shift_Ins
## Protein gain ROC= 0.522
## Protein loss ROC= 0.521

## De_novo_Start_OutOfFrame
## Protein gain ROC= 0.506
## Protein loss ROC= 0.520

## Splice_Site
## Protein gain ROC= 0.513
## Protein loss ROC= 0.501

## Start_Codon_SNP
## Protein gain ROC= 0.496
## Protein loss ROC= 0.528

## In_Frame_Del
## Protein gain ROC= 0.521
## Protein loss ROC= 0.509

## In_Frame_Ins
## Protein gain ROC= 0.490
## Protein loss ROC= 0.465

## De_novo_Start_OutOfFrame, Frame_Shift_Ins, Frame_Shift_Del, Nonsense_Mutation
## Protein gain ROC= 0.5035
## Protein loss ROC= 0.5123

## all
## Protein gain ROC= 0.5057
## Protein loss ROC= 0.5103

## missence nonsence
## Protein gain ROC= 0.5010
## Protein loss ROC= 0.5086

## frame shifts 
## Protein gain ROC= 0.5117
## Protein loss ROC= 0.5161


##### Step 3: Combine all factors into one dataset ####

# Step 1: Combine all factors into one dataset, add NA to non-info rows
#        Trim dataset: only keep useful collumns
#        Save dataset. 

#Factor.corr #9414 genes

colnames(Intrinsic_Disorder_MobiBD)<- c("acc", "feature", "start..end", "content_fraction", "Intrinsic_Disorder_MobiBD", "length")
colnames(lip_anchor)<- c("acc", "feature", "start..end", "content_fraction", "Loops.in.Protein(ANCHOR)", "length")
colnames(Intrinsic_Disorder_LowComplex)<- c("acc", "feature", "start..end", "content_fraction", "Intrinsic_Disorder_LowComplexity", "length")
colnames(Intrinsic_Disorder_Polarity)<- c("acc", "feature", "start..end", "content_fraction", "Intrinsic_Disorder_Polarity", "length")
colnames(Intrinsic_Disorder_Polyampholyte)<- c("acc", "feature", "start..end", "content_fraction", "Intrinsic_Disorder_Polyampholyte", "length")
colnames(Intrinsic_Disorder_homology)<- c("acc", "feature", "start..end", "content_fraction", "Intrinsic_Disorder_homology", "length")

#ProteinHL$Mean_HalfLife # 7125 data points
#Central.Dogma2$m # log10 abundance of mRNA, ~8k datapoints, 7029
#Central.Dogma2$p # log10 abundance of protein, ~8k datapoints
#Central.Dogma2$bm # log10 transcription rate, ~8k datapoints
#Central.Dogma2$bp # log10 translation rate, ~8k datapoints
#Central.Dogma2$lm # log10 pre-mRNA length, ~8k datapoints
#Central.Dogma2$lp # log10 protein length, ~8k datapoints
#mRNA_Decay3$Rate_3 # mRNA_Decay_rate, 2664 genes
#Protein.var.sd$log2.Protein.CoeffVar<-log2(Protein.var.sd$Variance.Mean.noAn) #intrinsic protein variance, 9256
#RNA.var.sd$log2.RNA.CoeffVar<-log2(RNA.var.sd$Variance.Mean.noAn) #instrinsic RNA variance, 9256
#GeneAmplification$AmplificationRatio #9401
#Acetylation4$Acetylation.Freq #4079, where not in database =0
#Methylation4$Methylation.Freq #3550, where not in database =0
#Phosphorylation4$Phosphorylation.Freq #8100, where not in database =0
#ubiquitination4$ubiquitination.Freq # 7035, where not in database =0
#Sumoylation4$Sumoylation.Freq # 1560, where not in database =0
#regulatory4$regulatory.Freq # 1754, where not in database =0
#CORUM.complex.subunits4$gene.Freq #2443, where not in database =0
#HIPPIE.PPI.all$PPI.Freq #8137 (confidence >0.6), where not in database =0
#Aggregation2$`Aggregation score`
#Dependency2$`Dependency.Score`
#NonExponentialDecay2$`Δ-score (Non-exponential decay delta)`
#mutations2$Hugo_Symbol
#mutations.FrameShift$Hugo_Symbol
#mutations.nonsense.missense$Hugo_Symbol

### !! Do not add 5'UTR and 3'UTR, they have multiple datapoints per gene and can mess up counts, correlations, etc. 


##Make dataframe
#All.Factors.Diff<-CN.Diff.xRNA.yProt.ThreeGroups2

## EDIT!!!:  add all factors one by one. merge via Gene_Symbol or Uniprot ID.
All.Factors.Diff<-merge(x=All.Factors.Diff, y= mutations.nonsense.missense, 
                        by.x="Gene_Symbol", by.y="Hugo_Symbol", 
                        sort = TRUE, all.x = TRUE)

## some factors only include genes if they have X, so all other genes are set to 0. 
All.Factors.Diff$gene.Freq[is.na(All.Factors.Diff$gene.Freq)]<-0 #all genes not in in CORUM dataset ==0 freq
All.Factors.Diff$Mutations.Freq[is.na(All.Factors.Diff$Mutations.Freq)]<-0 # all genes with no mutations, have mutation==0
All.Factors.Diff$mutations.FrameShift.Freq[is.na(All.Factors.Diff$mutations.FrameShift.Var1)]<-0 
All.Factors.Diff$mutations.nonsense.missense.Freq[is.na(All.Factors.Diff$mutations.nonsense.missense.Freq)]<-0 
All.Factors.Diff$PPI.Freq[is.na(All.Factors.Diff$PPI.Freq)]<-0 
All.Factors.Diff$ubiquitination.Freq[is.na(All.Factors.Diff$ubiquitination.Freq)]<-0 
All.Factors.Diff$Acetylation.Freq[is.na(All.Factors.Diff$Acetylation.Freq)]<-0 
All.Factors.Diff$Methylation.Freq[is.na(All.Factors.Diff$Methylation.Freq)]<-0 
All.Factors.Diff$Phosphorylation.Freq[is.na(All.Factors.Diff$Phosphorylation.Freq)]<-0 
All.Factors.Diff$Sumoylation.Freq[is.na(All.Factors.Diff$Sumoylation.Freq)]<-0 
All.Factors.Diff$regulatory.Freq[is.na(All.Factors.Diff$regulatory.Freq)]<-0 

All.Factors.Diff$Aggregation.score<- as.numeric(as.character(All.Factors.Diff$Aggregation.score)) #correct format

colnames(CN.Diff.xRNA.yProt.ThreeGroups2)
#write.csv(All.Factors.Diff, file= "Protein.AllFactors.csv")
## !! Or, you can download All factor values per gene from supplemental data: 
#All.Factors.Diff<- read.csv(file= "Protein.AllFactors.csv")

All.Factors.Diff.noAS.Gain<-subset(All.Factors.Diff, Protein.Diff.Gain> -0.1) # Buffer vs scaling
All.Factors.Diff.noAS.Loss<-subset(All.Factors.Diff, Protein.Diff.Loss< 0.1) # Buffer vs scaling

All.Factors.Diff.noS.Gain<-subset(All.Factors.Diff, Protein.Diff.Gain< 0.25) # Buffer vs anti-scaling
All.Factors.Diff.noS.Loss<-subset(All.Factors.Diff, Protein.Diff.Loss> -0.25) # Buffer vs anti-scaling

#Protein_ID: sp.Q9NQ94.A1CF_HUMAN  sp.P01023.A2MG_HUMAN  sp.A8K2U0.A2ML1_HUMAN
#RNA_ID: A1CF..29974.   A2M..2.        A2ML1..144568.
#RNA_Name: A1CF  A2M   A2ML1
#Uniprot_Acc: "A0AV96"   "A0AVF1"   "A0AVT1"   "A0FGR8-2"


##### Step 4: RUN ROC analysis as FUNCTION: calculate ROC AUC ####

## Define function for getting ROC data per dataset, per category: 


# make list of column names and what they are titles: 
NameTestCollumns= c("PPI.Freq", "gene.Freq", "regulatory.Freq", "Sumoylation.Freq", 
                    "ubiquitination.Freq", "Phosphorylation.Freq", "Methylation.Freq", "Acetylation.Freq", 
                    "AmplificationRatio", "log2.CoeffVar.RNA", "log2.Protein.CoeffVar", "Rate_3", 
                    "m", "p", "lp", "lm", "bm", "bp", 
                    "Mean_HalfLife", "Loops.in.Protein.ANCHOR.", "Intrinsic_Disorder_MobiBD", "Intrinsic_Disorder_LowComplexity", 
                    "Intrinsic_Disorder_Polarity", "Intrinsic_Disorder_Polyampholyte", "Intrinsic_Disorder_homology", 
                    "Dependency.Score", "Aggregation.score", "Non.exponential.decay.delta", "Mutations.Freq", "mutations.nonsense.missense.Freq", "mutations.nonsense.missense.Freq")

testdatanames= c("Protein-protein interaction", "Protein complex (CORUM)", "Protein regulatory sites", "Sumoylation sites", 
                 "Ubiquitination sites", "Phosphorylation sites", "Methylation sites", "Acetylation sites", 
                 "Percent gene amplification", "RNA neutral variance", "Protein neutral variance", "mRNA decay rate", 
                 "mRNA abundance", "Protein abundance", "Protein length", "mRNA length", "Transcription rate", "Translation rate", 
                 "Protein half life", "Loops in protein score", "Intrinsic protein disorder", "Low complexity score", 
                 "Protein polarity", "Protein polyampholyte score", "Homology score", 
                 "Dependency score", "Aggregation score", "Non-exponential decay delta", "Mutation count (all)", "Mutation count (frame shifts)", "Mutation count (nonsence and missense)")


dataset= All.Factors.Diff


FindROCAUC<- function(dataset, NameTestCollumns, testdatanames) {
  Corr.protein.loop<-data.frame(Data=as.character(), 
                                Chrm_CN=character(), 
                                Gene_Type=character(), 
                                Category=character(),
                                p_value=as.numeric(), 
                                Corr_coef=as.numeric(), 
                                significant=character(), 
                                ROCauc=as.numeric())
  
  
  for (i in 1:length(NameTestCollumns)){
    
    TestDataCol= NameTestCollumns[i]
    name= testdatanames[i]
    TestData<- as.numeric(as.character(dataset %>% pull(TestDataCol)))
    
    # Buffering, Gain: 
    Gain.Protein.corr<- cor.test(TestData, dataset$Protein.Diff.Gain, method="pearson")
    Corr.protein.loop<-rbind(Corr.protein.loop, 
                             data.frame(data=name, 
                                        Chrm_CN="Gain", 
                                        Gene_Type="Protein", 
                                        Category="Buffering",
                                        p_value=Gain.Protein.corr$p.value, 
                                        Corr_coef=Gain.Protein.corr$estimate, 
                                        significant=Gain.Protein.corr$p.value<0.05, 
                                        ROCauc=auc(dataset$Three.Protein.Gain=="Buffering", 
                                                   TestData, na.rm=TRUE ))) 
    
    # Buffering, Loss: 
    Loss.Protein.corr<-cor.test(TestData, dataset$Protein.Diff.Loss, method="pearson")
    Corr.protein.loop<-rbind(Corr.protein.loop, 
                             data.frame(data=name, 
                                        Chrm_CN="Loss", 
                                        Gene_Type="Protein", 
                                        Category="Buffering",
                                        p_value=Loss.Protein.corr$p.value, 
                                        Corr_coef=Loss.Protein.corr$estimate, 
                                        significant=Loss.Protein.corr$p.value<0.01, 
                                        ROCauc=auc(dataset$Three.Protein.Loss=="Buffering", 
                                                   TestData, na.rm=TRUE )))
    
    # AS, Gain: 
    Gain.Protein.corr<- cor.test(TestData, dataset$Protein.Diff.Gain, method="pearson")
    Corr.protein.loop<-rbind(Corr.protein.loop, 
                             data.frame(data=name, 
                                        Chrm_CN="Gain", 
                                        Gene_Type="Protein", 
                                        Category="Anti-Scaling",
                                        p_value=Gain.Protein.corr$p.value, 
                                        Corr_coef=Gain.Protein.corr$estimate, 
                                        significant=Gain.Protein.corr$p.value<0.05, 
                                        ROCauc=auc(dataset$Three.Protein.Gain=="Anti-Scaling", 
                                                   TestData, na.rm=TRUE ))) 
    
    # AS, Loss: 
    Loss.Protein.corr<-cor.test(TestData, dataset$Protein.Diff.Loss, method="pearson")
    Corr.protein.loop<-rbind(Corr.protein.loop, 
                             data.frame(data=name, 
                                        Chrm_CN="Loss", 
                                        Gene_Type="Protein", 
                                        Category="Anti-Scaling",
                                        p_value=Loss.Protein.corr$p.value, 
                                        Corr_coef=Loss.Protein.corr$estimate, 
                                        significant=Loss.Protein.corr$p.value<0.01, 
                                        ROCauc=auc(dataset$Three.Protein.Loss=="Anti-Scaling", 
                                                   TestData, na.rm=TRUE ))) 
    
    # Scaling, Gain: 
    Gain.Protein.corr<- cor.test(TestData, dataset$Protein.Diff.Gain, method="pearson")
    Corr.protein.loop<-rbind(Corr.protein.loop, 
                             data.frame(data=name, 
                                        Chrm_CN="Gain", 
                                        Gene_Type="Protein", 
                                        Category="Scaling",
                                        p_value=Gain.Protein.corr$p.value, 
                                        Corr_coef=Gain.Protein.corr$estimate, 
                                        significant=Gain.Protein.corr$p.value<0.05, 
                                        ROCauc=auc(dataset$Three.Protein.Gain=="Scaling", 
                                                   TestData, na.rm=TRUE ))) 
    
    # Scaling, Loss: 
    Loss.Protein.corr<-cor.test(TestData, dataset$Protein.Diff.Loss, method="pearson")
    Corr.protein.loop<-rbind(Corr.protein.loop, 
                             data.frame(data=name, 
                                        Chrm_CN="Loss", 
                                        Gene_Type="Protein", 
                                        Category="Scaling",
                                        p_value=Loss.Protein.corr$p.value, 
                                        Corr_coef=Loss.Protein.corr$estimate, 
                                        significant=Loss.Protein.corr$p.value<0.01, 
                                        ROCauc=auc(dataset$Three.Protein.Loss=="Scaling", 
                                                   TestData, na.rm=TRUE ))) 
    
  }
  return(Corr.protein.loop)
}
Corr.protein.ROCAUC<- FindROCAUC(dataset, NameTestCollumns, testdatanames)


## Now look at mean factor score per category: 
FindMeanFactorPerCategory<- function(dataset, NameTestCollumns, testdatanames) {
  Mean.Factor<-data.frame(Factor=as.character(),
                          Gain.AntiScaling= as.numeric(), 
                          Gain.Buffering= as.numeric(),  
                          Gain.Scaling= as.numeric(), 
                          Loss.AntiScaling= as.numeric(), 
                          Loss.Buffering= as.numeric(), 
                          Loss.Scaling= as.numeric() )
  
  for (i in 1:length(NameTestCollumns)){
    TestColName= NameTestCollumns[i]
    TestColData= dataset %>% pull(TestColName)
    name= testdatanames[i]
    # get mean expression of variable per group, 
    # normalize to 0-1 scale: (value-minimum)/max-minimum
    
    table.Gain <- data.frame( tapply(TestColData, dataset$Three.Protein.Gain, mean, na.rm=TRUE) )
    table.Gain$StandardData<- table.Gain[,1] 
    
    table.Loss <- data.frame( tapply(TestColData, dataset$Three.Protein.Loss, mean, na.rm=TRUE) )
    table.Loss$StandardData<- table.Loss[,1]
    
    Mean.Factor<-rbind(Mean.Factor, 
                       data.frame(Factor=name, 
                                  Gain.AntiScaling= table.Gain[1,2], 
                                  Gain.Buffering= table.Gain[2,2], 
                                  Gain.Scaling= table.Gain[3,2], 
                                  Loss.AntiScaling= table.Loss[1,2], 
                                  Loss.Buffering= table.Loss[2,2], 
                                  Loss.Scaling= table.Loss[3,2] ) )
  }
  return(Mean.Factor)
} # not standardized, just give values
FindMeanFactorPerCategory2<- function(dataset, NameTestCollumns, testdatanames) {
  Mean.Factor2<-data.frame(Factor=character(), 
                           Gain.AntiScaling= as.numeric(), 
                           Gain.Buffering= as.numeric(),  
                           Gain.Scaling= as.numeric(), 
                           Loss.AntiScaling= as.numeric(), 
                           Loss.Buffering= as.numeric(), 
                           Loss.Scaling= as.numeric() )
  
  for (i in 1:length(NameTestCollumns)){
    TestColName= NameTestCollumns[i]
    TestColData= dataset %>% pull(TestColName)
    name= testdatanames[i]
    # get mean expression of variable per group, 
    # normalize to 0-1 scale: (value-minimum)/max-minimum
    
    table.Gain <- data.frame( tapply(TestColData, dataset$Three.Protein.Gain, mean, na.rm=TRUE) )
    table.Gain$StandardData<- (table.Gain[,1] - mean(TestColData, na.rm=TRUE) )/ sd(TestColData, na.rm=TRUE)
    
    table.Loss <- data.frame( tapply(TestColData, dataset$Three.Protein.Loss, mean, na.rm=TRUE) )
    table.Loss$StandardData<- (table.Loss[,1] - mean(TestColData, na.rm=TRUE) )/ sd(TestColData, na.rm=TRUE)
    
    Mean.Factor2<-rbind(Mean.Factor2, 
                        data.frame(Factor=name, 
                                   Gain.AntiScaling= table.Gain[1,2], 
                                   Gain.Buffering= table.Gain[2,2], 
                                   Gain.Scaling= table.Gain[3,2], 
                                   Loss.AntiScaling= table.Loss[1,2], 
                                   Loss.Buffering= table.Loss[2,2], 
                                   Loss.Scaling= table.Loss[3,2] ) )
  }
  return(Mean.Factor2)
}# Standardized: x minus mean, divided by standard dev. 
##
Mean.Factor.PerCategory.notNorm<- FindMeanFactorPerCategory(dataset, NameTestCollumns, testdatanames)
#Mean.Factor.PerCategory<- FindMeanFactorPerCategory2(dataset, NameTestCollumns, testdatanames)


#setwd()#set working directory to where you downloaded this data
write.csv(Mean.Factor.PerCategory.notNorm, file="Mean.Factor.PerCategory.notStandardized")
#write.csv(Mean.Factor.PerCategory, file="Mean.Factor.PerCategory.Standardized")



#####         Plot ROC bargraph, mean expression boxplot ####

### boxplots of mean factor score/value in Buffered, scaling and antiscaling genes upon either gain and loss

y.lab="Protein Protein interactions" #Name of factor you are looking at
testcollumn<- All.Factors.Diff$PPI.Freq #get the collumn of factor you want to look at


## Chrm gain categories
ggplot(All.Factors.Diff, aes(x=Three.Protein.Gain, y=testcollumn, fill=Three.Protein.Gain))+
  geom_boxplot(outlier.shape = NA)+
  scale_fill_manual(values=c("darkorange2", "Gold2", "palegreen4"))+
  theme_classic()+
  ggtitle("Gain")+
  geom_hline(yintercept=0)+
  coord_cartesian(ylim=c(0,150))+
  xlab("")+
  ylab(y.lab)
# 5x4
# plot.mean.factor.Gain.boxplot.NED

#Chrm loss categories
ggplot(All.Factors.Diff, aes(x=Three.Protein.Loss, y=testcollumn, fill=Three.Protein.Loss))+
  geom_boxplot(outlier.shape = NA)+
  scale_fill_manual(values=c("darkorange2", "Gold2", "palegreen4"))+
  theme_classic()+
  ggtitle("Loss")+
  geom_hline(yintercept=0)+
  coord_cartesian(ylim=c(0,150))+
  xlab("")+
  ylab(y.lab)
# 5x4
# plot.mean.factor.Loss.boxplot.NED

G.AS<- subset(All.Factors.Diff, Three.Protein.Gain=="Anti-Scaling") #chrm gain, category anti scaling
G.B<- subset(All.Factors.Diff, Three.Protein.Gain=="Buffering") #chrm gain, category buffered
G.S<- subset(All.Factors.Diff, Three.Protein.Gain=="Scaling") #chrm gain, category scaling
L.AS<- subset(All.Factors.Diff, Three.Protein.Loss=="Anti-Scaling")
L.B<- subset(All.Factors.Diff, Three.Protein.Loss=="Buffering")
L.S<- subset(All.Factors.Diff, Three.Protein.Loss=="Scaling")

## Test significance between groups
# change the collumn to the factor you want to test. 
t.test(G.AS$PPI.Freq, G.B$PPI.Freq)
t.test(G.S$PPI.Freq,  G.B$PPI.Freq)
t.test(L.AS$PPI.Freq, L.B$PPI.Freq)
t.test(L.S$PPI.Freq,  L.B$PPI.Freq)



### Plot ROC bargraphs:  ###

## But first Add 5' and 3' UTR  AUC to data ***
NCBI.genedata.NM.5UTR$UTR5 #5' UTR data, unique entries, 33k
TestData<- NCBI.genedata.NM.5UTR
name<- "5' UTR length" 

Diff.test<-merge(x=TestData, y= CN.Diff.xRNA.yProt.ThreeGroups2, 
                 by.x="name2", by.y="Gene_Symbol", 
                 sort = TRUE)
TestDataCol<-Diff.test$UTR5

Gain.Protein.corr<- cor.test(TestDataCol, Diff.test$Protein.Diff.Gain, method="pearson")
Corr.protein.ROCAUC<-rbind(Corr.protein.ROCAUC, 
                                   data.frame(data=name, 
                                              Chrm_CN="Gain", 
                                              Gene_Type="Protein", 
                                              Category= "Buffering",
                                              p_value=Gain.Protein.corr$p.value, 
                                              Corr_coef=Gain.Protein.corr$estimate, 
                                              significant=Gain.Protein.corr$p.value<0.05, 
                                              ROCauc=auc(Diff.test$Three.Protein.Gain=="Buffering", 
                                                         TestDataCol ))) 
Loss.Protein.corr<- cor.test(TestDataCol, Diff.test$Protein.Diff.Loss, method="pearson")
Corr.protein.ROCAUC<-rbind(Corr.protein.ROCAUC, 
                                   data.frame(data=name, 
                                              Chrm_CN="Loss", 
                                              Gene_Type="Protein", 
                                              Category= "Buffering",
                                              p_value=Loss.Protein.corr$p.value, 
                                              Corr_coef=Loss.Protein.corr$estimate, 
                                              significant=Loss.Protein.corr$p.value<0.01, 
                                              ROCauc=auc(Diff.test$Three.Protein.Loss=="Buffering", 
                                                         TestDataCol ))) 


## Now add 3' UTR
NCBI.genedata.NM.3UTR$UTR3#3'UTR data, unique entries, 26k
TestData<- NCBI.genedata.NM.3UTR
name<- "3' UTR length" 

Diff.test<-merge(x=TestData, y= CN.Diff.xRNA.yProt.ThreeGroups2, 
                 by.x="name2", by.y="Gene_Symbol", 
                 sort = TRUE)
TestDataCol<-Diff.test$UTR3

Gain.Protein.corr<- cor.test(TestDataCol, Diff.test$Protein.Diff.Gain, method="pearson")
Corr.protein.ROCAUC<-rbind(Corr.protein.ROCAUC, 
                                   data.frame(data=name, 
                                              Chrm_CN="Gain", 
                                              Gene_Type="Protein", 
                                              Category= "Buffering",
                                              p_value=Gain.Protein.corr$p.value, 
                                              Corr_coef=Gain.Protein.corr$estimate, 
                                              significant=Gain.Protein.corr$p.value<0.05, 
                                              ROCauc=auc(Diff.test$Three.Protein.Gain=="Buffering", 
                                                         TestDataCol ))) 
Loss.Protein.corr<- cor.test(TestDataCol, Diff.test$Protein.Diff.Loss, method="pearson")
Corr.protein.ROCAUC<-rbind(Corr.protein.ROCAUC, 
                                   data.frame(data=name, 
                                              Chrm_CN="Loss", 
                                              Gene_Type="Protein", 
                                              Category= "Buffering",
                                              p_value=Loss.Protein.corr$p.value, 
                                              Corr_coef=Loss.Protein.corr$estimate, 
                                              significant=Loss.Protein.corr$p.value<0.01, 
                                              ROCauc=auc(Diff.test$Three.Protein.Loss=="Buffering", 
                                                         TestDataCol ))) 

write.csv(Corr.protein.ROCAUC, file="Factor.ROCAUC.correlationscore.csv")
## ! OR you can download the file from supplementary: 
#Corr.protein.ROCAUC<- read.csv(file="Factor.ROCAUC.correlationscore.csv")

## now that 5' and 3 UTR have been added, remove dataset specific factors: 
Corr.protein.noHighRNAVar2<-Corr.protein.ROCAUC
Corr.protein.noHighRNAVar2<-Corr.protein.noHighRNAVar2[-c(49:66, 151:156, 169:186),] #remove RNA/Protein neutral variance, & gene amp ratio


#Subset data by type: 
Corr.protein.noHighRNAVar2.Gain<- subset(Corr.protein.noHighRNAVar2, Chrm_CN=="Gain" & Gene_Type=="Protein" & Category=="Buffering")
Corr.protein.noHighRNAVar2.Gain<- Corr.protein.noHighRNAVar2.Gain[order(Corr.protein.noHighRNAVar2.Gain$ROCauc),]

Corr.protein.noHighRNAVar2.Loss<- subset(Corr.protein.noHighRNAVar2, Chrm_CN=="Loss" & Gene_Type=="Protein" & Category=="Buffering")
Corr.protein.noHighRNAVar2.Loss<- Corr.protein.noHighRNAVar2.Loss[order(Corr.protein.noHighRNAVar2.Loss$ROCauc),]


## Buffering ROC AUC bargraph
Corr.protein.noHighRNAVar2.Gain$data<- factor(Corr.protein.noHighRNAVar2.Gain$data, level=Corr.protein.noHighRNAVar2.Gain$data)
Corr.protein.noHighRNAVar2.Loss$data<- factor(Corr.protein.noHighRNAVar2.Loss$data, level=Corr.protein.noHighRNAVar2.Gain$data)

ggplot(Corr.protein.noHighRNAVar2.Gain, aes(x=data, y=ROCauc))+ #5x5 figure
  geom_bar(stat="identity")+
  ylab("ROC AUC")+
  xlab("")+ 
  theme_classic()+
  geom_hline(yintercept=0.5, color="black")+
  coord_flip(ylim=c(0.45, 0.6))+
  ggtitle("Buffering upon chrm gain, factors")
# 5x4
# plot.ROCauc.Protein.gain.Buffer_v4

ggplot(Corr.protein.noHighRNAVar2.Loss, aes(x=data, y=ROCauc))+ #5x5 figure
  geom_bar(stat="identity")+
  ylab("ROC AUC")+
  xlab("")+ 
  theme_classic()+
  geom_hline(yintercept=0.5, color="black")+
  coord_flip(ylim=c(0.45, 0.6))+
  ggtitle("Buffering upon chrm loss factors")
# 5x4
# plot.ROCauc.Protein.loss.Buffer_v4



#### plot dataset specific factors##
Corr.protein.ROCAUC.dataset<-Corr.protein.ROCAUC
Corr.protein.ROCAUC.dataset<-Corr.protein.ROCAUC.dataset[c(55:66, 151:155, 169:186),] #get RNA/Protein neutral variance, & gene amp ratio

#Subset data by type: 
Corr.protein.ROCAUC.dataset.Gain<- subset(Corr.protein.ROCAUC.dataset, Chrm_CN=="Gain" & Gene_Type=="Protein" & Category=="Buffering")
Corr.protein.ROCAUC.dataset.Gain<- Corr.protein.ROCAUC.dataset.Gain[order(Corr.protein.ROCAUC.dataset.Gain$ROCauc),]

Corr.protein.ROCAUC.dataset.Loss<- subset(Corr.protein.ROCAUC.dataset, Chrm_CN=="Loss" & Gene_Type=="Protein" & Category=="Buffering")
Corr.protein.ROCAUC.dataset.Loss<- Corr.protein.ROCAUC.dataset.Loss[order(Corr.protein.ROCAUC.dataset.Loss$ROCauc),]

## Buffering ROC AUC bargraph
Corr.protein.ROCAUC.dataset.Gain$data<- factor(Corr.protein.ROCAUC.dataset.Gain$data, level=Corr.protein.ROCAUC.dataset.Gain$data)
Corr.protein.ROCAUC.dataset.Loss$data<- factor(Corr.protein.ROCAUC.dataset.Loss$data, level=Corr.protein.ROCAUC.dataset.Gain$data)

ggplot(Corr.protein.ROCAUC.dataset.Gain, aes(x=data, y=ROCauc))+ #5x4 figure
  geom_bar(stat="identity")+
  ylab("ROC AUC")+
  xlab("")+ 
  theme_classic()+
  geom_hline(yintercept=0.5, color="black")+
  coord_flip(ylim=c(0.45, 0.65))+
  ggtitle("Buffering upon chrm gain, dataset specific factors")
# 5x4
# plot.ROCauc.Protein.gain.Buffer.datasetspecific_v2

ggplot(Corr.protein.ROCAUC.dataset.Loss, aes(x=data, y=ROCauc))+ #5x4 figure
  geom_bar(stat="identity")+
  ylab("ROC AUC")+
  xlab("")+ 
  theme_classic()+
  geom_hline(yintercept=0.5, color="black")+
  coord_flip(ylim=c(0.45, 0.65))+
  ggtitle("Buffering upon chrm loss, dataset specific factors")
# 5x4
# plot.ROCauc.Protein.loss.Buffer.datasetspecific_v2






### Anti Scaling ##
#Corr.protein.noHighRNAVar2<-Corr.protein.noHighRNAVar

Corr.protein.AS.Gain<- subset(Corr.protein.noHighRNAVar2, Chrm_CN=="Gain" & Gene_Type=="Protein" & Category=="Anti-Scaling")
Corr.protein.AS.Gain<- Corr.protein.AS.Gain[order(Corr.protein.AS.Gain$ROCauc),]

Corr.protein.AS.Loss<- subset(Corr.protein.noHighRNAVar2, Chrm_CN=="Loss" & Gene_Type=="Protein" & Category=="Anti-Scaling")
Corr.protein.AS.Loss<- Corr.protein.AS.Loss[order(Corr.protein.AS.Loss$ROCauc),]

ggplot(Corr.protein.AS.Gain, aes(x=reorder(data, ROCauc), y=ROCauc))+ #5x5 figure
  geom_bar(stat="identity")+
  ylab("ROC AUC")+
  xlab("")+ 
  theme_classic()+
  geom_hline(yintercept=0.5, color="black")+
  coord_flip(ylim=c(0.40, 0.65))+
  ggtitle("Anti-Scaling upon chrm gain, factors")
# 5x4
# plot.ROCauc.Protein.gain.AS_v1

ggplot(Corr.protein.AS.Loss, aes(x=reorder(data, ROCauc), y=ROCauc))+ #5x5 figure
  geom_bar(stat="identity")+
  ylab("ROC AUC")+
  xlab("")+ 
  theme_classic()+
  geom_hline(yintercept=0.5, color="black")+
  coord_flip(ylim=c(0.40, 0.65))+
  ggtitle("Anti-Scaling upon chrm loss factors")
# 5x4
# plot.ROCauc.Protein.loss.AS_v1




###Scaling ###
#Corr.protein.noHighRNAVar2<-Corr.protein.noHighRNAVar

Corr.protein.S.Gain<- subset(Corr.protein.noHighRNAVar2, Chrm_CN=="Gain" & Gene_Type=="Protein" & Category=="Scaling")
Corr.protein.S.Gain<- Corr.protein.S.Gain[order(Corr.protein.S.Gain$ROCauc),]

Corr.protein.S.Loss<- subset(Corr.protein.noHighRNAVar2, Chrm_CN=="Loss" & Gene_Type=="Protein" & Category=="Scaling")
Corr.protein.S.Loss<- Corr.protein.S.Loss[order(Corr.protein.S.Loss$ROCauc),]

ggplot(Corr.protein.S.Gain, aes(x=reorder(data, ROCauc), y=ROCauc))+ #5x5 figure
  geom_bar(stat="identity")+
  ylab("ROC AUC")+
  xlab("")+ 
  theme_classic()+
  geom_hline(yintercept=0.5, color="black")+
  coord_flip(ylim=c(0.44, 0.6))+
  ggtitle("Scaling upon chrm gain, factors")
# 5x4
# plot.ROCauc.Protein.gain.Scaling_v1

ggplot(Corr.protein.S.Loss, aes(x=reorder(data, ROCauc), y=ROCauc))+ #5x5 figure
  geom_bar(stat="identity")+
  ylab("ROC AUC")+
  xlab("")+ 
  theme_classic()+
  geom_hline(yintercept=0.5, color="black")+
  coord_flip(ylim=c(0.44, 0.6))+
  ggtitle("Scaling upon chrm loss factors")
# 5x4
# plot.ROCauc.Protein.loss.Scaling_v1



#####         PLOT ROC AUC curve for : buffer, Anti-Scaling, and Scale ####
# automatically plot all data. Put pdfs in below folder: 
#setwd()#set working directory to where you downloaded this data

### EDIT! change name and test collumn you want to test as needed: 
## Get ROC AUC for specific factor of interest. 

name <- "NED delta" #name of factor
TestDataCol <- All.Factors.Diff$Non.exponential.decay.delta #factor you want to see ROC curve of
# Category== Buffering
# Automatically make all the plots. 
# 3x3
pdf(file = paste0("Protein.Gain.ROC.Buffer.", name, ".pdf"),
    width = 3, 
    height = 3)
roc(All.Factors.Diff$Three.Protein.Gain=="Buffering", TestDataCol, 
    smoothed = TRUE,
    plot=TRUE, auc.polygon=FALSE, max.auc.polygon=FALSE, grid=FALSE,
    print.auc=TRUE, show.thres=TRUE)
dev.off()


pdf(file = paste0("Protein.Loss.ROC.Buffer.", name, ".pdf"),
    width = 3, 
    height = 3)
roc(All.Factors.Diff$Three.Protein.Loss=="Buffering", TestDataCol,
    smoothed = TRUE,
    plot=TRUE, auc.polygon=FALSE, max.auc.polygon=FALSE, grid=FALSE,
    print.auc=TRUE, show.thres=TRUE)
dev.off()




#Anti scaling
pdf(file = paste0("Protein.Gain.ROC.AS.", name, ".pdf"),
    width = 3, 
    height = 3)
roc(All.Factors.Diff$Three.Protein.Gain=="Anti-Scaling", TestDataCol, 
    smoothed = TRUE,
    plot=TRUE, auc.polygon=FALSE, max.auc.polygon=FALSE, grid=FALSE,
    print.auc=TRUE, show.thres=TRUE)
dev.off()


pdf(file = paste0("Protein.Loss.ROC.AS.", name, ".pdf"),
    width = 3, 
    height = 3)
roc(All.Factors.Diff$Three.Protein.Loss=="Anti-Scaling", TestDataCol,
    smoothed = TRUE,
    plot=TRUE, auc.polygon=FALSE, max.auc.polygon=FALSE, grid=FALSE,
    print.auc=TRUE, show.thres=TRUE)
dev.off()

# Scaling
pdf(file = paste0("Protein.Gain.ROC.Scaling.", name, ".pdf"),
    width = 3, 
    height = 3)
roc(All.Factors.Diff$Three.Protein.Gain=="Scaling", TestDataCol, 
    smoothed = TRUE,
    plot=TRUE, auc.polygon=FALSE, max.auc.polygon=FALSE, grid=FALSE,
    print.auc=TRUE, show.thres=TRUE)
dev.off()


pdf(file = paste0("Protein.Loss.ROC.Scaling.", name, ".pdf"),
    width = 3, 
    height = 3)
roc(All.Factors.Diff$Three.Protein.Loss=="Scaling", TestDataCol,
    smoothed = TRUE,
    plot=TRUE, auc.polygon=FALSE, max.auc.polygon=FALSE, grid=FALSE,
    print.auc=TRUE, show.thres=TRUE)
dev.off()
