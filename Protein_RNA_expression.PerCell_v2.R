###### Difference upon gain/loss dataset(s), and RNA and Protein expression vs DNA copy number boxplots ####
### 210127
### by Klaske Schukken
### Graphs for Depmap Protein & RNA data: 
### Graph: protein/RNA expression per cell in mono/disomy/trisomy categories
### Edited: 210204-- only in cell lines with both RNA and Protein data. 
##                -- and look at chromosome arm level gains/losses.

library('ggplot2')
library('tidyverse')
library('xlsx')
library(readxl)
library(reshape2)
library('BBmisc')
library('tidyr')
library(tidyverse) 

###### Step 1: Get data ####
### Step 1: Get data. 
## Get raw Protein and RNA Data. 

## Import Protein expression data and uri ben-david lab, Cohen-Sharir nature 2021 chromosome arm data
##

DataFileLocation= "/Documents" ## ! Set file path to place you want files to be writen to.

##New Protein data. Filtered for only cells with both RNA and Protein data: 
## generated in Protein_RNA_filtered_CellLine.R
setwd()
Protein.Expression.filtered<-read.delim2("Protein_Expression_filtered.csv", 
                                         dec=".", header = TRUE, sep=",")

#Name dataframe with protein_ID and corresponding Protein name, this way we have both gene ID and Protein name
Protein_ProID<-Protein_Expression.1[,c(1:2, 6)] 
Protein_ProID$Protein_Id<- str_replace_all(Protein_ProID$Protein_Id, "[|]", ".")
Protein_ProID$Protein_Id<- str_replace_all(Protein_ProID$Protein_Id, "[-]", ".")


#Protein Info. names, chrm location, etc. 
#setwd()
Protein_Info<-read_excel("HGNC_Protein_Info_edited.xlsx", sheet ="HGNC names")

#Add chrm arm data
Protein_Info4<-Protein_Info
Protein_Info4$arm <- gsub('[0-9]+', '', Protein_Info4$'Chromosome band') #Find test protein chromosome arm
Protein_Info4$arm <- gsub('[.]', '', Protein_Info4$arm) #also remove period, if needed
Protein_Info4$arm <- str_sub(Protein_Info4$arm, -1, -1) #start & end on last character. get only p or q
#filter cell data and get only aneuploidy scores for chrm arm of test protein location

###
#setwd()
aneuploid<- read.csv("arm_calls_data_bendavid.csv", header=TRUE)

####
# Import RNA expression data (no aneuploidy data). CCLE depmap data.

## RNA data. Filtered for only cells with both RNA and Protein data: 
#  RNA data for only cells that have both RNA and Protein Data available. (371 cell lines)
## a few extra cells are filtered out when we add aneuploidy data--> 367 cells total
## RNA expression differences generated in Protein_RNA_filtered_CellLine.R
#setwd()
RNA.Expression.filtered<-read.delim2("RNA_Expression_filtered.csv", 
                                     dec=".", header = TRUE, sep=",")

#RNA Info. names, chrm location, etc. 
#setwd()
RNA_Info<-read_excel("HGNC_Protein_Info_edited.xlsx", sheet ="HGNC names")
RNA_Info2<-RNA_Info[,c(2,3,11,12,22,14)]

RNA_Info3<-RNA_Info2
RNA_Info3$arm <- gsub('[0-9]+', '', RNA_Info3$'Chromosome band') #Find test RNA chromosome arm
RNA_Info3$arm <- gsub('[.]', '', RNA_Info3$arm) #also remove period, if needed
RNA_Info3$arm <- str_sub(RNA_Info3$arm, -1, -1) #start & end on last character. get only p or q


setwd(DataFileLocation) # set wd to where you want the graphs to go to.


###### Find RNA expression difference & t-test list, all ####
## Step 5: make list of all RNAs 
# MINIMUM OF 10 CELLS per cateogry (gain/neutral/loss), per gene. 
#colnames(RNA.Expression.filtered)<-sub("[..].*", "", as.character(colnames(RNA.Expression.filtered)))#12755 genes
#RNA.Expression.filtered<-RNA.Expression.filtered[,-c(1)]

t.test_RNA_category<- data.frame(RNA_ID=character(),
                                 RNA_Name=character(),
                                 Pvalue.Gain= numeric(),
                                 Diff.Gain= numeric(), 
                                 Pvalue.Loss= numeric(),
                                 Diff.Loss= numeric(),
                                 Pvalue.GainvLoss= numeric(),
                                 Diff.GainvLoss= numeric())

for (i in 2:length(RNA.Expression.filtered)){
  
  TestRNA_name=sub("[..].*", "", as.character(colnames(RNA.Expression.filtered[i])))
  TestRNA=colnames(RNA.Expression.filtered[i])
  
  testRNAChrm <- filter(RNA_Info3, RNA_Info3$"Approved Symbol"==TestRNA_name) #get test RNA data
  if (length(testRNAChrm$Chromosome)!=0){ #Make sure RNA data is in RNA_info3, else it crashes
    TestArm <- testRNAChrm$arm #Find test RNA chromosome arm
    TestChrm<- as.numeric(testRNAChrm$Chromosome) #Find test RNA chromosome
    #filter cell data and get only aneuploidy scores for chrm arm of test RNA location
    testchrm.percell <- filter(aneuploid, chrom==TestChrm) 
    testchrm.percell <- filter(testchrm.percell, arm==TestArm)
    
    RNA_Expression_TestRNA <- RNA.Expression.filtered %>% select(Cell_line, all_of(TestRNA))
    
    #Combine RNA data with aneuploidy data for test RNA location
    RNA.effect.Chrm<-merge(y= RNA_Expression_TestRNA, x= testchrm.percell, 
                           by.y="Cell_line", by.x="DepMap_ID", 
                           sort = TRUE)# 368 cells, 12757 RNAs 
    
    ## put data into categories based on chrm arm number
    Chrm.tri<- filter(RNA.effect.Chrm, arm_call==1) #all cells with chrm gain for  testRNA chrm
    Chrm.di<- filter(RNA.effect.Chrm, arm_call==0) #all cells with chrm neutral ploidy for  testRNA chrm
    Chrm.mono<- filter(RNA.effect.Chrm, arm_call==-1) #all cells with chrm loss for  testRNA chrm
    
    ##Make data frame with t-test info about RNA expression per arm number category.
    ## Return this data at the end of the function. 
    if (colSums(!is.na(Chrm.tri[8]))>=10 &  #minimum 10 cells in "tri" category
        #Added if statement so only genes with 2+ values per condition are analyzed
        #if I don't do this, the t-test crashes and I get no values. 
        colSums(!is.na(Chrm.di[8]))>=10 &
        colSums(!is.na(Chrm.mono[8]))>=10) {
      # Trisomy vs disomy
      di.tri<-t.test(Chrm.tri[,8], # [,8] because that is column with RNA data
                     Chrm.di[,8], # [,8] because that is column with RNA data
                     mu = 0, 
                     alt = "two.sided",
                     conf.level = 0.99) #get p-value
      Diff.Tri.Di<-mean(Chrm.tri[,8], na.rm=TRUE) - mean(Chrm.di[,8], na.rm=TRUE) #get difference
      #Disomy vs monosome
      Mono.Di<-t.test(Chrm.mono[,8], # [,8] because that is column with RNA data
                      Chrm.di[,8], # [,8] because that is column with RNA data
                      mu = 0, 
                      alt = "two.sided",
                      conf.level = 0.99) #get p-value
      Diff.Mono.Di<- mean(Chrm.mono[,8], na.rm=TRUE) - mean(Chrm.di[,8], na.rm=TRUE) #get difference
      # Tri vs monosomy
      Tri.Mono<-t.test(Chrm.mono[,8], # [,8] because that is column with RNA data
                       Chrm.tri[,8], # [,8] because that is column with RNA data
                       mu = 0, 
                       alt = "two.sided",
                       conf.level = 0.99) #get p-value
      Diff.Tri.Mono<-mean(Chrm.tri[,8], na.rm=TRUE) - mean(Chrm.mono[,8], na.rm=TRUE) #get difference
      
      t.test_RNA_category<- rbind(t.test_RNA_category, data.frame(
        RNA_ID=colnames(RNA.effect.Chrm[8]),
        RNA_Name=TestRNA_name,
        Pvalue.Gain= di.tri$p.value,
        Diff.Gain= Diff.Tri.Di, 
        Pvalue.Loss= Mono.Di$p.value,
        Diff.Loss= Diff.Mono.Di, 
        Pvalue.GainvLoss= Tri.Mono$p.value,
        Diff.GainvLoss= Diff.Tri.Mono))
    }
  }
}
t.test_RNA_category<-distinct(t.test_RNA_category)#17774 genes
t.test_RNA_category<-t.test_RNA_category[order(t.test_RNA_category$Pvlaue.Tri.Mono),]

#17621 RNAs measured. 
setwd(DataFileLocation)
write.csv(t.test_RNA_category, 
          file =paste("RNA_Loss.Neutral.Gain_Difference_Pvalue_min10points.csv", sep=','), 
          row.names = TRUE)

t.test_RNA_category<-read.delim2("RNA_Loss.Neutral.Gain_Difference_Pvalue_min10points.csv", 
                                 dec=".", header = TRUE, sep=",")

###### Find Protein expression difference & t-test list, all ####
## make list of pvalue and difference between mono- di-tri cells per gene 
## for each gene in Protein data. using only cells that have both RNA and Protein data

##First step is to combine Protein Uniprot ID info with Depmap Protein_ID
## Using information given by initial mmcr protein database. see above. 
## Get Protein_ID as dataframe, combine with uniprot id and Protein symbol data. 
Protein.ID<-data.frame(Protein_ID=
                         colnames(Protein.Expression.filtered[3:length(Protein.Expression.filtered)]))

##Now merge depmap column name info with Protein info. 
## First add Uniprot ID and gene symbol info to depmap data. 
Depmap.Protein.info2<-merge(x= Protein.ID, y= Protein_ProID, 
                            by.x="Protein_ID", by.y="Protein_Id", 
                            sort = TRUE)# 3 collumns, 12755. all Proteins given Uniprot IDs.
# Now combine all genes with same uniprot ID: 
Depmap.Protein.info3<-merge(x= Depmap.Protein.info2, y= Protein_Info4, 
                            by.x="Uniprot_Acc", by.y="UniProt ID(supplied by UniProt)", 
                            sort = TRUE)# 10594 Proteins
# Find those genes without uniprot id: 
No_UniprotID <- anti_join(Depmap.Protein.info2, Protein_Info4, #finding genes_Symbol without match
                          by = c("Uniprot_Acc" = "UniProt ID(supplied by UniProt)"))
#...find the genes with matching gene symbols
Depmap.Protein.info4<-merge(x= No_UniprotID, y= Protein_Info4, 
                            by.x="Gene_Symbol", by.y="Approved Symbol", 
                            sort = TRUE)# 10 collumns, 253 genes, only no gene-symbol genes with uniprot_ID

# Merge genes with gene symbol and genes with only uniprot ID. 
Depmap.Protein.info5<-merge(x= Depmap.Protein.info3, y= Depmap.Protein.info4, 
                            all=TRUE)# 11 collumns, 12100 genes


t.test_prot_category<- data.frame(Protein_ID=character(),
                                  Protein_Name=character(),
                                  Pvlaue.Tri.Di= numeric(),
                                  Diff.Gain= numeric(), 
                                  Pvlaue.Di.Mono= numeric(),
                                  Diff.Di_Mono= numeric(),
                                  Pvlaue.Tri.Mono= numeric(),
                                  Diff.GainvLoss= numeric())



for (i in 3:length(Protein.Expression.filtered)){
  TestProt = colnames(Protein.Expression.filtered[i])
  
  testProtChrm <- filter(Depmap.Protein.info5, Protein_ID==TestProt) #get test protein data
  
  
  #use If statement to check that Protein Info has Protein of interest. 
  if (length(testProtChrm[,1])!=0){
    TestArm <- testProtChrm$arm #Find test protein chromosome arm
    TestChrm<- as.numeric(testProtChrm$Chromosome) #Find test protein chromosome
    #filter cell data and get only aneuploidy scores for chrm arm of test protein location
    testchrm.percell <- filter(aneuploid, chrom==TestChrm) 
    testchrm.percell <- filter(testchrm.percell, arm==TestArm)
    
    Protein_Expression_TestProt <- Protein.Expression.filtered %>% select(Cell_line, all_of(TestProt))
    
    #Combine protein data with aneuploidy data for test protein location
    Protein.effect.Chrm<-merge(y= Protein_Expression_TestProt, x= testchrm.percell, 
                               by.y="Cell_line", by.x="DepMap_ID", 
                               sort = TRUE)# 368 cells, 12757 proteins 
    
    ## put data into categories based on chrm arm number
    Chrm.tri<- filter(Protein.effect.Chrm, arm_call==1) #all cells triploid for testProt chrm
    Chrm.di<- filter(Protein.effect.Chrm, arm_call==0) #all cells diploid for testProt chrm
    Chrm.mono<- filter(Protein.effect.Chrm, arm_call==-1) #all cells mono for testProt chrm
    
    ##Make data frame with t-test info about protein expression per arm number category.
    ## Return this data at the end of the function. 
    if (colSums(!is.na(Chrm.tri[8]))>=10 & #[8] is column with protein expression data. check it has values. 
        #Added if statement so only genes with 2+ values per condition are analyzed
        #if I don't do this, the t-test crashes and I get no values. 
        colSums(!is.na(Chrm.di[8]))>=10 &
        colSums(!is.na(Chrm.mono[8]))>=10) {
      # Trisomy vs disomy
      di.tri<-t.test(Chrm.tri[,8], # [,8] because that is column with Protein data
                     Chrm.di[,8], # [,8] because that is column with Protein data
                     mu = 0, 
                     alt = "two.sided",
                     conf.level = 0.99) #get p-value
      Diff.Tri.Di<-mean(Chrm.tri[,8], na.rm=TRUE) - mean(Chrm.di[,8], na.rm=TRUE) #get difference
      #Disomy vs monosome
      Di.Mono<-t.test(Chrm.mono[,8], # [,8] because that is column with Protein data
                      Chrm.di[,8], # [,8] because that is column with Protein data
                      mu = 0, 
                      alt = "two.sided",
                      conf.level = 0.99) #get p-value
      Diff.Di.Mono<-mean(Chrm.di[,8], na.rm=TRUE) - mean(Chrm.mono[,8], na.rm=TRUE) #get difference
      Di.Mono.ploidy<- mean(Chrm.di[,6], na.rm=TRUE)
      # Tri vs monosomy
      Tri.Mono<-t.test(Chrm.mono[,8], # [,8] because that is column with Protein data
                       Chrm.tri[,8], # [,8] because that is column with Protein data
                       mu = 0, 
                       alt = "two.sided",
                       conf.level = 0.99) #get p-value
      Diff.Tri.Mono<-mean(Chrm.tri[,8], na.rm=TRUE) - mean(Chrm.mono[,8], na.rm=TRUE) #get difference
      Tri.Mono.ploidy<- mean(Chrm.tri[,6], na.rm=TRUE)
      
      t.test_prot_category<- rbind(t.test_prot_category, 
                                   data.frame(Protein_ID=TestProt,
                                              Protein_Name=testProtChrm$Gene_Symbol,
                                              Pvlaue.Tri.Di= di.tri$p.value,
                                              Diff.Gain= Diff.Tri.Di, 
                                              Pvlaue.Di.Mono= Di.Mono$p.value,
                                              Diff.Di_Mono= Diff.Di.Mono, 
                                              Pvlaue.Tri.Mono= Tri.Mono$p.value,
                                              Diff.GainvLoss= Diff.Tri.Mono))
    } else {
      t.test_prot_category<-t.test_prot_category
    }
  }
}


t.test_prot_category<-t.test_prot_category[order(t.test_prot_category$Pvlaue.Tri.Mono),]#11458 genes. only the proteins with 3+ cells per mono/di.tri category. lost ~1k proteins. 
t.test_prot_category<-distinct(t.test_prot_category)
t.test_prot_category$Diff.Di_Mono<-t.test_prot_category$Diff.Di_Mono*-1 #make this into Aneu-Diploid.

setwd(DataFileLocation)
write.csv(t.test_prot_category, 
          file =paste("Protein_Loss.Neutral.Gain_Difference_Pvalue_min10points.csv", sep=','), 
          row.names = TRUE)

#t.test_prot_category<-read.delim2("Protein_Loss.Neutral.Gain_Difference_Pvalue_min10points.csv", 
#                                  dec=".", header = TRUE, sep=",")


### Combine RNA and Protein Data: 
## This is with filtered cell lines (only cells with RNA and Protein and aneuploid data:)
## This is with filtered genes: only genes with 10+ cells in all categories: RNA & protein gain, neutral and loss
CN.Diff.xRNA.yProt<-merge(x=t.test_RNA_category, y=t.test_prot_category, by.x="RNA_Name", by.y="Protein_Name")
#CN.Diff.xRNA.yProt<- subset(CN.Diff.xRNA.yProt, select=-c(X,X.x, X.y, X.1,X.2))

colnames(CN.Diff.xRNA.yProt)<-c("RNA_Name", "RNA_ID", 
                                "RNA.Pvalue.Gain", "RNA.Diff.Gain", 
                                "RNA.Pvalue.Loss", "RNA.Diff.Loss", 
                                "RNA.Pvalue.GainvLoss", "RNA.Diff.GainvLoss", 
                                "Protein_ID", "Protein.Pvalue.Gain", 
                                "Protein.Diff.Gain", "Protein.Pvalue.Loss", 
                                "Protein.Diff.Loss", "Protein.Pvalue.GainvLoss", 
                                "Protein.Diff.GainvLoss")
#CN.Diff.xRNA.yProt<-CN.Diff.xRNA.yProt[,-c(2,10:11)]

write.csv(CN.Diff.xRNA.yProt, 
          file =paste("RNA.Protein_loss.neutral.gain_Difference_Pvalue_min10points.csv", sep=','), 
          row.names = TRUE)

#read.csv(CN.Diff.xRNA.yProt, 
#          file =paste("RNA.Protein_loss.neutral.gain_Difference_Pvalue_min10points.csv", sep=','), 
#          row.names = TRUE)

write.csv(CN.Diff.xRNA.yProt.ThreeGroups, 
          file =paste("RNA.Protein_loss.neutral.gain_Difference_Pvalue_min10points_3cat.csv", sep=','), 
          row.names = TRUE)

###### Categorize genes by significance and difference #####
## Categorize genes by if they scale or buffer upon chrm gain/loss. plot. 
## Did not end up using these graphs in the paper. but can use this to find highly 
##    significant scaling/buffering proteins/RNAs. 

## Find genes per category: Per chrm gain: 
## Categories based on Significance and +/- 0 Difference,
## can be used to find genes that are significantly Scaling/Buffering (or NS) for gain & loss

CN.Diff.xRNA.yProt<-merge(x=t.test_RNA_category, y=t.test_prot_category, by.x="RNA_ID", by.y="Protein_Name")

CN.gain.RNAbuff.Protbuff<- subset(CN.Diff.xRNA.yProt, 
                                  (RNA.Pvalue.Gain>0.05 | RNA.Diff.Gain<0) & 
                                    (Protein.Pvalue.Gain>0.05 | Protein.Diff.Gain<0) )
CN.gain.RNAbuff.Protbuff<- CN.gain.RNAbuff.Protbuff[order(CN.gain.RNAbuff.Protbuff$Protein.Pvalue.Gain),]
CN.gain.RNAscale.Protscale<- subset(CN.Diff.xRNA.yProt, 
                                    (RNA.Pvalue.Gain<0.05 & RNA.Diff.Gain>0) & 
                                      (Protein.Pvalue.Gain<0.05 & Protein.Diff.Gain>0) )
CN.gain.RNAscale.Protscale<- CN.gain.RNAscale.Protscale[order(CN.gain.RNAscale.Protscale$Protein.Pvalue.Gain),]
CN.gain.RNAscale.Protbuff<- subset(CN.Diff.xRNA.yProt, 
                                   (RNA.Pvalue.Gain<0.05 & RNA.Diff.Gain>0) & 
                                     (Protein.Pvalue.Gain>0.05 | Protein.Diff.Gain<0) )
CN.gain.RNAscale.Protbuff<- CN.gain.RNAscale.Protbuff[order(CN.gain.RNAscale.Protbuff$Protein.Pvalue.Gain),]
CN.gain.RNAbuff.Protscale<- subset(CN.Diff.xRNA.yProt, 
                                   (RNA.Pvalue.Gain>0.05 | RNA.Diff.Gain<0) &
                                     (Protein.Pvalue.Gain<0.05 & Protein.Diff.Gain>0) )
CN.gain.RNAbuff.Protscale<- CN.gain.RNAbuff.Protscale[order(CN.gain.RNAbuff.Protscale$Protein.Pvalue.Gain),]

length(CN.gain.RNAbuff.Protbuff$RNA_ID) 
CN.gain.RNAbuff.Protbuff$RNA_Name[1:5]
length(CN.gain.RNAscale.Protscale$RNA_ID) 
CN.gain.RNAscale.Protscale$RNA_Name[1:5]
length(CN.gain.RNAscale.Protbuff$RNA_ID) 
CN.gain.RNAscale.Protbuff$RNA_Name[1:5]
length(CN.gain.RNAbuff.Protscale$RNA_ID) 
CN.gain.RNAbuff.Protscale$RNA_Name[1:5]

## Upon chromosome gain: 
##Minimum of 3 data points per condition (-1, 0, +1):
#RNA buff, Protein Buffered: 4737 genes 
#     top 5: IRF1    DTNA    CYBRD1  ZMYND19 FCGRT  
#RNA scale, Protein scales: 2309 genes 
#     top 5 most significant: PURB   NUDCD3 STAU1  RPRD1B CDK5
#RNA scale, Protein buff: 3450 genes 
#     top 5 most significant: ACTR5   MLKL    C9orf85 TMEM14C MTCH1  
#RNA buff, Protein scales: 445 genes 
#     top 5 most significant: IP6K2  MUM1   ABHD10 DOK3   RNF223

## Upon chromosome gain: 
##Minimum of 10 data points per condition (-1, 0, +1):
#RNA buff, Protein Buffered: 3849 genes 
#     top 5: LAMC2 AP1M2 DOK3  PCDH9 LAMB3  
#RNA scale, Protein scales: 2205 genes 
#     top 5 most significant: PURB   NUDCD3 STAU1  RPRD1B CDK5
#RNA scale, Protein buff: 3013 genes 
#     top 5 most significant: TNNC2    TRAPPC12 RFC2     HAUS4    PSMB2     
#RNA buff, Protein scales: 346 genes 
#     top 5 most significant: ABHD10   TMEM126A RABEP2   ALDOA    PGAM5


####
##Find categories for chromosome loss: 
CN.loss.RNAbuff.Protbuff<- subset(CN.Diff.xRNA.yProt, 
                                  (RNA.Pvalue.Loss>0.05 | RNA.Diff.Loss>0) & 
                                    (Protein.Pvalue.Loss>0.05 | Protein.Diff.Loss<0) )
CN.loss.RNAbuff.Protbuff<- CN.loss.RNAbuff.Protbuff[order(CN.loss.RNAbuff.Protbuff$Protein.Pvalue.Loss),]

CN.loss.RNAscale.Protscale<- subset(CN.Diff.xRNA.yProt, 
                                    (RNA.Pvalue.Loss<0.05 & RNA.Diff.Loss<0) & 
                                      (Protein.Pvalue.Loss<0.05 & Protein.Diff.Loss<0) )
CN.loss.RNAscale.Protscale<- CN.loss.RNAscale.Protscale[order(CN.loss.RNAscale.Protscale$Protein.Pvalue.Loss),]

CN.loss.RNAscale.Protbuff<- subset(CN.Diff.xRNA.yProt, 
                                   (RNA.Pvalue.Loss<0.05 & RNA.Diff.Loss<0) & 
                                     (Protein.Pvalue.Loss>0.05 | Protein.Diff.Loss>0) )
CN.loss.RNAscale.Protbuff<- CN.loss.RNAscale.Protbuff[order(CN.loss.RNAscale.Protbuff$Protein.Pvalue.Loss),]

CN.loss.RNAbuff.Protscale<- subset(CN.Diff.xRNA.yProt, 
                                   (RNA.Pvalue.Loss>0.05 | RNA.Diff.Loss>0) & 
                                     (Protein.Pvalue.Loss<0.05 & Protein.Diff.Loss<0) )
CN.loss.RNAbuff.Protscale<- CN.loss.RNAbuff.Protscale[order(CN.loss.RNAbuff.Protscale$Protein.Pvalue.Loss),]

length(CN.loss.RNAbuff.Protbuff$RNA_ID) 
CN.loss.RNAbuff.Protbuff$RNA_Name[1:5]
length(CN.loss.RNAscale.Protscale$RNA_ID) 
CN.loss.RNAscale.Protscale$RNA_Name[1:5]
length(CN.loss.RNAscale.Protbuff$RNA_ID) 
CN.loss.RNAscale.Protbuff$RNA_Name[1:5]
length(CN.loss.RNAbuff.Protscale$RNA_ID) 
CN.loss.RNAbuff.Protscale$RNA_Name[1:5]

## Upon chromosome loss: 
##Minimum of 3 data points per condition (-1, 0, +1):
#RNA buff, Protein Buffered: 4059 genes 
#     top 5: GLI2   EDIL3  GMPS   TXNRD3 CEBPB
#RNA scale, Protein scales: 3021 genes 
#     top 5 most significant: PDE12  TXNL1  HDHD2  NARS   PPP2CB
#RNA scale, Protein buff: 3503 genes 
#     top 5 most significant: RASSF2 ZNF622 ADD3   ATMIN  NHSL1 
#RNA buff, Protein scales: 358 genes 
#     top 5 most significant: DNAJB1   MEX3D    UBE2E2   CHMP2B   ARHGEF18

## Upon chromosome loss: 
##Minimum of 3 data points per condition (-1, 0, +1):
#RNA buff, Protein Buffered: 3386 genes 
#     top 5: DSG2     DSG3     CYBRD1   SERPINB5 S100A10 
#RNA scale, Protein scales: 2855 genes 
#     top 5 most significant: PDE12  TXNL1  HDHD2  NARS   PPP2CB
#RNA scale, Protein buff: 2911 genes 
#     top 5 most significant: CTNNB1 TAF13  TRIM11 STARD9 SPCS1  
#RNA buff, Protein scales: 261 genes 
#     top 5 most significant: DNAJB1   MEX3D    UBE2E2   CHMP2B   ARHGEF18

###Now to look at genes that are buffeted/scaling in BOTH gain and loss 
CN.aneu.RNAbuff.Protbuff<- subset(CN.Diff.xRNA.yProt, 
                                  ((RNA.Pvalue.Loss>0.05 | RNA.Diff.Loss>0) & (RNA.Pvalue.Gain>0.05 |RNA.Diff.Gain<0) & (RNA.Pvalue.GainvLoss>0.05| RNA.Diff.GainvLoss<0)) & # RNA not sig OR wrong direction
                                    ((Protein.Pvalue.Loss>0.05 | Protein.Diff.Loss>0) & (Protein.Pvalue.Gain>0.05 |Protein.Diff.Gain<0) & (Protein.Pvalue.GainvLoss>0.05| Protein.Diff.GainvLoss<0)))  # Protein not sig OR wrong direction
CN.aneu.RNAbuff.Protbuff<- CN.aneu.RNAbuff.Protbuff[order(CN.aneu.RNAbuff.Protbuff$Protein.Pvalue.GainvLoss),]

CN.aneu.RNAscale.Protscale<- subset(CN.Diff.xRNA.yProt, 
                                    ((RNA.Pvalue.Loss<0.05 & RNA.Diff.Loss<0) & (RNA.Pvalue.Gain<0.05 & RNA.Diff.Gain>0) & (RNA.Pvalue.GainvLoss<0.05 & RNA.Diff.GainvLoss>0)) & # RNA sig AND scales directionally
                                      ((Protein.Pvalue.Loss<0.05 & Protein.Diff.Loss<0) & (Protein.Pvalue.Gain<0.05 & Protein.Diff.Gain>0) & (Protein.Pvalue.GainvLoss<0.05 & Protein.Diff.GainvLoss>0)))# Protein sig AND scales directionally
CN.aneu.RNAscale.Protscale<- CN.aneu.RNAscale.Protscale[order(CN.aneu.RNAscale.Protscale$Protein.Pvalue.GainvLoss),]

CN.aneu.RNAscale.Protbuff<- subset(CN.Diff.xRNA.yProt, 
                                   ((RNA.Pvalue.Loss<0.05 & RNA.Diff.Loss<0) & (RNA.Pvalue.Gain<0.05 & RNA.Diff.Gain>0) & (RNA.Pvalue.GainvLoss<0.05 & RNA.Diff.GainvLoss>0)) & # RNA sig AND scales directionally
                                     ((Protein.Pvalue.Loss>0.05 | Protein.Diff.Loss>0) & (Protein.Pvalue.Gain>0.05 |Protein.Diff.Gain<0) & (Protein.Pvalue.GainvLoss>0.05| Protein.Diff.GainvLoss<0)))  # Protein not sig OR wrong direction
CN.aneu.RNAscale.Protbuff<- CN.aneu.RNAscale.Protbuff[order(CN.aneu.RNAscale.Protbuff$Protein.Pvalue.GainvLoss),]

CN.aneu.RNAbuff.Protscale<- subset(CN.Diff.xRNA.yProt, 
                                   ((RNA.Pvalue.Loss>0.05 | RNA.Diff.Loss>0) & (RNA.Pvalue.Gain>0.05 |RNA.Diff.Gain<0) & (RNA.Pvalue.GainvLoss>0.05| RNA.Diff.GainvLoss<0)) & # RNA not sig OR wrong direction
                                     ((Protein.Pvalue.Loss<0.05 & Protein.Diff.Loss<0) & (Protein.Pvalue.Gain<0.05 & Protein.Diff.Gain>0) & (Protein.Pvalue.GainvLoss<0.05 & Protein.Diff.GainvLoss>0)))# Protein sig AND scales directionally
CN.aneu.RNAbuff.Protscale<- CN.aneu.RNAbuff.Protscale[order(CN.aneu.RNAbuff.Protscale$Protein.Pvalue.GainvLoss),]

length(CN.aneu.RNAbuff.Protbuff$RNA_ID) 
CN.aneu.RNAbuff.Protbuff$RNA_Name[1:5]
length(CN.aneu.RNAscale.Protscale$RNA_ID) 
CN.aneu.RNAscale.Protscale$RNA_Name[1:5]
length(CN.aneu.RNAscale.Protbuff$RNA_ID) 
CN.aneu.RNAscale.Protbuff$RNA_Name[1:5]
length(CN.aneu.RNAbuff.Protscale$RNA_ID) 
CN.aneu.RNAbuff.Protscale$RNA_Name[1:5]

## Buffered upon BOTH gain and loss, scale with either one or other: 
##Minimum of 3 data points per condition (-1, 0, +1):
#RNA buff, Protein Buffered: 1489 genes 
#     top 5: IRF1   IAH1   NOS1AP RAB21  HPCAL4
#RNA scale, Protein scales: 884 genes 
#     top 5 most significant: PTPN2  SMCHD1 RNMT   USP14  PPP4R1
#RNA scale, Protein buff: 1425 genes 
#     top 5 most significant: DICER1  CARD8   RAPGEF1 DCLRE1C MDM4 
#RNA buff, Protein scales: 16 genes 
#     top 5 most significant: SIX5     STMN4    GOLGA3   GLB1L3   ARHGEF28

## Buffered upon BOTH gain and loss, scale with either one or other: 
##Minimum of 10 data points per condition (-1, 0, +1):
#RNA buff, Protein Buffered: 1705 genes 
#     top 5: DSG3     AP1M2    DSC2     SERPINB5 S100A10
#RNA scale, Protein scales: 860 genes 
#     top 5 most significant: PTPN2  SMCHD1 RNMT   USP14  PPP4R1
#RNA scale, Protein buff: 1014 genes 
#     top 5 most significant: RAB11B COPG1  HDAC3  TRIM11 HAUS4 
#RNA buff, Protein scales: 3 genes 
#     top 5 most significant: GLB1L3  ST3GAL2 SLC27A2 

ProtExp.ChrmCN.filtered("GLB1L3")
RNAExp.ChrmCN.filtered("GLB1L3")
ProtExp.ChrmCN.filtered("ST3GAL2")
RNAExp.ChrmCN.filtered("ST3GAL2")
ProtExp.ChrmCN.filtered("SLC27A2")
RNAExp.ChrmCN.filtered("SLC27A2")


###Some other genes that I had for some reason
ProtExp.ChrmCN.filtered("DSG3")
RNAExp.ChrmCN.filtered("DSG3")
ProtExp.ChrmCN.filtered("ABL2")
RNAExp.ChrmCN.filtered("ABL2")
ProtExp.ChrmCN.filtered("STMN4")
RNAExp.ChrmCN.filtered("STMN4")
ProtExp.ChrmCN.filtered("TROAP")
RNAExp.ChrmCN.filtered("TROAP")
ProtExp.ChrmCN.filtered("ERBB4")
RNAExp.ChrmCN.filtered("ERBB4")


#Make category vectors for buffer/scale upon gain, loss or both: 

CN.Diff.xRNA.yProt$Gain<-NA 
for (i in 1:length(CN.Diff.xRNA.yProt$RNA_ID)){
  if (CN.Diff.xRNA.yProt$RNA_ID[i] %in% CN.gain.RNAbuff.Protbuff$RNA_ID){
    CN.Diff.xRNA.yProt$Gain[i]<-"RNA buffered, Protein buffered"
  } else if (CN.Diff.xRNA.yProt$RNA_ID[i] %in% CN.gain.RNAscale.Protscale$RNA_ID){
    CN.Diff.xRNA.yProt$Gain[i]<-"RNA scales, Protein scales"
  } else if (CN.Diff.xRNA.yProt$RNA_ID[i] %in% CN.gain.RNAscale.Protbuff$RNA_ID){
    CN.Diff.xRNA.yProt$Gain[i]<-"RNA scales, Protein buffered"
  } else if (CN.Diff.xRNA.yProt$RNA_ID[i] %in% CN.gain.RNAbuff.Protscale$RNA_ID){
    CN.Diff.xRNA.yProt$Gain[i]<-"RNA buffered, Protein scales"
  }
}
CN.Diff.xRNA.yProt$Gain<-factor(CN.Diff.xRNA.yProt$Gain, levels=c("RNA buffered, Protein buffered", "RNA scales, Protein buffered", "RNA scales, Protein scales", "RNA buffered, Protein scales"))

CN.Diff.xRNA.yProt$Loss<-NA 
for (i in 1:length(CN.Diff.xRNA.yProt$RNA_ID)){
  if (CN.Diff.xRNA.yProt$RNA_ID[i] %in% CN.loss.RNAbuff.Protbuff$RNA_ID){
    CN.Diff.xRNA.yProt$Loss[i]<-"RNA buffered, Protein buffered"
  } else if (CN.Diff.xRNA.yProt$RNA_ID[i] %in% CN.loss.RNAscale.Protscale$RNA_ID){
    CN.Diff.xRNA.yProt$Loss[i]<-"RNA scales, Protein scales"
  } else if (CN.Diff.xRNA.yProt$RNA_ID[i] %in% CN.loss.RNAscale.Protbuff$RNA_ID){
    CN.Diff.xRNA.yProt$Loss[i]<-"RNA scales, Protein buffered"
  } else if (CN.Diff.xRNA.yProt$RNA_ID[i] %in% CN.loss.RNAbuff.Protscale$RNA_ID){
    CN.Diff.xRNA.yProt$Loss[i]<-"RNA buffered, Protein scales"
  }
}
CN.Diff.xRNA.yProt$Loss<-factor(CN.Diff.xRNA.yProt$Loss, levels=c("RNA buffered, Protein buffered", "RNA scales, Protein buffered", "RNA scales, Protein scales", "RNA buffered, Protein scales"))

CN.Diff.xRNA.yProt$Gain.Loss<-NA 
for (i in 1:length(CN.Diff.xRNA.yProt$RNA_ID)){
  if (CN.Diff.xRNA.yProt$RNA_ID[i] %in% CN.aneu.RNAbuff.Protbuff$RNA_ID){
    CN.Diff.xRNA.yProt$Gain.Loss[i]<-"RNA buffered, Protein buffered"
  } else if (CN.Diff.xRNA.yProt$RNA_ID[i] %in% CN.aneu.RNAscale.Protscale$RNA_ID){
    CN.Diff.xRNA.yProt$Gain.Loss[i]<-"RNA scales, Protein scales"
  } else if (CN.Diff.xRNA.yProt$RNA_ID[i] %in% CN.aneu.RNAscale.Protbuff$RNA_ID){
    CN.Diff.xRNA.yProt$Gain.Loss[i]<-"RNA scales, Protein buffered"
  } else if (CN.Diff.xRNA.yProt$RNA_ID[i] %in% CN.aneu.RNAbuff.Protscale$RNA_ID){
    CN.Diff.xRNA.yProt$Gain.Loss[i]<-"RNA buffered, Protein scales"
  }
}

CN.Diff.xRNA.yProt$Gain.Loss<-factor(CN.Diff.xRNA.yProt$Gain.Loss, levels=c("RNA buffered, Protein buffered", "RNA scales, Protein buffered", "RNA scales, Protein scales", "RNA buffered, Protein scales"))

table.diff.scale<-table(CN.Diff.xRNA.yProt$Gain)
table.diff.scale<-rbind(table.diff.scale,table(CN.Diff.xRNA.yProt$Loss))
table.diff.scale<-rbind(table.diff.scale,table(CN.Diff.xRNA.yProt$Gain.Loss))
rownames(table.diff.scale)<-c("Gain", "Loss", "Gain and Loss")

melt.table.diff.scale<- melt(table.diff.scale)

ggplot(data=melt.table.diff.scale, aes(x=Var1, y=value))+
  geom_bar(stat="identity", position="dodge", aes(fill=Var2))+
  xlab("Gene scaling or buffering upon chromosome aneuploidy")+
  ylab("Gene count")+ 
  scale_fill_manual(values=c("grey20", "grey40", "grey60", "grey80"))+
  labs(fill=c("Gene regulation type (significance):"))+
  theme_classic()




###### Categorize genes: Scaling, buffering. anti-scaling ####
## Make THREE Categories: scaling, anti-scale and buffering, based on difference. 
## Three groups: Scaling, buffering and anti-scaling
## Three categories, upon gain: sclaing (>0.25), buffering (-0.1 to 0.25), anti-scaling (< -0.1)
## 0.25 because thats halfway between .5 fold change, expected from 1.5-fold increase. 
## -0.1 to account for non-significant downregulation upon gain
## vice versa for chrm loss. 
CN.Diff.xRNA.yProt.ThreeGroups<- CN.Diff.xRNA.yProt
CN.Diff.xRNA.yProt.ThreeGroups$Three.RNA.Gain<- cut(CN.Diff.xRNA.yProt.ThreeGroups$RNA.Diff.Gain,
                                                    breaks=c(-Inf,-0.1,0.25,Inf),
                                                    include.lowest=TRUE,
                                                    labels=c("Anti-Scaling","Buffering","Scaling"))
CN.Diff.xRNA.yProt.ThreeGroups$Three.RNA.Loss<- cut(CN.Diff.xRNA.yProt.ThreeGroups$RNA.Diff.Loss,
                                                    breaks=c(-Inf,-0.25,0.1,Inf),
                                                    include.lowest=TRUE,
                                                    labels=c("Scaling", "Buffering","Anti-Scaling"))
CN.Diff.xRNA.yProt.ThreeGroups$Three.Protein.Gain<- cut(CN.Diff.xRNA.yProt.ThreeGroups$Protein.Diff.Gain,
                                                        breaks=c(-Inf,-0.1,0.25,Inf),
                                                        include.lowest=TRUE,
                                                        labels=c("Anti-Scaling","Buffering","Scaling"))
CN.Diff.xRNA.yProt.ThreeGroups$Three.Protein.Loss<- cut(CN.Diff.xRNA.yProt.ThreeGroups$Protein.Diff.Loss,
                                                        breaks=c(-Inf,-0.25,0.1,Inf),
                                                        include.lowest=TRUE,
                                                        labels=c("Scaling", "Buffering","Anti-Scaling"))
CN.Diff.xRNA.yProt.ThreeGroups$Three.Protein.Loss<-factor(CN.Diff.xRNA.yProt.ThreeGroups$Three.Protein.Loss, 
                                                          levels=c("Anti-Scaling","Buffering","Scaling"))
CN.Diff.xRNA.yProt.ThreeGroups$Three.RNA.Loss<-factor(CN.Diff.xRNA.yProt.ThreeGroups$Three.RNA.Loss, 
                                                          levels=c("Anti-Scaling","Buffering","Scaling"))

### count percentages
sum(CN.Diff.xRNA.yProt.ThreeGroups$Three.Protein.Loss=="Anti-Scaling")/length(CN.Diff.xRNA.yProt.ThreeGroups$RNA_Name)
x<-subset(CN.Diff.xRNA.yProt.ThreeGroups, Three.RNA.Gain=="Buffering")
sum(x$Three.Protein.Gain=="Buffering")/length(x$RNA_Name)
x<-subset(CN.Diff.xRNA.yProt.ThreeGroups, Three.Protein.Loss=="Buffering")
sum(x$Three.Protein.Gain=="Buffering")/length(x$RNA_Name)


##Scatterplots with categories colored in: 
#RNA Chrm gain graphs
ggplot(CN.Diff.xRNA.yProt.ThreeGroups, aes(x=RNA.Diff.Gain, y=-log2(RNA.Pvalue.Gain), color=Three.RNA.Gain))+
  geom_point(size=2)+
  scale_color_manual(values=c("salmon2", "Gold2", "palegreen3"))+
  xlab("Difference in RNA expression")+
  ylab("log2(p-value)")+
  labs(color = "Quartiles:")+ #legend title
  theme_classic()+
  coord_cartesian(xlim=c(-2.7, 2.7), ylim=c(0,100))+
  ggtitle("RNA Gain scatterplot quartiles")

#RNA Loss graphs
ggplot(CN.Diff.xRNA.yProt.ThreeGroups, aes(x=RNA.Diff.Loss, y=-log2(RNA.Pvalue.Loss), color=Three.RNA.Loss))+
  geom_point(size=2)+
  scale_color_manual(values=c("salmon2", "Gold2", "palegreen3"))+
  xlab("Difference in RNA expression")+
  ylab("log2(p-value)")+
  labs(color = "Quartiles:")+ #legend title
  theme_classic()+
  coord_cartesian(xlim=c(-2.7, 2.7), ylim=c(0,100))+
  ggtitle("RNA Loss scatterplot quartiles")

#Protein Chrm gain graphs
ggplot(CN.Diff.xRNA.yProt.ThreeGroups, aes(x=Protein.Diff.Gain, y=-log2(Protein.Pvalue.Gain), color=Three.Protein.Gain))+
  geom_point(size=2)+
  scale_color_manual(values=c("salmon2", "Gold2", "palegreen3"))+
  xlab("Difference in protein expression")+
  ylab("log2(p-value)")+
  labs(color = "Quartiles:")+ #legend title
  theme_classic()+
  coord_cartesian(xlim=c(-2.7, 2.7), ylim=c(0,100))+
  ggtitle("Protein Gain scatterplot quartiles")

#Protein Loss graphs
ggplot(CN.Diff.xRNA.yProt.ThreeGroups, aes(x=Protein.Diff.Loss, y=-log2(Protein.Pvalue.Loss), color=Three.Protein.Loss))+
  geom_point(size=2)+
  scale_color_manual(values=c("salmon2", "Gold2", "palegreen3"))+
  xlab("Difference in protein expression")+
  ylab("log2(p-value)")+
  labs(color = "Quartiles:")+ #legend title
  theme_classic()+
  coord_cartesian(xlim=c(-2.7, 2.7), ylim=c(0,100))+
  ggtitle("protein Loss scatterplot quartiles")


###Bar graph of RNA and Protein quantiles Scaling/Buffered

#format data so I can plot all 4 groups: RNA/Prot gain/loss
dat.m <- melt(CN.Diff.xRNA.yProt.ThreeGroups,id.vars='RNA_ID', measure.vars=c('Three.RNA.Gain','Three.Protein.Gain','Three.RNA.Loss', 'Three.Protein.Loss'))

## RNA & Protein Gain
ggplot(data= dat.m, aes(x=variable, fill=value)) + 
  geom_bar(position="fill")+ #dodge (next to eachother)
  scale_fill_manual(values=c("salmon2", "Gold2", "palegreen3"))+
  xlab("RNA Difference upon chrm arm gain: category")+
  ylab("Percent of genes")+
  labs(fill = "Protein Difference:\nCategory")+ #legend title
  theme_classic()+
  #coord_cartesian(ylim=c(0,6000))+
  theme(axis.text.x = element_text(angle = 45, hjust=1))+
  ggtitle("Relationship between RNA and Protein Difference\nupon chromosome arm gain: per Category")
#4x4
# plot.Protein.RNA.Bargraph.ThreeCat.Gain


####
###Bar graph of RNA and Protein quantiles Scaling/Buffered
## RNA & Protein Gain
ggplot(data= CN.Diff.xRNA.yProt.ThreeGroups, aes(x=Three.RNA.Gain, fill=Three.Protein.Gain)) + 
  geom_bar(position="fill")+
  scale_fill_manual(values=c("salmon2", "Gold2", "palegreen3"))+
  xlab("RNA Difference upon chrm arm gain: category")+
  ylab("Percent of genes")+
  labs(fill = "Protein Difference:\nCategory")+ #legend title
  theme_classic()+
  theme(axis.text.x = element_text(angle = 45, hjust=1))+
  ggtitle("Relationship between RNA and Protein Difference\nupon chromosome arm gain: per Category")
#4x4
# plot.Protein.RNA.Bargraph.ThreeCat.Gain

#RNA & Protein Loss
ggplot(data= CN.Diff.xRNA.yProt.ThreeGroups, aes(x=Three.RNA.Loss, fill=Three.Protein.Loss)) + 
  geom_bar(position="fill")+
  scale_fill_manual(values=c("salmon2", "Gold2", "palegreen3"))+
  xlab("RNA Difference upon chrm arm loss: category")+
  ylab("Percent of genes")+
  labs(fill = "Protein Difference:\nCategory")+ #legend title
  theme_classic()+
  theme(axis.text.x = element_text(angle = 45, hjust=1))+
  ggtitle("Relationship between RNA and Protein Difference\nupon chromosome arm loss: per Category")
#4x4
# plot.Protein.RNA.Bargraph.ThreeCat.Loss

## RNA Gain & Loss
ggplot(data= CN.Diff.xRNA.yProt.ThreeGroups, aes(x=Three.RNA.Gain, fill=Three.RNA.Loss)) + 
  geom_bar(position="fill")+
  scale_fill_manual(values=c("salmon2", "Gold2", "palegreen3"))+
  xlab("RNA Difference upon chrm arm gain: category")+
  ylab("Percent of genes")+
  labs(fill = "RNA Difference upon loss:\nCategory")+ #legend title
  theme_classic()+
  theme(axis.text.x = element_text(angle = 45, hjust=1))+
  ggtitle("Relationship between RNA genes \nupon chromosome arm gain or loss")
#4x4
# plot.Protein.RNA.Bargraph.ThreeCat.Gain

## Protein Gain & Loss
ggplot(data= CN.Diff.xRNA.yProt.ThreeGroups, aes(x=Three.Protein.Gain, fill=Three.Protein.Loss)) + 
  geom_bar(position="fill")+
  scale_fill_manual(values=c("salmon2", "Gold2", "palegreen3"))+
  xlab("Protein difference upon chrm arm gain: category")+
  ylab("Percent of genes")+
  labs(fill = "Protein Difference upon loss:\nCategory")+ #legend title
  theme_classic()+
  theme(axis.text.x = element_text(angle = 45, hjust=1))+
  ggtitle("Relationship between protein genes \nupon chromosome arm gain or loss")
#4x4
# plot.Protein.RNA.Bargraph.ThreeCat.Gain


###### Protein vs RNA gain/loss correlation graphs ######
## density plot of change in RNA & Protein upon chrm gain/Loss 
## using only data with 10+ datapoints per condition. 

CN.Diff.xRNA.yProt.ThreeGroups
#9413 genes analyzed
# CHROMOSOME GAIN: correlate RNA and Protein expression difference
# Mimimum of 10 datapoints per condition (gain, no aneu, loss)
ggplot(CN.Diff.xRNA.yProt.ThreeGroups, aes(x=CN.Diff.xRNA.yProt.ThreeGroups$RNA.Diff.Gain, y=CN.Diff.xRNA.yProt.ThreeGroups$Protein.Diff.Gain))+
  stat_density_2d(alpha=0.3, geom= "polygon", color="black",
                  aes() )+
  xlab("Difference in RNA expression")+
  ylab("Difference in protein expression")+
  theme_classic()+
  geom_hline(yintercept=0.0)+
  geom_vline(xintercept=0.0)+
  geom_smooth(method="lm", color="Red")+
  coord_cartesian(xlim=c(-0.7, 0.7), ylim=c(-0.7,0.7))+
  ggtitle("Chromosome gain")
#4x4
# plot.Protein.RNA.expression.density.Gain
Chrm_Gain_min10<-cor.test(CN.Diff.xRNA.yProt.ThreeGroups$Protein.Diff.Gain, CN.Diff.xRNA.yProt.ThreeGroups$RNA.Diff.Gain, method="pearson")
#cor=0.546
#P < 2E-16


# CHOMOSOME LOSS: correlate RNA and Protein expression difference
ggplot(CN.Diff.xRNA.yProt.ThreeGroups, aes(x=RNA.Diff.Loss, y=Protein.Diff.Loss))+
  stat_density_2d(alpha=0.3, geom= "polygon", color="black",
                  aes() )+
  xlab("Difference in RNA expression")+
  ylab("Difference in protein expression")+
  theme_classic()+
  geom_hline(yintercept=0.0)+
  geom_vline(xintercept=0.0)+
  geom_smooth(method="lm", color="Red")+
  coord_cartesian(xlim=c(-0.7, 0.7), ylim=c(-0.7,0.7))+
  ggtitle("Chromosome loss")
#4x4
#plot.Protein.RNA.expression.density.Loss
Chrm_Loss_min10<-cor.test(CN.Diff.xRNA.yProt.ThreeGroups$RNA.Diff.Loss, CN.Diff.xRNA.yProt.ThreeGroups$Protein.Diff.Loss, method="pearson")
Chrm_Loss_min10
#cor=0.554
#P<2E-16




###### Define function: plot Protein expression by chrm CN ####
### Step 2: make function to plot protein expression in mono, di and triploid cells 
## ProtExp.ChrmCN.filtered
## for Protein Data:
ProtExp.ChrmCN.filtered<- function(Protein){
  
  TestProt=Protein
  
  testProtChrm <- filter(Protein_Info4, Protein_Info4$"Approved Symbol"==TestProt) #get test protein data
  TestArm <- gsub('[0-9]+', '', testProtChrm$'Chromosome band') #Find test protein chromosome arm
  TestArm <- gsub('[.]', '', TestArm) #also remove period, if needed
  TestChrm<- as.numeric(testProtChrm$Chromosome) #Find test protein chromosome
  #filter cell data and get only aneuploidy scores for chrm arm of test protein location
  testchrm.percell <- filter(aneuploid, chrom==TestChrm) 
  testchrm.percell <- filter(testchrm.percell, arm==TestArm)
  
  
  TestProteinInfo <- filter(Protein_ProID, Gene_Symbol==TestProt)
  TestProteinID<- TestProteinInfo$Protein_Id
  
  Protein_Expression_TestProt <- Protein.Expression.filtered %>% select(Cell_line, TestProteinID)
  
  #Combine protein data with aneuploidy data for test protein location
  Protein.effect.Chrm<-merge(y= Protein_Expression_TestProt, x= testchrm.percell, 
                             by.y="Cell_line", by.x="DepMap_ID", 
                             sort = TRUE)# 368 cells, 12757 proteins 
  
  ## put data into categories based on chrm arm number
  Chrm.tri<- filter(Protein.effect.Chrm, arm_call==1) #all cells triploid for testProt chrm
  Chrm.di<- filter(Protein.effect.Chrm, arm_call==0) #all cells diploid for testProt chrm
  Chrm.mono<- filter(Protein.effect.Chrm, arm_call==-1) #all cells mono for testProt chrm
  
  ##Make data frame with t-test info about protein expression per arm number category.
  ## Return this data at the end of the function. 
  if (length(Chrm.tri[!is.na(Chrm.tri)])>=3 & 
      #Added if statement so only genes with 2+ values per condition are analyzed
      #if I don't do this, the t-test crashes and I get no values. 
      length(Chrm.di[!is.na(Chrm.di)])>=3 &
      length(Chrm.di[!is.na(Chrm.mono)])>=3) {
    # Trisomy vs disomy
    di.tri<-t.test(Chrm.tri[,8], # [,8] because that is column with Protein data
                   Chrm.di[,8], # [,8] because that is column with Protein data
                   mu = 0, 
                   alt = "two.sided",
                   conf.level = 0.99) #get p-value
    Diff.Tri.Di<-mean(Chrm.tri[,8], na.rm=TRUE) - mean(Chrm.di[,8], na.rm=TRUE) #get difference
    #Disomy vs monosome
    Mono.Di<-t.test(Chrm.mono[,8], # [,8] because that is column with Protein data
                    Chrm.di[,8], # [,8] because that is column with Protein data
                    mu = 0, 
                    alt = "two.sided",
                    conf.level = 0.99) #get p-value
    Diff.Mono.Di<-mean(Chrm.mono[,8], na.rm=TRUE) - mean(Chrm.di[,8], na.rm=TRUE)  #get difference
    # Tri vs monosomy
    Tri.Mono<-t.test(Chrm.mono[,8], # [,8] because that is column with Protein data
                     Chrm.tri[,8], # [,8] because that is column with Protein data
                     mu = 0, 
                     alt = "two.sided",
                     conf.level = 0.99) #get p-value
    Diff.Tri.Mono<-mean(Chrm.tri[,8], na.rm=TRUE) - mean(Chrm.mono[,8], na.rm=TRUE) #get difference
    
    t.test_prot_category<- data.frame(Protein_ID=colnames(Protein.effect.Chrm[8]), 
                                      Protein_Name=TestProt,
                                      Pvlaue.Tri.Di= di.tri$p.value,
                                      Diff.Gain= Diff.Tri.Di, 
                                      Pvlaue.Loss= Mono.Di$p.value,
                                      Diff.Loss= Diff.Mono.Di, 
                                      Pvlaue.Tri.Mono= Tri.Mono$p.value,
                                      Diff.GainvLoss= Diff.Tri.Mono)
  } else {
    t.test_prot_category<- c("Not enough cells to perform t-test on all conditions") 
  }
  
  
  ##Make Plot of protein expression per aneuploid category
  setwd(DataFileLocation) 
  pdf(paste("Plot.",TestProt,".Protein.Expression.per.Tri.Di.Mono.pdf", sep=''), width=4, height=4)## Save as PDF
  Plot.testProt<- ggplot(Protein.effect.Chrm,
                         aes(x = as_factor(arm_call), y=Protein.effect.Chrm[,8])) + 
    geom_boxplot(fill= c("-1"="dodgerblue3", "0"="grey80", "1"="red"), outlier.shape = NA) +
    geom_jitter() +
    xlab(paste("Chromosome arm call of ", TestProt," location: Chrm",TestChrm, TestArm))+
    ylab("Protein Expression per cell line") +
    theme_classic()+
    ggtitle(paste(TestProt," protein expression per cell\n per corresponding chromosome copy number"))
  print(Plot.testProt)
  dev.off()##stop saving PDF
  
  return(t.test_prot_category)
}

###### Define function: plot RNA expression by chrm CN #####
### Step 3: make function to plot RNA expression in mono, di and triploid cells 
## for RNA Data:
RNAExp.ChrmCN.filtered<- function(RNA){
  
  TestRNA=RNA
  
  testRNAChrm <- filter(RNA_Info3, RNA_Info3$"Approved Symbol"==TestRNA) #get test RNA data
  TestArm <- gsub('[0-9]+', '', testRNAChrm$'Chromosome band') #Find test RNA chromosome arm
  TestArm <- gsub('[.]', '', TestArm) #also remove period, if needed
  TestChrm<- as.numeric(testRNAChrm$Chromosome) #Find test RNA chromosome
  #filter cell data and get only aneuploidy scores for chrm arm of test RNA location
  testchrm.percell <- filter(aneuploid, chrom==TestChrm) 
  testchrm.percell <- filter(testchrm.percell, arm==TestArm)
  
  colnames(RNA.Expression.filtered)<-sub("[..].*", "", as.character(colnames(RNA.Expression.filtered)))#12755 genes
  
  RNA_Expression_TestRNA <- RNA.Expression.filtered %>% select(Cell_line, all_of(TestRNA))
  
  #Combine RNA data with aneuploidy data for test RNA location
  RNA.effect.Chrm<-merge(y= RNA_Expression_TestRNA, x= testchrm.percell, 
                         by.y="Cell_line", by.x="DepMap_ID", 
                         sort = TRUE)# 368 cells, 12757 RNAs 
  
  ## put data into categories based on chrm arm number
  Chrm.tri<- filter(RNA.effect.Chrm, arm_call==1) #all cells triploid for testRNA chrm
  Chrm.di<- filter(RNA.effect.Chrm, arm_call==0) #all cells diploid for testRNA chrm
  Chrm.mono<- filter(RNA.effect.Chrm, arm_call==-1) #all cells mono for testRNA chrm
  
  ##Make data frame with t-test info about RNA expression per arm number category.
  ## Return this data at the end of the function. 
  if (length(Chrm.tri[!is.na(Chrm.tri)])>=3 & 
      #Added if statement so only genes with 2+ values per condition are analyzed
      #if I don't do this, the t-test crashes and I get no values. 
      length(Chrm.di[!is.na(Chrm.di)])>=3 &
      length(Chrm.di[!is.na(Chrm.mono)])>=3) {
    # Trisomy vs disomy
    di.tri<-t.test(Chrm.tri[,8], # [,8] because that is column with RNA data
                   Chrm.di[,8], # [,8] because that is column with RNA data
                   mu = 0, 
                   alt = "two.sided",
                   conf.level = 0.99) #get p-value
    Diff.Tri.Di<-mean(Chrm.tri[,8], na.rm=TRUE) - mean(Chrm.di[,8], na.rm=TRUE) #get difference
    #Disomy vs monosome
    Mono.Di<-t.test(Chrm.mono[,8], # [,8] because that is column with RNA data
                    Chrm.di[,8], # [,8] because that is column with RNA data
                    mu = 0, 
                    alt = "two.sided",
                    conf.level = 0.99) #get p-value
    Diff.Mono.Di<- mean(Chrm.mono[,8], na.rm=TRUE) - mean(Chrm.di[,8], na.rm=TRUE)  #get difference
    # Tri vs monosomy
    Tri.Mono<-t.test(Chrm.mono[,8], # [,8] because that is column with RNA data
                     Chrm.tri[,8], # [,8] because that is column with RNA data
                     mu = 0, 
                     alt = "two.sided",
                     conf.level = 0.99) #get p-value
    Diff.Tri.Mono<-mean(Chrm.tri[,8], na.rm=TRUE) - mean(Chrm.mono[,8], na.rm=TRUE) #get difference
    
    t.test_RNA_category<- data.frame(RNA_ID=colnames(RNA.effect.Chrm[8]),
                                     RNA_Name=TestRNA,
                                     Pvlaue.Tri.Di= di.tri$p.value,
                                     Diff.Gain= Diff.Tri.Di, 
                                     Pvlaue.Loss= Mono.Di$p.value,
                                     Diff.Loss= Diff.Mono.Di, 
                                     Pvlaue.Tri.Mono= Tri.Mono$p.value,
                                     Diff.GainvLoss= Diff.Tri.Mono)
  } else {
    t.test_RNA_category<- c("Not enough cells to perform t-test on all conditions") 
  }
  
  
  ##Make Plot of RNA expression per aneuploid category
  setwd(DataFileLocation)
  pdf(paste("Plot.",TestRNA,".RNA.Expression.per.Tri.Di.Mono.pdf", sep=''), width=4, height=4)## Save as PDF
  Plot.testRNA<- ggplot(RNA.effect.Chrm,
                        aes(x = as_factor(arm_call), y=RNA.effect.Chrm[,8])) + 
    geom_boxplot(fill= c("-1"="dodgerblue3", "0"="grey80", "1"="red"), outlier.shape = NA) +
    geom_jitter() +
    xlab(paste("Chromosome arm call of ", TestRNA," location: Chrm",TestChrm, TestArm))+
    ylab("RNA Expression per cell line") +
    theme_classic()+
    ggtitle(paste(TestRNA," RNA expression per cell\n per corresponding chromosome copy number"))
  print(Plot.testRNA)
  dev.off()##stop saving PDF
  
  return(t.test_RNA_category)
}

###### Boxplots of specific gene expression per gain/neutral/loss ####

## plot specific protein expression per category
## plot and get p-values for Protein/RNA expression changes by chromosome CN
p.P53<-ProtExp.ChrmCN.filtered("TP53")
p.MDM2<-ProtExp.ChrmCN.filtered("MDM2")
p.MDM4<-ProtExp.ChrmCN.filtered("MDM4") #Significant Difference! 
p.CDK1<-ProtExp.ChrmCN.filtered("CDK1")
p.CDK13<-ProtExp.ChrmCN.filtered("CDK13") 
p.CDK12<-ProtExp.ChrmCN.filtered("CDK12")

## RNA expression
r.P53<-RNAExp.ChrmCN.filtered("TP53") #Significant Difference!
r.MDM2<-RNAExp.ChrmCN.filtered("MDM2") #Significant Difference!
r.MDM4<-RNAExp.ChrmCN.filtered("MDM4") #Significant Difference!
r.CDK1<-RNAExp.ChrmCN.filtered("CDK1") #Significant Difference!
r.CDK13<-RNAExp.ChrmCN.filtered("CDK13") #NS
r.CDK12<-RNAExp.ChrmCN.filtered("CDK12") #NS

## Found top 5 deregulated Proteins and RNA, by p-value, betwaan Mono and Tri. 
##plot Top 5 Protein
t.test_prot_category[1:5,2]
p.PTPN2<-ProtExp.ChrmCN.filtered("PTPN2")
r.PTPN2<-RNAExp.ChrmCN.filtered("PTPN2")
p.SMCHD1<-ProtExp.ChrmCN.filtered("SMCHD1")
r.SMCHD1<-RNAExp.ChrmCN.filtered("SMCHD1")
p.RNMT<-ProtExp.ChrmCN.filtered("RNMT")
r.RNMT<-RNAExp.ChrmCN.filtered("RNMT")
p.USP14<-ProtExp.ChrmCN.filtered("USP14")
r.USP14<-RNAExp.ChrmCN.filtered("USP14")
p.PPP4R1<-ProtExp.ChrmCN.filtered("PPP4R1")
r.PPP4R1<-RNAExp.ChrmCN.filtered("PPP4R1")

##plot Top 5 RNA
t.test_RNA_category[1:5,1]
p.NDUFAF5<-ProtExp.ChrmCN.filtered("NDUFAF5")
r.NDUFAF5<-RNAExp.ChrmCN.filtered("NDUFAF5")
p.NDUFV2<-ProtExp.ChrmCN.filtered("NDUFV2")
r.NDUFV2<-RNAExp.ChrmCN.filtered("NDUFV2")
p.MGME1<-ProtExp.ChrmCN.filtered("MGME1")
r.MGME1<-RNAExp.ChrmCN.filtered("MGME1")
p.ESF1<-ProtExp.ChrmCN.filtered("ESF1")
r.ESF1<-RNAExp.ChrmCN.filtered("ESF1")
p.MAVS<-ProtExp.ChrmCN.filtered("MAVS")
r.MAVS<-RNAExp.ChrmCN.filtered("MAVS")

## Plot TSG scaling w/gain loss at RNA
## Onco and TSG genes should be in Bailey lists

##Tumor supressor gene (TSG) and oncogene (OG) list
## Download Bailey et al. 2018 supplemental figure 8, upload table 1. example below: 
# Bailey et al. Comprehensive Characterization of Cancer Driver Genes and Mutations, cell, 2018

Oncogene.Bailey.list #get all oncogenes in Bailey et al 
TSG.Bailey.list # get all TSG in Bailey et al 

## Oncogenes
## Bailey
Oncog.Protein.ttest<-subset(t.test_prot_category, Protein_Name %in% Oncogene.Bailey.list)
Oncog.Protein.ttest<- Oncog.Protein.ttest[order(Oncog.Protein.ttest$Diff.Di_Mono),]
# most buffeted to most scaling (Tri and mono)--No sig, 
# Most scaling to most buffeted (Tri and mono)-- No Sig
# setwd(DataFileLocation)
# write_csv(Oncog.Protein.ttest, "Devoli.Onco.csv")
# Q1, Q4- No Sig
oncog.RNA.ttest<- subset(t.test_RNA_category, RNA_ID %in% Oncogene.Bailey.list)
oncog.RNA.ttest<- oncog.RNA.ttest[order(oncog.RNA.ttest$Diff.Gain),]
# 25 out of 39 (64%) genes not differently expressed between chrm gain or loss (protein level)
# 9 out of 37 (24%) genes not differenetly expressed between chrm gain or loss (RNA level)

p.BRAF<-ProtExp.ChrmCN.filtered("BRAF") #Sig (Most sig onco protein)
r.BRAF<-RNAExp.ChrmCN.filtered("BRAF") #Sig (tri)
p.MAPK1<-ProtExp.ChrmCN.filtered("MAPK1") #SIG
r.MAPK1<-RNAExp.ChrmCN.filtered("MAPK1") #SIG 
p.RRAS2<-ProtExp.ChrmCN.filtered("RRAS2") #SIG (Tri)
r.RRAS2<-RNAExp.ChrmCN.filtered("RRAS2") #SIG (Tri)
p.MTOR<-ProtExp.ChrmCN.filtered("MTOR") #NS
r.MTOR<-RNAExp.ChrmCN.filtered("MTOR") #SIG!
p.MYC<-ProtExp.ChrmCN.filtered("MYC") #NS
r.MYC<-RNAExp.ChrmCN.filtered("MYC") #NS

## TSG
# Bailey
#TSG.list2<-TSG.list[-c(34,35)] #Remove HLA-A and HLA-B

TSG.Protein.ttest<-subset(t.test_prot_category, Protein_Name %in% TSG.Devoli.list)
TSG.Protein.ttest<- TSG.Protein.ttest[order(TSG.Protein.ttest$Diff.Di_Mono),]
# most buffeted to most scaling (Tri and mono)-- NS, NS
# Most scaling to most buffeted (Tri and mono)-- NS, NS

TSG.RNA.ttest<- subset(t.test_RNA_category, RNA_ID %in% TSG.Bailey.list)
TSG.RNA.ttest<- TSG.RNA.ttest[order(TSG.RNA.ttest$Diff.Gain),]


# 35 of 66 genes (53%) are Not differently different expressed between gain or loss (protein level)
# 24 of 69 genes (35%) are Not differently different expressed between gain or loss (RNA level)
p.SMAD4<-ProtExp.ChrmCN.filtered("SMAD4") #Sig
r.SMAD4<-RNAExp.ChrmCN.filtered("SMAD4") #Sig
p.ATM<-ProtExp.ChrmCN.filtered("ATM") #Sig
r.ATM<-RNAExp.ChrmCN.filtered("ATM") #Sig
p.TP53<-ProtExp.ChrmCN.filtered("TP53") #NS
r.TP53<-RNAExp.ChrmCN.filtered("TP53") # Sig (mono)
p.MAP3K1<-ProtExp.ChrmCN.filtered("MAP3K1") #NS
r.MAP3K1<-RNAExp.ChrmCN.filtered("MAP3K1") #slight sig tri
p.CTNND1<-ProtExp.ChrmCN.filtered("CTNND1") #NS
r.CTNND1<-RNAExp.ChrmCN.filtered("CTNND1") #NS
p.KANSL1<-ProtExp.ChrmCN.filtered("KANSL1") #NS
r.KANSL1<-RNAExp.ChrmCN.filtered("KANSL1") #Sig
p.CDKN1A<-ProtExp.ChrmCN.filtered("CDKN1A") #NS
r.CDKN1A<-RNAExp.ChrmCN.filtered("CDKN1A") #NS
p.BRCA1<-ProtExp.ChrmCN.filtered("BRCA1") #NS
r.BRCA1<-RNAExp.ChrmCN.filtered("BRCA1") #Sig on mono
p.SMAD4<-ProtExp.ChrmCN.filtered("SMAD4") #NS
r.SMAD4<-RNAExp.ChrmCN.filtered("SMAD4") #Sig on mono

GainNS.LossNS<- subset(TSG.Protein.ttest, 
                       (Pvlaue.Tri.Di>0.05 &
                          Diff.Gain<0.2 &
                          Pvlaue.Di.Mono>0.05 &
                          Diff.Di_Mono> -0.2 ) )




### find gene increase upon gain, ns upon loss
GainSig.LossNS<- subset(CN.Diff.xRNA.yProt.ThreeGroups, 
                        (Protein.Pvalue.Gain<0.05 & #0r 0.05/9414
                           Protein.Diff.Gain>0.25 &
                           Protein.Pvalue.Loss>0.05 &
                           Protein.Diff.Loss>0 ) )
GainSig.LossNS<-GainSig.LossNS[order(GainSig.LossNS$Protein.Diff.Gain),]
p.KANSL1<-ProtExp.ChrmCN.filtered("TUBA3C") # gain sig, loss ns
r.KANSL1<-RNAExp.ChrmCN.filtered("TUBA3C") # gain sig, loss ns
p.KANSL1<-ProtExp.ChrmCN.filtered("GOLGA2") # gain sig, loss ns
r.KANSL1<-RNAExp.ChrmCN.filtered("GOLGA2") # gain sig, loss ns

FindPercent<- subset(CN.Diff.xRNA.yProt.ThreeGroups, 
                     (Protein.Pvalue.Gain<0.05 & #Sig 
                        Protein.Diff.Gain>0 & #increase
                        Protein.Pvalue.Loss<0.05 &
                        Protein.Diff.Loss<0 &
                        RNA.Pvalue.Gain<0.05 & 
                        RNA.Diff.Gain>0 &
                        RNA.Pvalue.Loss<0.05 &
                        RNA.Diff.Loss<0 ) )
(length(FindPercent$RNA_Name)*100)/length(CN.Diff.xRNA.yProt.ThreeGroups$RNA_Name)
# 19% no sig RNA nor Protein, gain and loss
# 9 % sig RNA and protein, gain and loss
# 13% sig RNA, not protein, gain and loss
# 6 % sig gain only, RNA and Protein
# 9 % sig loss only, RNA and Protein
# 7 % sig RNA, Protein gain only
# 11% sig RNA, Protein loss only
# 7 % sig RNA gain only
# 9 % sig RNA loss only
# 1.2 % sig Protein gain only
# 1.3 % sig Protein loss only
# 3.4 % Sig anti-scaling at RNA or protein, gain or loss

# 16% Protein gain sig, loss not (regardless of RNA)
# 22% Protein loss sig, gain not (regardless of RNA)


###### Make RNA/Protein difference & t-test list per Diploid, triploid, tetraploid ####

### subset 367 cells into ploidy
# change this for diploid, triploid or tetraploid subsets: 
FilteredCellPloidyData <- subset(aneuploid, DepMap_ID %in% Protein.Expression.filtered$Cell_line)
FilteredCellPloidyData <- FilteredCellPloidyData[,c(2,6)]
FilteredCellPloidyData <- unique(FilteredCellPloidyData)

#MonoploidCells <- subset(FilteredCellPloidyData, ploidy<1.5) # 7, 2%
DiploidCells <- subset(FilteredCellPloidyData, ploidy>=1.5 & ploidy<2.5) # 185, 50%
TriploidCells <- subset(FilteredCellPloidyData, ploidy>=2.5 & ploidy<3.5) # 137, 37%
TetraploidCells <- subset(FilteredCellPloidyData, ploidy>=3.5 & ploidy<4.5) # 24, 7%
#HighPloidyCells <- subset(FilteredCellPloidyData, ploidy>4.5) # 14, mean ploidy =5.7, 4%



### RNA difference & p-value in near-diploid, near-triploid and near-tetraploid aneuploidy
colnames(RNA.Expression.filtered)<-sub("[..].*", "", as.character(colnames(RNA.Expression.filtered)))#12755 genes
RNA.Expression.ByPloidy<-RNA.Expression.filtered[,-c(1)]

#set cell subgroup by ploidy: 
RNA.Expression.ByPloidy<- subset(RNA.Expression.ByPloidy, Cell_line %in% TetraploidCells$DepMap_ID)

t.test_RNA_category_byPloidy<- data.frame(RNA_ID=character(),
                                 RNA_Name=character(),
                                 Pvalue.Gain= numeric(),
                                 Diff.Gain= numeric(), 
                                 Pvalue.Loss= numeric(),
                                 Diff.Loss= numeric(),
                                 Pvalue.GainvLoss= numeric(),
                                 Diff.GainvLoss= numeric())

for (i in 2:length(RNA.Expression.ByPloidy)){
  
  TestRNA_name=sub("[..].*", "", as.character(colnames(RNA.Expression.ByPloidy[i])))
  TestRNA=colnames(RNA.Expression.ByPloidy[i])
  
  testRNAChrm <- filter(RNA_Info3, RNA_Info3$"Approved Symbol"==TestRNA_name) #get test RNA data
  if (length(testRNAChrm$Chromosome)!=0){ #Make sure RNA data is in RNA_info3, else it crashes
    TestArm <- testRNAChrm$arm #Find test RNA chromosome arm
    TestChrm<- as.numeric(testRNAChrm$Chromosome) #Find test RNA chromosome
    #filter cell data and get only aneuploidy scores for chrm arm of test RNA location
    testchrm.percell <- filter(aneuploid, chrom==TestChrm) 
    testchrm.percell <- filter(testchrm.percell, arm==TestArm)
    
    RNA_Expression_TestRNA <- RNA.Expression.ByPloidy %>% select(Cell_line, all_of(TestRNA))
    
    #Combine RNA data with aneuploidy data for test RNA location
    RNA.effect.Chrm<-merge(y= RNA_Expression_TestRNA, x= testchrm.percell, 
                           by.y="Cell_line", by.x="DepMap_ID", 
                           sort = TRUE)# 368 cells, 12757 RNAs 
    
    ## put data into categories based on chrm arm number
    Chrm.tri<- filter(RNA.effect.Chrm, arm_call==1) #all cells triploid for testRNA chrm
    Chrm.di<- filter(RNA.effect.Chrm, arm_call==0) #all cells diploid for testRNA chrm
    Chrm.mono<- filter(RNA.effect.Chrm, arm_call==-1) #all cells mono for testRNA chrm
    
    ##Make data frame with t-test info about RNA expression per arm number category.
    ## Return this data at the end of the function. 
    if (colSums(!is.na(Chrm.tri[8]))>=3 &  
        #Added if statement so only genes with 2+ values per condition are analyzed
        #if I don't do this, the t-test crashes and I get no values. 
        colSums(!is.na(Chrm.di[8]))>=3 & 
        colSums(!is.na(Chrm.mono[8]))>=3) { 
      # Trisomy vs disomy
      di.tri<-t.test(Chrm.tri[,8], # [,8] because that is column with RNA data
                     Chrm.di[,8], # [,8] because that is column with RNA data
                     mu = 0, 
                     alt = "two.sided",
                     conf.level = 0.99) #get p-value
      Diff.Tri.Di<-mean(Chrm.tri[,8], na.rm=TRUE) - mean(Chrm.di[,8], na.rm=TRUE) #get difference
      #Disomy vs monosome
      Mono.Di<-t.test(Chrm.mono[,8], # [,8] because that is column with RNA data
                      Chrm.di[,8], # [,8] because that is column with RNA data
                      mu = 0, 
                      alt = "two.sided",
                      conf.level = 0.99) #get p-value
      Diff.Mono.Di<- mean(Chrm.mono[,8], na.rm=TRUE) - mean(Chrm.di[,8], na.rm=TRUE) #get difference
      # Tri vs monosomy
      Tri.Mono<-t.test(Chrm.mono[,8], # [,8] because that is column with RNA data
                       Chrm.tri[,8], # [,8] because that is column with RNA data
                       mu = 0, 
                       alt = "two.sided",
                       conf.level = 0.99) #get p-value
      Diff.Tri.Mono<-mean(Chrm.tri[,8], na.rm=TRUE) - mean(Chrm.mono[,8], na.rm=TRUE) #get difference
      
      t.test_RNA_category_byPloidy<- rbind(t.test_RNA_category_byPloidy, data.frame(
        RNA_ID=colnames(RNA.effect.Chrm[8]),
        RNA_Name=TestRNA_name,
        Pvalue.Gain= di.tri$p.value,
        Diff.Gain= Diff.Tri.Di, 
        Pvalue.Loss= Mono.Di$p.value,
        Diff.Loss= Diff.Mono.Di, 
        Pvalue.GainvLoss= Tri.Mono$p.value,
        Diff.GainvLoss= Diff.Tri.Mono))
    } else {
      t.test_RNA_category_byPloidy<-t.test_RNA_category_byPloidy
    }
  }
}
t.test_RNA_category_byPloidy<-distinct(t.test_RNA_category_byPloidy) 
# diploid: 10841 RNA (17774 genes in all data)
# triploid: 16434 RNA
# tetraploid: 9579 RNA (min 3)
t.test_RNA_category_byPloidy<-t.test_RNA_category_byPloidy[order(t.test_RNA_category_byPloidy$Pvalue.GainvLoss),]


### Protein difference & p-value in near-diploid, near-triploid and near-tetraploid aneuploidy 
## make list of pvalue and difference between mono- di-tri cells per gene 
## for each gene in Protein data. using only cells that have both RNA and Protein data

Protein.Expression.byPloidy<- subset(Protein.Expression.filtered, 
                                     Cell_line %in% TetraploidCells$DepMap_ID)

Protein.ID<-data.frame(Protein_ID=
                         colnames(Protein.Expression.filtered[3:length(Protein.Expression.filtered)]))

##Now merge depmap column name info with Protein info. 
## First add Uniprot ID and gene symbol info to depmap data. 
Depmap.Protein.info2<-merge(x= Protein.ID, y= Protein_ProID, 
                            by.x="Protein_ID", by.y="Protein_Id", 
                            sort = TRUE)# 3 collumns, 12755. all Proteins given Uniprot IDs.
# Now combine all genes with same uniprot ID: 
Depmap.Protein.info3<-merge(x= Depmap.Protein.info2, y= Protein_Info4, 
                            by.x="Uniprot_Acc", by.y="UniProt ID(supplied by UniProt)", 
                            sort = TRUE)# 10594 Proteins
# Find those genes without uniprot id: 
No_UniprotID <- anti_join(Depmap.Protein.info2, Protein_Info4, #finding genes_Symbol without match
                          by = c("Uniprot_Acc" = "UniProt ID(supplied by UniProt)"))
#...find the genes with matching gene symbols
Depmap.Protein.info4<-merge(x= No_UniprotID, y= Protein_Info4, 
                            by.x="Gene_Symbol", by.y="Approved Symbol", 
                            sort = TRUE)# 10 collumns, 253 genes, only no gene-symbol genes with uniprot_ID

# Merge genes with gene symbol and genes with only uniprot ID. 
Depmap.Protein.info5<-merge(x= Depmap.Protein.info3, y= Depmap.Protein.info4, 
                            all=TRUE)# 11 collumns, 12100 genes


t.test_prot_category_byPloidy<- data.frame(Protein_ID=character(),
                                  Protein_Name=character(),
                                  Pvlaue.Tri.Di= numeric(),
                                  Diff.Gain= numeric(), 
                                  Pvlaue.Di.Mono= numeric(),
                                  Diff.Di_Mono= numeric(),
                                  Pvlaue.Tri.Mono= numeric(),
                                  Diff.GainvLoss= numeric())


for (i in 3:length(Protein.Expression.byPloidy)){
  TestProt = colnames(Protein.Expression.byPloidy[i])
  
  testProtChrm <- filter(Depmap.Protein.info5, Protein_ID==TestProt) #get test protein data
  
  
  #use If statement to check that Protein Info has Protein of interest. 
  if (length(testProtChrm[,1])!=0){
    TestArm <- testProtChrm$arm #Find test protein chromosome arm
    TestChrm<- as.numeric(testProtChrm$Chromosome) #Find test protein chromosome
    #filter cell data and get only aneuploidy scores for chrm arm of test protein location
    testchrm.percell <- filter(aneuploid, chrom==TestChrm) 
    testchrm.percell <- filter(testchrm.percell, arm==TestArm)
    
    Protein_Expression_TestProt <- Protein.Expression.byPloidy %>% select(Cell_line, all_of(TestProt))
    
    #Combine protein data with aneuploidy data for test protein location
    Protein.effect.Chrm<-merge(y= Protein_Expression_TestProt, x= testchrm.percell, 
                               by.y="Cell_line", by.x="DepMap_ID", 
                               sort = TRUE)# 368 cells, 12757 proteins 
    
    ## put data into categories based on chrm arm number
    Chrm.tri<- filter(Protein.effect.Chrm, arm_call==1) #all cells triploid for testProt chrm
    Chrm.di<- filter(Protein.effect.Chrm, arm_call==0) #all cells diploid for testProt chrm
    Chrm.mono<- filter(Protein.effect.Chrm, arm_call==-1) #all cells mono for testProt chrm
    
    ##Make data frame with t-test info about protein expression per arm number category.
    ## Return this data at the end of the function. 
    if (colSums(!is.na(Chrm.tri[8]))>=3 & #[8] is column with protein expression data. check it has values. 
        #Added if statement so only genes with 2+ values per condition are analyzed
        #if I don't do this, the t-test crashes and I get no values. 
        colSums(!is.na(Chrm.di[8]))>=3 & #usually 10, now 3**
        colSums(!is.na(Chrm.mono[8]))>=3) { #usually 10, now 3**
      # Trisomy vs disomy
      di.tri<-t.test(Chrm.tri[,8], # [,8] because that is column with Protein data
                     Chrm.di[,8], # [,8] because that is column with Protein data
                     mu = 0, 
                     alt = "two.sided",
                     conf.level = 0.99) #get p-value
      Diff.Tri.Di<-mean(Chrm.tri[,8], na.rm=TRUE) - mean(Chrm.di[,8], na.rm=TRUE) #get difference
      #Disomy vs monosome
      Di.Mono<-t.test(Chrm.mono[,8], # [,8] because that is column with Protein data
                      Chrm.di[,8], # [,8] because that is column with Protein data
                      mu = 0, 
                      alt = "two.sided",
                      conf.level = 0.99) #get p-value
      Diff.Di.Mono<- mean(Chrm.mono[,8], na.rm=TRUE) - mean(Chrm.di[,8], na.rm=TRUE)  #get difference
      # Tri vs monosomy
      Tri.Mono<-t.test(Chrm.mono[,8], # [,8] because that is column with Protein data
                       Chrm.tri[,8], # [,8] because that is column with Protein data
                       mu = 0, 
                       alt = "two.sided",
                       conf.level = 0.99) #get p-value
      Diff.Tri.Mono<-mean(Chrm.tri[,8], na.rm=TRUE) - mean(Chrm.mono[,8], na.rm=TRUE) #get difference
      
      t.test_prot_category_byPloidy<- rbind(t.test_prot_category_byPloidy, 
                                   data.frame(Protein_ID=TestProt,
                                              Protein_Name=testProtChrm$Gene_Symbol,
                                              Pvlaue.Tri.Di= di.tri$p.value,
                                              Diff.Gain= Diff.Tri.Di, 
                                              Pvlaue.Di.Mono= Di.Mono$p.value,
                                              Diff.Di_Mono= Diff.Di.Mono, 
                                              Pvlaue.Tri.Mono= Tri.Mono$p.value,
                                              Diff.GainvLoss= Diff.Tri.Mono))
    } else {
      t.test_prot_category_byPloidy<-t.test_prot_category_byPloidy
    }
  }
}

t.test_prot_category_byPloidy<-distinct(t.test_prot_category_byPloidy) 
# diploid= 4685 proteins
# triploid = 7738 proteins
# tetraploid = 4147 proteins
t.test_prot_category_byPloidy<-t.test_prot_category_byPloidy[order(t.test_prot_category_byPloidy$Pvlaue.Tri.Mono),]# *** genes


######          Combine RNA and Protein Data, save and calculate diff per ploidy group ####
## This is with filtered cell lines (only cells with RNA and Protein and aneuploid data:)
## This is with filtered genes: only genes with 10+ cells in all categories: RNA & protein gain, neutral and loss
setwd(DataFileLocation)

CN.Diff.RNA.Prot_byPloidy<-merge(x=t.test_RNA_category_byPloidy, 
                                 y=t.test_prot_category_byPloidy, 
                                 by.x="RNA_Name", by.y="Protein_Name") #4147

colnames(CN.Diff.RNA.Prot_byPloidy)<-c("RNA_Name", "RNA_ID", 
                                "RNA.Pvalue.Gain", "RNA.Diff.Gain", 
                                "RNA.Pvalue.Loss", "RNA.Diff.Loss", 
                                "RNA.Pvalue.GainvLoss", "RNA.Diff.GainvLoss", 
                                "Protein_ID", "Protein.Pvalue.Gain", 
                                "Protein.Diff.Gain", "Protein.Pvalue.Loss", 
                                "Protein.Diff.Loss", "Protein.Pvalue.GainvLoss", 
                                "Protein.Diff.GainvLoss")

write.csv(CN.Diff.RNA.Prot_byPloidy, #change name as needed ! for tetraploid triploid diploid
          file =paste("RNA.Protein_Loss.Neutral.Gain_Difference_Pvalue_Triploid_min3cells.csv", sep=','), 
          row.names = TRUE)

### Filter for same gene set 
x<- merge(x=CN.Diff.RNA.Prot_Triploid, y=CN.Diff.RNA.Prot_Diploid, 
          by.x="RNA_Name", by.y="RNA_Name")
Min10GeneSet_dividedbyPloidy<-x$RNA_Name


### Diploid 
#CN.Diff.RNA.Prot_Diploid<-CN.Diff.RNA.Prot_byPloidy
CN.Diff.RNA.Prot_Diploid2<-subset(CN.Diff.RNA.Prot_Diploid, RNA_Name %in% Min10GeneSet_dividedbyPloidy)
# cells= 185
# genes (>10 points/category, in RNA and Protein)= 4489
# Mean Prot gain: 0.188957, 14%
# Mean RNA  gain: 0.328812, 26%
# Mean Prot loss: -0.1717, -11%
# Mean RNA  loss: -0.3156, -19%

length(CN.Diff.RNA.Prot_Diploid$Protein.Diff.Gain)# 4489
sum(CN.Diff.RNA.Prot_Diploid$Protein.Diff.Gain >= log2(3/2)) #318, 7.1% >= DNA CN change upon gain
sum(CN.Diff.RNA.Prot_Diploid$Protein.Diff.Gain <= log2(1/2)) #5, 0.1% <= DNA CN change upon loss



### Triploid 
#CN.Diff.RNA.Prot_Triploid<-CN.Diff.RNA.Prot_byPloidy
#CN.Diff.RNA.Prot_Triploid2<-subset(CN.Diff.RNA.Prot_Triploid, RNA_Name %in% Min10GeneSet_dividedbyPloidy)

# cells= 137
# genes (>10 points/category)= 7411
# Mean Prot gain:  0.134, 9.7%
# Mean RNA  gain:  0.245, 19%
# Mean Prot loss: -0.151, -10%
# Mean RNA  loss: -0.288, -18%

# number of genes >= log2(4/3) expression (>= DNA CN change) **
length(CN.Diff.RNA.Prot_Triploid$Protein.Diff.Gain) #7411 proteins
sum(CN.Diff.RNA.Prot_Triploid$Protein.Diff.Gain >= log2(4/3)) #801, 11% >= DNA CN change upon gain
sum(CN.Diff.RNA.Prot_Triploid$Protein.Diff.Gain <= log2(2/3)) #82, 1.1% <= DNA CN change upon loss


### Tetraploid 
#CN.Diff.RNA.Prot_Tetraploid<-CN.Diff.RNA.Prot_byPloidy
CN.Diff.RNA.Prot_Tetraploid2<-subset(CN.Diff.RNA.Prot_Tetraploid, RNA_Name %in% Min10GeneSet_dividedbyPloidy)

# cells= 24
# genes (>3 points/category)= 3977
# Mean Prot gain:  10%
# Mean RNA  gain:  16%
# Mean Prot loss: -9%
# Mean RNA  loss: -18%

# number of genes >= log2(4/3) expression (>= DNA CN change) **
length(CN.Diff.RNA.Prot_Tetraploid2$Protein.Diff.Gain)# 2200
sum(CN.Diff.RNA.Prot_Tetraploid2$Protein.Diff.Gain >= log2(5/4)) #659, 30% >= DNA CN change upon gain
sum(CN.Diff.RNA.Prot_Tetraploid2$Protein.Diff.Gain <= log2(3/4)) #215, 9.7% <= DNA CN change upon loss


### All 
# CN.Diff.xRNA.yProt.ThreeGroups
CN.Diff.xRNA.yProt.ThreeGroups2<-subset(CN.Diff.xRNA.yProt.ThreeGroups, RNA_Name %in% Min10GeneSet_dividedbyPloidy)

# cells= 367
# genes (>10 points/category)= 9414
# Mean Prot gain: 12%
# Mean RNA  gain: 22%
# Mean Prot loss: -8.4%
# Mean RNA  loss: -15%

####Heatmap for mean expression upon gain and loss, protein/RNA
#calculate mean gain loss in DNA copy number upon chrm gain and loss, account for ploidy
log2(3/2)*(185/367) + log2(4/3)*(137/367) + log2(5/4)*(24/367) + log2(2/1)*(7/367) + log2(6.7/5.7)*(14/367)
log2(1/2)*(185/367) + log2(2/3)*(137/367) + log2(3/4)*(24/367) + log2(4.7/5.7)*(14/367)


#make dataframe for heatmap
Di.Tri.all.meanDiff<-data.frame(
  Cells=as.factor(c("Diploid", "Diploid", "Diploid", "Diploid", "Diploid", "Diploid", 
                       "Triploid", "Triploid", "Triploid", "Triploid", "Triploid", "Triploid", 
                       #"Tetraploid", "Tetraploid", "Tetraploid", "Tetraploid", "Tetraploid", "Tetraploid", 
                       "All", "All", "All", "All", "All", "All")),
  Condition=as.factor(c("Gain.DNA","Gain.RNA", "Gain.Protein", "Loss.DNA", "Loss.RNA", "Loss.Protein", 
                           "Gain.DNA","Gain.RNA", "Gain.Protein", "Loss.DNA", "Loss.RNA", "Loss.Protein", 
                           #"Gain.DNA","Gain.RNA", "Gain.Protein", "Loss.DNA", "Loss.RNA", "Loss.Protein", 
                           "Gain.DNA","Gain.RNA", "Gain.Protein", "Loss.DNA", "Loss.RNA", "Loss.Protein")), 
  values=as.numeric(c(log2(3/2), mean(CN.Diff.RNA.Prot_Diploid$RNA.Diff.Gain), mean(CN.Diff.RNA.Prot_Diploid$Protein.Diff.Gain), 
                      log2(1/2), mean(CN.Diff.RNA.Prot_Diploid$RNA.Diff.Loss), mean(CN.Diff.RNA.Prot_Diploid$Protein.Diff.Loss), 
                      log2(4/3), mean(CN.Diff.RNA.Prot_Triploid$RNA.Diff.Gain), mean(CN.Diff.RNA.Prot_Triploid$Protein.Diff.Gain), 
                      log2(2/3), mean(CN.Diff.RNA.Prot_Triploid$RNA.Diff.Loss), mean(CN.Diff.RNA.Prot_Triploid$Protein.Diff.Loss), 
                      #log2(5/4), mean(CN.Diff.RNA.Prot_Tetraploid$RNA.Diff.Gain), mean(CN.Diff.RNA.Prot_Tetraploid$Protein.Diff.Gain), 
                      #log2(3/4), mean(CN.Diff.RNA.Prot_Tetraploid$RNA.Diff.Loss), mean(CN.Diff.RNA.Prot_Tetraploid$Protein.Diff.Loss), 
                      0.4988, mean(CN.Diff.xRNA.yProt.ThreeGroups$RNA.Diff.Gain), mean(CN.Diff.xRNA.yProt.ThreeGroups$Protein.Diff.Gain), 
                      -0.76, mean(CN.Diff.xRNA.yProt.ThreeGroups$RNA.Diff.Loss), mean(CN.Diff.xRNA.yProt.ThreeGroups$Protein.Diff.Loss) 
  )
  )
)

Di.Tri.all.meanDiff$Condition<- factor(Di.Tri.all.meanDiff$Condition, levels=c("Loss.Protein", "Loss.RNA", "Loss.DNA", "Gain.Protein", "Gain.RNA", "Gain.DNA"))
Di.Tri.all.meanDiff$Cells<- factor(Di.Tri.all.meanDiff$Cells, levels=c("Diploid", "Triploid", "Tetraploid", "All"))

#plot heatmap of mean difference by ploidy: 
ggplot(Di.Tri.all.meanDiff, aes(x=Cells, y=Condition))+
  geom_raster(aes(fill = values), hjust=0.5, vjust=0.5, interpolate=FALSE)+
  xlab("Cancer cell line ploidy")+
  ylab("Gene expression difference")+
  theme_classic()+
  theme(axis.text.x = element_text(angle = 45, hjust=1, color="black"), 
        axis.text.y = element_text(color="black"))+
  scale_fill_gradient2(low="dodgerblue3", mid="white", high="Red", midpoint=0, 
                       name="Difference", limits=c(-1, 1)) +
  ggtitle("Mean difference per cell ploidy")
# 5x4
# plot.heatmap.meanDiff.Ploidy
# plot.heatmap.meanDiff.Ploidy_withTetra


###### Low/high Aneuploid cells only, still see buffering of ribosomes? ####

### Protein difference & p-value in near-diploid, near-triploid and near-tetraploid aneuploidy 
## make list of pvalue and difference between mono- di-tri cells per gene 
## for each gene in Protein data. using only cells that have both RNA and Protein data

### sub-step 1: Get protein and RNA expression data set for only cell lines with low aneuploidy: 

#List of low aneuploid cells from protein_expression_data_GeneScore.R
#quantile(Protein_Scores$Gene_ploidy_Score) #0, 1773, 2668, 3309, 6985
#low.aneuploid<- subset(Protein_Scores, Gene_ploidy_Score<= 1773) #only 94 cells
#low.aneuploid.cells<-low.aneuploid$Broad_ID
#high.aneuploid<- subset(Protein_Scores, Gene_ploidy_Score> 3309) #only 94 cells
#high.aneuploid.cells<-high.aneuploid$Broad_ID


Protein.Expression_lowaneuploidCells<- subset(Protein.Expression.filtered_min10Cells, Cell_line %in% low.aneuploid.cells)
#Protein.Expression_lowaneuploidCells<- subset(Protein.Expression.filtered_min10Cells, Cell_line %in% high.aneuploid.cells)

RNA.Expression_lowaneuploidCells<- subset(Protein.Expression.filtered_min10Cells, Cell_line %in% low.aneuploid.cells)
#RNA.Expression_lowaneuploidCells<- subset(Protein.Expression.filtered_min10Cells, Cell_line %in% high.aneuploid.cells)



### sub-step 2: get protein difference in cells with low aneuploidy
Protein.ID<-data.frame(Protein_ID=
                         colnames(Protein.Expression_lowaneuploidCells[3:length(Protein.Expression_lowaneuploidCells)]))

# Now merge depmap column name info with Protein info. 
# First add Uniprot ID and gene symbol info to depmap data. 
Depmap.Protein.info2<-merge(x= Protein.ID, y= Protein_ProID, 
                            by.x="Protein_ID", by.y="Protein_Id", 
                            sort = TRUE)# 3 collumns, 12755. all Proteins given Uniprot IDs.
# Now combine all genes with same uniprot ID: 
Depmap.Protein.info3<-merge(x= Depmap.Protein.info2, y= Protein_Info4, 
                            by.x="Uniprot_Acc", by.y="UniProt ID(supplied by UniProt)", 
                            sort = TRUE)# 10594 Proteins
# Find those genes without uniprot id: 
No_UniprotID <- anti_join(Depmap.Protein.info2, Protein_Info4, #finding genes_Symbol without match
                          by = c("Uniprot_Acc" = "UniProt ID(supplied by UniProt)"))
#...find the genes with matching gene symbols
Depmap.Protein.info4<-merge(x= No_UniprotID, y= Protein_Info4, 
                            by.x="Gene_Symbol", by.y="Approved Symbol", 
                            sort = TRUE)# 10 collumns, 253 genes, only no gene-symbol genes with uniprot_ID

# Merge genes with gene symbol and genes with only uniprot ID. 
Depmap.Protein.info5<-merge(x= Depmap.Protein.info3, y= Depmap.Protein.info4, 
                            all=TRUE)# 11 collumns, 12100 genes


t.test_prot_category_LowPloidy<- data.frame(Protein_ID=character(),
                                           Protein_Name=character(),
                                           Pvlaue.Tri.Di= numeric(),
                                           Diff.Gain= numeric(), 
                                           Pvlaue.Di.Mono= numeric(),
                                           Diff.Di_Mono= numeric(),
                                           Pvlaue.Tri.Mono= numeric(),
                                           Diff.GainvLoss= numeric())


for (i in 3:length(Protein.Expression_lowaneuploidCells)){
  TestProt = colnames(Protein.Expression_lowaneuploidCells[i])
  
  testProtChrm <- filter(Depmap.Protein.info5, Protein_ID==TestProt) #get test protein data
  
  
  #use If statement to check that Protein Info has Protein of interest. 
  if (length(testProtChrm[,1])!=0){
    TestArm <- testProtChrm$arm #Find test protein chromosome arm
    TestChrm<- as.numeric(testProtChrm$Chromosome) #Find test protein chromosome
    #filter cell data and get only aneuploidy scores for chrm arm of test protein location
    testchrm.percell <- filter(aneuploid, chrom==TestChrm) 
    testchrm.percell <- filter(testchrm.percell, arm==TestArm)
    
    Protein_Expression_TestProt <- Protein.Expression_lowaneuploidCells %>% select(Cell_line, all_of(TestProt))
    
    #Combine protein data with aneuploidy data for test protein location
    Protein.effect.Chrm<-merge(y= Protein_Expression_TestProt, x= testchrm.percell, 
                               by.y="Cell_line", by.x="DepMap_ID", 
                               sort = TRUE)# 368 cells, 12757 proteins 
    
    ## put data into categories based on chrm arm number
    Chrm.tri<- filter(Protein.effect.Chrm, arm_call==1) #all cells triploid for testProt chrm
    Chrm.di<- filter(Protein.effect.Chrm, arm_call==0) #all cells diploid for testProt chrm
    Chrm.mono<- filter(Protein.effect.Chrm, arm_call==-1) #all cells mono for testProt chrm
    
    ##Make data frame with t-test info about protein expression per arm number category.
    ## Return this data at the end of the function. 
    if (colSums(!is.na(Chrm.tri[8]))>=3 & #[8] is column with protein expression data. check it has values. 
        #Added if statement so only genes with 2+ values per condition are analyzed
        #if I don't do this, the t-test crashes and I get no values. 
        colSums(!is.na(Chrm.di[8]))>=3 & #usually 10, now 3**
        colSums(!is.na(Chrm.mono[8]))>=3) { #usually 10, now 3**
      # Trisomy vs disomy
      di.tri<-t.test(Chrm.tri[,8], # [,8] because that is column with Protein data
                     Chrm.di[,8], # [,8] because that is column with Protein data
                     mu = 0, 
                     alt = "two.sided",
                     conf.level = 0.99) #get p-value
      Diff.Tri.Di<-mean(Chrm.tri[,8], na.rm=TRUE) - mean(Chrm.di[,8], na.rm=TRUE) #get difference
      #Disomy vs monosome
      Di.Mono<-t.test(Chrm.mono[,8], # [,8] because that is column with Protein data
                      Chrm.di[,8], # [,8] because that is column with Protein data
                      mu = 0, 
                      alt = "two.sided",
                      conf.level = 0.99) #get p-value
      Diff.Di.Mono<- mean(Chrm.mono[,8], na.rm=TRUE) - mean(Chrm.di[,8], na.rm=TRUE)  #get difference
      # Tri vs monosomy
      Tri.Mono<-t.test(Chrm.mono[,8], # [,8] because that is column with Protein data
                       Chrm.tri[,8], # [,8] because that is column with Protein data
                       mu = 0, 
                       alt = "two.sided",
                       conf.level = 0.99) #get p-value
      Diff.Tri.Mono<-mean(Chrm.tri[,8], na.rm=TRUE) - mean(Chrm.mono[,8], na.rm=TRUE) #get difference
      
      t.test_prot_category_LowPloidy<- rbind(t.test_prot_category_LowPloidy, 
                                            data.frame(Protein_ID=TestProt,
                                                       Protein_Name=testProtChrm$Gene_Symbol,
                                                       Pvlaue.Tri.Di= di.tri$p.value,
                                                       Diff.Gain= Diff.Tri.Di, 
                                                       Pvlaue.Di.Mono= Di.Mono$p.value,
                                                       Diff.Di_Mono= Diff.Di.Mono, 
                                                       Pvlaue.Tri.Mono= Tri.Mono$p.value,
                                                       Diff.GainvLoss= Diff.Tri.Mono))
    } else {
      t.test_prot_category_LowPloidy<-t.test_prot_category_LowPloidy
    }
  }
}

t.test_prot_category_LowPloidy<-distinct(t.test_prot_category_LowPloidy) 

t.test_prot_category_LowPloidy<-t.test_prot_category_LowPloidy[order(t.test_prot_category_LowPloidy$Pvlaue.Tri.Mono),]# *** genes

# show ribosomes still buffered
# compare low aneuploid and high aneuploid cell lines difference expressions are correlated. 



### sub-step 3: get RNA difference in cells with low aneuploidy
# Low ploidy analysis: RNA expression difference: 
# prep data frame: 
colnames(RNA.Expression.filtered)<-sub("[..].*", "", as.character(colnames(RNA.Expression.filtered)))#12755 genes
RNA.Expression.LowPloidy<-RNA.Expression.filtered[,-c(1)]

RNA.Expression.LowPloidy<- subset(RNA.Expression.LowPloidy, Cell_line %in% high.aneuploid.cells)

# Run loop to get difference in experssion per gene in only low aneuploid cells
t.test_RNA_category_LowPloidy<- data.frame(RNA_ID=character(),
                                          RNA_Name=character(),
                                          Pvalue.Gain= numeric(),
                                          Diff.Gain= numeric(), 
                                          Pvalue.Loss= numeric(),
                                          Diff.Loss= numeric(),
                                          Pvalue.GainvLoss= numeric(),
                                          Diff.GainvLoss= numeric())

for (i in 2:length(RNA.Expression.LowPloidy)){
  
  TestRNA_name=sub("[..].*", "", as.character(colnames(RNA.Expression.LowPloidy[i])))
  TestRNA=colnames(RNA.Expression.LowPloidy[i])
  
  testRNAChrm <- filter(RNA_Info3, RNA_Info3$"Approved Symbol"==TestRNA_name) #get test RNA data
  if (length(testRNAChrm$Chromosome)!=0){ #Make sure RNA data is in RNA_info3, else it crashes
    TestArm <- testRNAChrm$arm #Find test RNA chromosome arm
    TestChrm<- as.numeric(testRNAChrm$Chromosome) #Find test RNA chromosome
    #filter cell data and get only aneuploidy scores for chrm arm of test RNA location
    testchrm.percell <- filter(aneuploid, chrom==TestChrm) 
    testchrm.percell <- filter(testchrm.percell, arm==TestArm)
    
    RNA_Expression_TestRNA <- RNA.Expression.LowPloidy %>% select(Cell_line, all_of(TestRNA))
    
    #Combine RNA data with aneuploidy data for test RNA location
    RNA.effect.Chrm<-merge(y= RNA_Expression_TestRNA, x= testchrm.percell, 
                           by.y="Cell_line", by.x="DepMap_ID", 
                           sort = TRUE)# 368 cells, 12757 RNAs 
    
    ## put data into categories based on chrm arm number
    Chrm.tri<- filter(RNA.effect.Chrm, arm_call==1) #all cells triploid for testRNA chrm
    Chrm.di<- filter(RNA.effect.Chrm, arm_call==0) #all cells diploid for testRNA chrm
    Chrm.mono<- filter(RNA.effect.Chrm, arm_call==-1) #all cells mono for testRNA chrm
    
    ##Make data frame with t-test info about RNA expression per arm number category.
    ## Return this data at the end of the function. 
    if (colSums(!is.na(Chrm.tri[8]))>=3 &  
        #Added if statement so only genes with 2+ values per condition are analyzed
        #if I don't do this, the t-test crashes and I get no values. 
        colSums(!is.na(Chrm.di[8]))>=3 & 
        colSums(!is.na(Chrm.mono[8]))>=3) { 
      # Trisomy vs disomy
      di.tri<-t.test(Chrm.tri[,8], # [,8] because that is column with RNA data
                     Chrm.di[,8], # [,8] because that is column with RNA data
                     mu = 0, 
                     alt = "two.sided",
                     conf.level = 0.99) #get p-value
      Diff.Tri.Di<-mean(Chrm.tri[,8], na.rm=TRUE) - mean(Chrm.di[,8], na.rm=TRUE) #get difference
      #Disomy vs monosome
      Mono.Di<-t.test(Chrm.mono[,8], # [,8] because that is column with RNA data
                      Chrm.di[,8], # [,8] because that is column with RNA data
                      mu = 0, 
                      alt = "two.sided",
                      conf.level = 0.99) #get p-value
      Diff.Mono.Di<- mean(Chrm.mono[,8], na.rm=TRUE) - mean(Chrm.di[,8], na.rm=TRUE) #get difference
      # Tri vs monosomy
      Tri.Mono<-t.test(Chrm.mono[,8], # [,8] because that is column with RNA data
                       Chrm.tri[,8], # [,8] because that is column with RNA data
                       mu = 0, 
                       alt = "two.sided",
                       conf.level = 0.99) #get p-value
      Diff.Tri.Mono<-mean(Chrm.tri[,8], na.rm=TRUE) - mean(Chrm.mono[,8], na.rm=TRUE) #get difference
      
      t.test_RNA_category_LowPloidy<- rbind(t.test_RNA_category_LowPloidy, data.frame(
        RNA_ID=colnames(RNA.effect.Chrm[8]),
        RNA_Name=TestRNA_name,
        Pvalue.Gain= di.tri$p.value,
        Diff.Gain= Diff.Tri.Di, 
        Pvalue.Loss= Mono.Di$p.value,
        Diff.Loss= Diff.Mono.Di, 
        Pvalue.GainvLoss= Tri.Mono$p.value,
        Diff.GainvLoss= Diff.Tri.Mono))
    } else {
      t.test_RNA_category_LowPloidy<-t.test_RNA_category_LowPloidy
    }
  }
}

t.test_RNA_category_LowPloidy<-distinct(t.test_RNA_category_LowPloidy) 
t.test_RNA_category_LowPloidy<-t.test_RNA_category_LowPloidy[order(t.test_RNA_category_LowPloidy$Pvalue.GainvLoss),]


### Sub step 4: combine low ploidy RNA and protein data
CN.Diff.RNA.Prot_HighPloidy<-merge(x=t.test_RNA_category_LowPloidy, 
                                 y=t.test_prot_category_LowPloidy, 
                                 by.x="RNA_Name", by.y="Protein_Name") #4147

colnames(CN.Diff.RNA.Prot_HighPloidy)<-c("RNA_Name", "RNA_ID", 
                                       "RNA.Pvalue.Gain", "RNA.Diff.Gain", 
                                       "RNA.Pvalue.Loss", "RNA.Diff.Loss", 
                                       "RNA.Pvalue.GainvLoss", "RNA.Diff.GainvLoss", 
                                       "Protein_ID", "Protein.Pvalue.Gain", 
                                       "Protein.Diff.Gain", "Protein.Pvalue.Loss", 
                                       "Protein.Diff.Loss", "Protein.Pvalue.GainvLoss", 
                                       "Protein.Diff.GainvLoss")

## add categories based on cutoffs.  (-Inf,-0.1,0.25,Inf)
CN.Diff.RNA.Prot_HighPloidy$Three.RNA.Gain<- cut(CN.Diff.RNA.Prot_HighPloidy$RNA.Diff.Gain,
                                                    breaks=c(-Inf,-0.1,0.25,Inf),
                                                    include.lowest=TRUE,
                                                    labels=c("Anti-Scaling","Buffering","Scaling"))
CN.Diff.RNA.Prot_HighPloidy$Three.RNA.Loss<- cut(CN.Diff.RNA.Prot_HighPloidy$RNA.Diff.Loss,
                                                    breaks=c(-Inf,-0.25,0.1,Inf),
                                                    include.lowest=TRUE,
                                                    labels=c("Scaling", "Buffering","Anti-Scaling"))
CN.Diff.RNA.Prot_HighPloidy$Three.Protein.Gain<- cut(CN.Diff.RNA.Prot_HighPloidy$Protein.Diff.Gain,
                                                        breaks=c(-Inf,-0.1,0.25,Inf),
                                                        include.lowest=TRUE,
                                                        labels=c("Anti-Scaling","Buffering","Scaling"))
CN.Diff.RNA.Prot_HighPloidy$Three.Protein.Loss<- cut(CN.Diff.RNA.Prot_HighPloidy$Protein.Diff.Loss,
                                                        breaks=c(-Inf,-0.25,0.1,Inf),
                                                        include.lowest=TRUE,
                                                        labels=c("Scaling", "Buffering","Anti-Scaling"))
CN.Diff.RNA.Prot_HighPloidy$Three.Protein.Loss<-factor(CN.Diff.RNA.Prot_HighPloidy$Three.Protein.Loss, 
                                                          levels=c("Anti-Scaling","Buffering","Scaling"))
CN.Diff.RNA.Prot_HighPloidy$Three.RNA.Loss<-factor(CN.Diff.RNA.Prot_HighPloidy$Three.RNA.Loss, 
                                                      levels=c("Anti-Scaling","Buffering","Scaling"))


#setwd(DataFileLocation)
write.csv(CN.Diff.RNA.Prot_LowPloidy, #change name as needed
          file =paste("RNA.Protein_Loss.Neutral.Gain_Difference_Pvalue_LowPloidy.min3cells.csv", sep=','), 
          row.names = TRUE)

#CN.Diff.RNA.Prot_LowPloidy<- read.csv("RNA.Protein_Loss.Neutral.Gain_Difference_Pvalue_LowPloidy.min10cells.csv")
#CN.Diff.RNA.Prot_LowPloidy<-CN.Diff.RNA.Prot_LowPloidy[,-1]

### Repeat above code with "High ploidy cells" to get difference upon gain in the top 1/4th high aneuploid cells 
write.csv(CN.Diff.RNA.Prot_HighPloidy, #change name as needed
          file =paste("RNA.Protein_Loss.Neutral.Gain_Difference_Pvalue_HighPloidy.min3cells.csv", sep=','), 
          row.names = TRUE)

#write.csv(Protein.effect.Chrm$DepMap_ID, #save a list of the cell lines used in analysis. 
#          file =paste("CellLines_ProdeinDosageCompensationManuscript.csv", sep=','), 
#          row.names = TRUE)
