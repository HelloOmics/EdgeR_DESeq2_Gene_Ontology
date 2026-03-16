
# Day 2 -------------------------------------------------------------------

# used via import as shown by Paulo
#library(readr)
#Matrix1 <- read_delim("Matrix1.txt", delim = "\t", 
#                      escape_double = FALSE, trim_ws = TRUE)

View(Matrix1)

Table1=read.table("C:/Users/roxan/OneDrive/Desktop/Masters_Biology_FU/Bioinf/Bioinf/Matrix2.txt", header=TRUE, sep="\t",
                  quote="", stringsAsFactor=FALSE)
Table1

#plotting:
# first define x and y achsis
x=Table1[,3]
y=Table1[,4]
plot(x,y)

plot(x,y, xlab="Expression1", ylab="Expression2")

plot(x,y, xlab="Expression1", ylab="Expression2", type="p",
     pch=20, col="red")

?plot

boxplot(y)
boxplot(y, xlab="Gene 2", ylab="Expression level") 


# Exercise VI -------------------------------------------------------------

#1. Some people prefer to have slides with a black background. Then it would 
#look nicer, if also your plot had a black background and if your axes and axis 
#labels would be in white. Hint: First, to change the color of your graphics 
#background use par(bg="black"). Then checkthe help pages for par() to see how 
#you can change the color for your axes and axis labels.2. Optional: Have a 
#look at the help pages for plot() and par() and play with your graphics to
#modify your plots a bit more, e.g. by changing colors, fonts etc. You can also 
#have a look atggplot to get an idea of what else is possible with R:
#  http://docs.ggplot2.org/current/

par(bg="black") # get black background!
par(mfrow = c(2,2), pty = "s") #get 4 panels!! 
?plot
plot(x,y, xlab="Expression1", ylab="Expression2", type="p",
     pch=20, col="red")
boxplot(y, xlab="Gene 2", ylab="Expression level")
plot(sin, from = -2*pi, to = 2*pi)
plot(1:10, 1:10)

par(mfrow = c(1,1), pty = "s")  # go back to one panel
plot(x,y, xlab="Expression1", ylab="Expression2", type="p",
     pch=20, col="red")

par(bg="white") # background back to white!
plot(x,y, xlab="Expression1", ylab="Expression2", type="p",
     pch=20, col="red")



# Introduction to Gene expression and network analysis --------------------

# dataset = social status alters immune reulation and response to infection in 
# macaques
# the lower your social rank the more stress you have 
# => the worse your immune system?
# took isolated blood cells from different individuals
# treated these blood cells with LPS to mimic an infection 
# the used RNE-seq to check gene expression

# quantification of gene expr4ssion
# amplification of RNA, the fragmentation of all RNA
# Cobversion of RNA to cDNA => get sequenced (short reads of about 150 nt long)
# => mapping to reference genome => exons/introns
# counts of reads => quantification of gene expression level
# can count the exons, transcipts or genes
# amplification process is the same so gene expressions will eb proportional
# if we just sum up the reads we would find that short genes have less reads 
# and addume that they are expressed less, but that not true = Normalisation!!
# RPKM => reads per kb of exon model per million mapped reads
# we no lobnger use RPKM in within species or indivduals
# however we do normalise for gene length when comparing between species i.e 
# humans and mice
# library sizes will also be slightly different because RNA pipetted more in 
# sample a than sample b
# how many RNAs and diversity of RNAs? also critical to potentially normalise
# in our case no problem comparing blood samples, but trasncriptome size a lot
# larger in brain vs liver
# 7 commonly used normalization methods for RNA-seq
# all methods come up with different results!!!
# compared normalisations
# RPKM is worse if we compare across samples, only within a sample!!
# what si false positive rate of differential gene expression?
# calculated via false positive rate, 
# false positive rate lowest with Dseq or TMM
# assumption that gene expression isnt that diffent
# how Dseq works? all pairwise comparison, checks median, normalising read count
#anylzses the medians across all samples...
# edgeR pick one reference sample and then compares all other to this
# refernce read comparison, uses weighted mean, normalises differently...

# why use logged data?
# up and down regulation
# log ratios symmetrical around 0, 

getwd()

# SRR36606 => individual
# Ensemmu => gene name?
# values = expression levels

readcounts=read.table("Macaques_readcounts.txt", header=TRUE, row.names=1, sep="\t", stringsAsFactors=FALSE)
dim(readcounts) #check table dimensions

head(readcounts)
# shows 50 samples from 25 different individfual, LPS and control

coldata=read.table("Macaques_ID.txt", header=TRUE, row.names=1, sep="\t", stringsAsFactors=FALSE)
dim(coldata) #check table dimensions

if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("DESeq2")

library("DESeq2")
a

coldata$condition = factor(coldata$condition) # need to set factors to compare them
# now we set DESeq object, can contain multiple tables/matrices
# => can work with mutiple tables
# we call it dds, DESeqDataSet

dds=DESeqDataSetFromMatrix(countData = readcounts, colData = coldata, design = ~ condition)
dds # shows is some info about file
# column names = sample name / individual
# row names = gene name
# coldta names = gives us the condition i.e Control or LPS

# we now want to filter reads
# common is taking out all genes and expressions unde a threshold => noise 
keep = rowSums(counts(dds)) >= 10 #keep rows that have more than 10 in all rows
dds = dds[keep,]
dds # smaller  number of genes now

# specify factors for analysis
dds$condition = factor(dds$condition, levels = c("NC","LPS"))

# now we wrun the normalisation
dds = DESeq(dds) # will get a status report
# per default does the differential gene expression
# will give up p values for the differential expressed genes?
# wil also adjust the p values via  benjamin hochberg method

result = results(dds) 
mcols(result)$description 

## [1] "mean of normalized counts for all samples" => baseMEAN
## [2] "log2 fold change (MLE): condition LPS vs NC" =
## [3] "standard error: condition LPS vs NC"      => standard error  
## [4] "Wald statistic: condition LPS vs NC"        
## [5] "Wald test p-value: condition LPS vs NC"     
## [6] "BH adjusted p-values

# p value adjustment
# p value 0,01 => one in 100 will be false positive
# very bad when you test for 21k genes
# benferroni would destroy pour p value for big data sets
# use benjamini hochberg => compromise between
# we reduce the p values to compensate for potential values

# when analyse big data set
# shrinkage = reduce the noide of dataset

resultLFC = lfcShrink(dds, coef="condition_LPS_vs_NC", type="normal")
resultLFC
# specify that we are comparing LPS and NC

summary(result)
# shows us more info, result here an object (more complicated matrix?)
# could still exclude the low count reads

# sort result via p value

resultOrdered = result[order(result$pvalue),]
sum(result$padj < 0.1, na.rm=TRUE)
sum(result$padj < 0.05, na.rm=TRUE) #the significant ones at 5% level!!

# 5% threshold
# of 1000 we expect 5% to be false positives and not significantly different 
# from the reference
# the up and down regulated ones are the ones that falls under the 5%
# can adjust alpha

#Filtering of genes with False Discovery Rate (FDR) smaller than 0.05 and 
#assigning them to a new object:
result005 = results(dds, alpha=0.05)
summary(result005)

#plotMA() is used for plotting log2 Fold changes vs. mean expression 
#values (baseMean) for all genes. Genes with p<0.1 are shown in red. 
# Points with smaller or larger values than the y-axis are shown as triangles.
plotMA(result, ylim=c(-2,2))

plotMA(resultLFC, ylim=c(-2,2)) 
#shows the shrunken set
# shows left how far up or down regulated
# x aches shows number of reads (many or not?) ??
# can also plot information for indiv. genes
# here take the gene with the lowest p value

plotCounts(dds, gene=which.min(result$padj), intgroup="condition")

#Exporting results:
#Writing a result file with gene names and p values:
write.csv(as.data.frame(resultOrdered), file="DESeq2_DEgenes_condition_LPS_NC.csv")
# saved in working directory
# getwd() setwd()


# Edge R ------------------------------------------------------------------

# both normalise but edgeR compares the reads to one reference
# uses the weighted mean and library size?
library("edgeR")

count_edgeR_obj=DGEList(counts=readcounts, group=coldata$condition)
count_edgeR_obj
# lib size => corrects library size is mean is not 1?

# 2 Normalisation steps
count_edgeR_obj=estimateCommonDisp(count_edgeR_obj) #need to do both
count_edgeR_obj=estimateTagwiseDisp(count_edgeR_obj)
#edgeR uses the quantile-adjusted conditional maximum likelihood (qCML) 
#method to estimate the dispersion(s) before comparing two groups of 
#samples. It first determines the common dispersion and then the 
#dispersion for each gene. For the gene-wise dispersion, it implements 
#an empirical Bayes strategy for squeezing the gene-wise dispersions 
#towards the common dispersion.
#Takes a few seconds!

# Result of normalisation
edgeR_DEgenes=exactTest(count_edgeR_obj)
# negative values are more expressed in controls
# positive values are more expressed in LPS samples
# p values sorted from smallest to largest
# genes called tags here, remnant from times where we counted gene tags not
# genes themselves

#The function topTags() is used to show the top differentially 
#expressed genes (default: based on p-value).
topTags(edgeR_DEgenes)

edgeR_DEgenes
# both edgeR and DESeq2 agree on which genes p value is smallest

#To show the top differentially expressed genes based on fold change use:
topTags(edgeR_DEgenes, sort.by = "logFC")

#As seen above, the edgeR_DEgenes object contains multiple elements. 
#The first one is the table with logFC, logCPM, and p-values for each gene. 
#To get access to this table and assign it to a new variable, call:
edgeR_DEgenesTable=edgeR_DEgenes$table
head(edgeR_DEgenesTable)

#Now you can extract significant genes.
signedgeR_DEgenes=edgeR_DEgenesTable[edgeR_DEgenesTable[,3]<0.05,]

#Write a result file with genes sorted by p-value.
edgeROrdered <- edgeR_DEgenesTable[order(edgeR_DEgenesTable$PValue),]
write.csv(as.data.frame(edgeR_DEgenesTable), file="edgeR_DEgenes_condition_LPS_NC.csv")


# Homework ----------------------------------------------------------------

#Excercises

#Which method calls more significant genes, DESeq2 or edgeR?
  
#How many genes were called differentially expressed with both methods? 
#Hint1: check out %in% Hint2: are p-values adjusted for both methods? 
#Check overlap!!
  
#A venn diagram is also a good way to visualise overlap. 
#Check out the eulerr package (Note, this package is not installed on 
#the evop server, so work with your neighbor who has a PC or Mac)

# can go through rest of tutorial
# how does rank indluence expression?

if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("biomaRt")

# 400 significant genes, what now?
# gene list is the beginning of analysis process
# can now see if genes have a function, associated with different diseases?
# how are they integrated in network?
# are they expressed in other tissues?
# connected pathways? 
# which are specific to cell species,

# Gene ontology Data base, founded 1988
# database about gene function
# controlled vocabulary to help search
# 3 taxonomies of functional categories
# biological process (signal transduction, digestion etc)
# molecular function (transcription factor, receptor, transmembrane transporter)
# cellular component (mitochondria, cytosol, nucleus etc)

#! multiple testing problems

# to attach function to genes must use biomart
# only then can we use topGO to check for enrichment
# biomaRt => can be queried using SQL
# => online data basde than can be intergrated through R?
# give me genes and orthologues, give me etc....

a

library(biomaRt)

listMarts() # shows categories
#We first have to choose a Mart and a Dataset.
#The function listMarts() shows all available Marts. We need ensembl

ensembl=useMart("ensembl") #  makes object

listDatasets(ensembl)
mart = useDataset("mmulatta_gene_ensembl", useMart("ensembl"))
#We want to convert the Ensembl Gene IDs of all expressed genes 
#into GeneSymbols (“external gene names” document).

expressedGenes=row.names(result)
listAttributes(mart)
# listattributes shiws all dfferent,... with gene name
# getbm = get biomart = get information and put into 

GeneSymbols = getBM(filters= "ensembl_gene_id", attributes= c("ensembl_gene_id",
            'external_gene_name'),values=expressedGenes,mart= mart)

# find all the tables we need then use getbm to put everything in one
# file?
# need internet to access this biomaRt database

dim(GeneSymbols) 
# can be very slow if many scientist using it, if doesnt work we can 
# use a different server mirror
# we use ensembl gene ID to find out what the gene is called

head (GeneSymbols)
# boom we now have the names thanks to the IDs


# Getting GO (gene ontology) Information --------------------------------------------------

listMarts()
ensembl=useMart("ensembl")
listDatasets(ensembl)
mart <- useDataset("hsapiens_gene_ensembl", useMart("ensembl"))

head(GeneSymbols) 
# already got the table that shows all expressed 
# genes with ther ID and Gene names

# remove out duplicates
uniqGenesymbols = unique(GeneSymbols[,2])
head(uniqGenesymbols)
GeneGONames = getBM(filters= "hgnc_symbol", attributes= c("hgnc_symbol", "go_id"),
                    values=uniqGenesymbols, mart= mart)
dim()
#attributes is from => got it from listattributes function 

# we habe now attached to GO information to our gene ID, so we can now go
# use topGO o see if they are gene enriched
# need some information for this analysis
# we need a list of all the genese as a background/reference to compare to
# background = expressed genes
# we need all the significant genes
# need the mapping for the genes => biomaRt

library("topGO")
result
DE_genes = subset(row.names(result), result$padj<0.01)
head(GeneSymbols)

DEGenesymbols = subset(GeneSymbols, GeneSymbols$ensembl_gene_id %in% DE_genes)
#%in% checks overlap, which genes are in both?
uniqDEGenesymbols = unique(DEGenesymbols[,2])
length(uniqDEGenesymbols)

geneList = factor(as.integer(uniqGenesymbols %in% uniqDEGenesymbols))
names(geneList) = uniqGenesymbols
str(geneList)
#gives us a vector that has 1 or not => 1 is differentially expressed and 0 is 
# not differentially expressed

# need the mapping of the genes to the GO group in a particular format
GeneGONames[GeneGONames[,1]=="PDPN",1:2]

#Now we reformat the entire table, gene by gene, and save the reformated table 
#into a file. This takes some seconds.

Genes2GO = matrix(,length(uniqGenesymbols),2)
for (i in 1:length(uniqGenesymbols))
{
  temp=GeneGONames[GeneGONames[,1]==uniqGenesymbols[i],1:2]
  tempGOs = paste(temp[,2], collapse=",")
  Genes2GO[i,1]=temp[1,1]
  Genes2GO[i,2]=tempGOs
}
write.table(Genes2GO, "Genes2Go.txt", quote=FALSE, row.names=FALSE, 
            col.names=FALSE, se="\t")
Genes2GOmap=readMappings(file = "Genes2Go.txt")
str(head(Genes2GOmap))

GOdata = new("topGOdata", ontology = "MF", allGenes = geneList,
             annot = annFUN.gene2GO, gene2GO = Genes2GOmap)

#The topGOdata object holds information about the ontology to analyze, 
#number of genes and DE genes. Note that not all genes have GO information, 
#which is why the number of feasible genes is smaller. Also the information 
#about the GO graph is contained in the topGOdata object:
GOdata
# ontology MF => molecular function

# fisher exact test:
resultFisher = runTest(GOdata, algorithm = "classic", statistic = "fisher")
# tested over 4000 GO groups
# 78 terms has a p value under 0,0

resultFisher

sigterms = resultFisher@geneData["SigTerms"] #option within object?
# list is a number of different variables (tables, vectors etv)
# @ is used like a dollar sign to reference a variable but for a list not a table?
sigterms

sigGOIDs = GenTable(GOdata, classicFisher = resultFisher, topNodes = sigterms)
head(sigGOIDs)

# need to adjust the p values for multiple testing = benjamini hochberg
# take the GO groups that are still significant after 
#multiple testing and p adjustment

qval = p.adjust(sigGOIDs$classicFisher, met='BH')
sigGOIDscorrected = cbind(sigGOIDs, qval)
head(sigGOIDscorrected, n=30) # default is 6 rows but can change to 30
#         GO.ID                                        Term Annotated Significant
#1  GO:0005126                   cytokine receptor binding       167          35

#    Expected classicFisher         qval
#1      8.93       2.1e-12 2.797200e-09

library("Rgraphviz")

showSigOfNodes(GOdata, score(resultFisher), firstSigNodes = 5, useInfo = 'all')
# the redder the more significant

# Homework ----------------------------------------------------------------
# do excatly the same however find the GO MAp/ traits that describe the 
# genes function in biological processess

#GOData =.... ontology = MF => change to BP (Biological processes)

library(biomaRt)

listMarts() # shows categories
#We first have to choose a Mart and a Dataset.
#The function listMarts() shows all available Marts. We need ensembl

ensembl=useMart("ensembl") # we select enssembl database and make an object

listDatasets(ensembl) # list what ensembl has to offer = 200+ categories
mart = useDataset("mmulatta_gene_ensembl", useMart("ensembl"))
# we select mmulatta because we are studying macaques
#We want to convert the Ensembl Gene IDs of all expressed genes 
#into GeneSymbols (“external gene names” document).

expressedGenes=row.names(result) 
# need result from previous collating of two documents
listAttributes(mart)
# listattributes shiws all dfferent,... with gene name
# getbm = get biomart = get information and put into 

GeneSymbols = getBM(filters= "ensembl_gene_id", attributes= c("ensembl_gene_id",
                                                              'external_gene_name'),values=expressedGenes,mart= mart)

# find all the tables we need then use getbm to put everything in one
# file?
# need internet to access this biomaRt database

dim(GeneSymbols) 
# can be very slow if many scientist using it, if doesnt work we can 
# use a different server mirror
# we use ensembl gene ID to find out what the gene is called

head (GeneSymbols)
# boom we now have the names thanks to the IDs


# Getting GO (gene ontology) Information --------------------------------------------------

listMarts()
ensembl=useMart("ensembl")
listDatasets(ensembl)
mart <- useDataset("hsapiens_gene_ensembl", useMart("ensembl"))

head(GeneSymbols) 
# already got the table that shows all expressed 
# genes with ther ID and Gene names

# remove out duplicates
uniqGenesymbols = unique(GeneSymbols[,2])
head(uniqGenesymbols)
GeneGONames = getBM(filters= "hgnc_symbol", attributes= c("hgnc_symbol", "go_id"),
                    values=uniqGenesymbols, mart= mart)
dim(GeneGONames)
#attributes is from => got it from listattributes function 

# we habe now attached to GO information to our gene ID, so we can now go
# use topGO o see if they are gene enriched
# need some information for this analysis
# we need a list of all the genese as a background/reference to compare to
# background = expressed genes
# we need all the significant genes
# need the mapping for the genes => biomaRt

library("topGO")
result
DE_genes = subset(row.names(result), result$padj<0.01)
head(GeneSymbols)

DEGenesymbols = subset(GeneSymbols, GeneSymbols$ensembl_gene_id %in% DE_genes)
#%in% checks overlap, which genes are in both?
uniqDEGenesymbols = unique(DEGenesymbols[,2])
length(uniqDEGenesymbols)

geneList = factor(as.integer(uniqGenesymbols %in% uniqDEGenesymbols))
names(geneList) = uniqGenesymbols
str(geneList)
#gives us a vector that has 1 or not => 1 is differentially expressed and 0 is 
# not differentially expressed

# need the mapping of the genes to the GO group in a particular format
GeneGONames[GeneGONames[,1]=="PDPN",1:2]

#Now we reformat the entire table, gene by gene, and save the reformated table 
#into a file. This takes some seconds.

Genes2GO = matrix(,length(uniqGenesymbols),2)
for (i in 1:length(uniqGenesymbols))
{
  temp=GeneGONames[GeneGONames[,1]==uniqGenesymbols[i],1:2]
  tempGOs = paste(temp[,2], collapse=",")
  Genes2GO[i,1]=temp[1,1]
  Genes2GO[i,2]=tempGOs
}
write.table(Genes2GO, "Genes2Go.txt", quote=FALSE, row.names=FALSE, 
            col.names=FALSE, se="\t")
Genes2GOmap=readMappings(file = "Genes2Go.txt")
str(head(Genes2GOmap))

GOdata = new("topGOdata", ontology = "BP", allGenes = geneList,
             annot = annFUN.gene2GO, gene2GO = Genes2GOmap)

#The topGOdata object holds information about the ontology to analyze, 
#number of genes and DE genes. Note that not all genes have GO information, 
#which is why the number of feasible genes is smaller. Also the information 
#about the GO graph is contained in the topGOdata object:
GOdata
# ontology MF => molecular function

# fisher exact test:
resultFisher = runTest(GOdata, algorithm = "classic", statistic = "fisher")
# tested over 4000 GO groups
# 78 terms has a p value under 0,0

resultFisher

sigterms = resultFisher@geneData["SigTerms"] #option within object?
# list is a number of different variables (tables, vectors etv)
# @ is used like a dollar sign to reference a variable but for a list not a table?
sigterms

sigGOIDs = GenTable(GOdata, classicFisher = resultFisher, topNodes = sigterms)
head(sigGOIDs)

# need to adjust the p values for multiple testing = benjamini hochberg
# take the GO groups that are still significant after 
#multiple testing and p adjustment

qval = p.adjust(sigGOIDs$classicFisher, met='BH')
sigGOIDscorrected = cbind(sigGOIDs, qval)
head(sigGOIDscorrected, n=30) # default is 6 rows but can change to 30
#         GO.ID                                        Term Annotated Significant
#1  GO:0005126                   cytokine receptor binding       167          35

#    Expected classicFisher         qval
#1      8.93       2.1e-12 2.797200e-09

library("Rgraphviz")

showSigOfNodes(GOdata, score(resultFisher), firstSigNodes = 5, useInfo = 'all')


# Correction --------------------------------------------------------------

# exercise 1
# Dseq2
# had made a resultordered file => can subset that for the differentially
# expressed genes
# can also use dim on signDESe2_DEgenes
# length[subsetting index 3rd column from signDESe2_DEgenes] 
# => gives Anzahl of columns, ie number of differentially expressed
# ---------------------
# edge R
# use function toptags, the top 3500 genes (n=3500) with p value=0,05
# => makes a subset of first 3500 genes for those with p value under 0,05
# can use function to subset also by using sort.by (cutoff only p value though)
# !!!!! really useful
# => make a new variable with these significant gene values, need this
# to compare DESeq2 and edgeR
# --------------------
# the gene with the lowest p value was the same
#-----------------------------------
# overlap:
# we make a vector with the significantly (underp 0.05) differentiall 
# expressed genes from DESeq2 file and edgeR file => take rownames()
# give us the length of vector %in% vector (i%in% is overlap)
# almost all overlap?
# if want to be sure about false positives, use both DESeq2 and edgeR!!
#-----------
# Venn diagramm with eulerr package
# select two colours
# adjusted edge R has a lot less genes and more overlap = better  version?
# adjusted p values (p values reduced) => corrected for less false positives
# benjamini-hochberg to adjust p values
# chat gpt answer about false discovery rate, false positive rate over gene
#----------
# we initially used the macaque database over biomaRt, however
# the Gene ontology (GO) information we used the human database
# GO data is not available for every species and the human is most studied
# we assume that orthologue in another species has the same function however

# venn diagramm
install.packages("eulerr")
fit = eulerr::euler() #....
#... need to make a significant gene file (pvalue under 0,05)

