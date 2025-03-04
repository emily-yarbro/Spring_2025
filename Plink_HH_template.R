
# Intro to PLINK ----------------------------------------------------------

#PLINK is an open source toolset for genome association analyses. What makes it so useful is its' ability to run large-scale analyses in a computationally efficient way.

#We are going to be using PLINK to run a GWAS. In order to do so, we'll need to make sure we have the correct files. I have provided you with three files: "newdgrp2.bed", "newdgrp2.bim", and newdgrp2.fam". These make up the binary fileset that we are going to tell PLINK to reference

#Normally, we would run PLINK in our terminal. However, because some of the other work we are going to do is more R-related, I'm going to have you use the system() function to access PLINK from RStudio:

# PCA --------------------------------------------------------------------

#One of the first things we'll need to consider before running our GWAS is how our covariates come into play. To do this, we'll run a PCA. Thankfully, it's really easy to do this in PLINK.

#--bflie: causes the binary fileset (in this case "newdgrp2" to be referenced)
#--pca: extracts top 20 PCs of variance-standardized relationship matrix
# var-wts: generates an eigenvec file that we can use to help us visualize relationships
#--out: Indicates what we want our output file to be called

#Let's put it all together and run our PCA analysis
system("./plink --bfile newdgrp2 --pca var-wts --out DGRP2_pca")

#If you take a look at your data folder, you should have three new files. One is a log file that contains the same information that popped up in your console when you ran the PCA. The other two are .eigenval and .eigenvec files that we will use to visualize trends in our data.

#Let's start with reading in our eigenvector file and naming it 'pcs'


#The first two rows are line information, so we don't want to include this when we are visualizing our data. Try changing colnames of the first 5 PCs for the plot we are going to make.


#I also think it's fun to add a bit of whimsy to my plots, so let's set up a list of colors to use.


#ploting everything
plot(pcs[,3:7], col = cols, lower.panel = NULL)
#lower.panel = NULL just keeps it from adding redundant panels on the bottom.

#I sometimes find these panels a little difficult to interpret/not as informative about some of the other PC values I should consider including in my GWAS. Let's take a look at our eigenval file to see if we can come up with something better.


#Now we'll read in the eigenvalues and save it as an object 'eigen'


#Let's save the first column of eigen as a list to help with plotting later on

#Create a string of 1-20 for the x-axis of the plot


#Use plot() to generate a scree plot (plot that uses eigenvalues of factors/PCs in a multivariate analysis)



#As we increase the number of clusters, we decrease the variance. We can use our scree plot to visualize where our largest reduction is to inform what an ideal K would be based on the data. In this case, I chose 5.

#saving the PCs so I can add it to my covariate file
write.table(pcs, "DGRP2_pcs_both_sexes.txt", col.names = F, quote = F, row.names = F, sep = "\t")

#Another thing we need to consider involves a little biology. Many of the DGRP lines are infected with Wolbachia, a bacterial parasite whose presence has been shown to impact reproduction, survival, and development (among other things), within Drosophila. Because of this, we need to include infection status in our covariate file. I did this for you already by adding a column to the file "DGRP2_covariates" that I shared with you. In this case, '1' = infected and '0' = non-infected.

#Let's take a look at the layout of the file. We won't save the file as an object, this will be purely visual.
read.table("DGRP2_covariates.txt")


# GWAS (Finally!) ---------------------------------------------------------

#What is a GWAS? It's an association test we can use to map the genetic basis to complex traits.

#Before we get into running the GWAS, let's get familiar with the syntax of our command.

#--bflie: causes the binary fileset (in this case "newdgrp2" to be referenced)
#--pheno: causes phenotype values to be read from third column from specified space/tab-delimited file we've provided (in this case "DT_SD_pheno.txt)
#--pheno-merge: tells PLINK to use phenotype values in .fam file if no value is present in --pheno file. If we don't include this, the phenotype will be considered 'missing'
#--linear: writes a linear regression report to "filename".assoc.linear
#hide-covar: removes covariate-specific lines from main report
#--covar: designates file from which to load the covariates (in this case "DGRP2_covariates.txt)
#--covar-number 1-6: Let's us tell PLINK we only want to focus on values for the first 6 covariates. (NOTE: it ignores the FID and IID from the first two columns)
#--adjust: causes an .adjusted file to be generated with each report that holds multiple testing corrections for the raw p-values
#--out: Indicates what we want our output file to be called

# SD across all concentrations

system("./plink --bfile newdgrp2 --pheno DT_SD_pheno.txt --all-pheno --pheno-merge --linear hide-covar --covar DGRP2_covariates.txt --covar-number 1-6 --adjust --out DT_SD_GWAS")


# Manhattan Plots ---------------------------------------------------------

library(tidyr)
library(ggplot2)


#Reading in the results of our GWAS. This is really really big and may take a long time to load
results<-read.table("DT_SD_GWAS.P1.assoc.linear", header = TRUE, stringsAsFactors = FALSE)

#Let's take a look at the str() of 'results'


#Notice how the "SNP" column holds chromosome, position, and 'type' information. We are going to want to separate these values into different columns within 'results'.

#This will take a little while to run
results<-separate(results,SNP,c("chrom","pos","type"),sep = "_")
#there should now be 11 variables

#convert pos to numeric



#convert to -log(P) so that smaller P values are now higher on the plot


#store results$pos in a new object called 'POS'


#store results$log in a new object called 'logp'


#convert chrom to a factor and store in new object called 'chm'


#Now let's use ggplot() to plot it! We are going to store the plot as an object 'g', and we want to plot results with POS on the x-axis, logp on the y-axis, and we want to color by chromosome. We also want to label the y-axis "-log10 (p-value)", the x-axis "Base Pair Position", and the title "DT SD Manhattan plot".
#Note: I'm asking you to save this as an object, 'g', so we can play around with the orientation of the plot!

#Horizontal stacking for chromosomes:

#NOTE, this takes a very long time to run. You will also get a warning message about removing rows due to missing values - don't panic, that's fine!

#Vertical stacking for chromosomes



# Annotations -------------------------------------------------------------

#After we run our GWAS, we can take the GWAS results with the annotation files and figures out which SNPs meet our p-value cutoff and add annotations to those

#########################
#Naming the input file. Here, we are assigning the character string "DT_SD_GWAS.P1.assoc.linear" to the object "input file". This stores the name of the input file and allows us to call it for later processing.
input.file<-"DT_SD_GWAS.P1.assoc.linear"
#Naming our eventual output file
output.file<-"DT_SD_annot.csv"
###################################

#installing "data.table". This package is essentially an extension of "data.frame" but it allows us to do faster aggregations with large data.
install.packages("data.table")
library(data.table)

#setting our p-value cutoff
Pcut<-0.00001

#Loading the SNPs from our GWAS
Pvalues<-read.table(input.file, header = TRUE)

#load the annotation file. NOTE: don't try to open this one on your computer (it won't like it). This will take a while
annot<-read.table('dgrp.fb557.annot.txt')

#Checking the structure of both data frames
str(Pvalues)
str(annot)

#subset to only those with P<0.00001
sig.snps<-subset(Pvalues,Pvalues$P<Pcut)

#add annot column and other column to sig.snps
sig.snps$annotation<-NA
sig.snps$other<-NA

#Here we are mapping annotations from annot$V3 to corresponding SNPS in sig.snps$SNP based on matching SNP IDs. 
#match(sig.snps$SNP,annot$V1): finding the position of the values of sig.snps$SNP within annot$V1. If a value exists that meets that criteria, match() returns the index of that match. If no match, we get NA.
sig.snps$annotation<-annot$V3[match(sig.snps$SNP,annot$V1)]

#writing our output file!
write.csv(sig.snps,output.file)