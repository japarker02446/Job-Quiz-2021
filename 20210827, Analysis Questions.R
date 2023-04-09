# 20210827, Analysis Questions.R
# 
# Analysis of data and questions in support of application for
# Principle Scientist, Bioinformatics.
# 
# Data gathered from various sources:
# - immune_genes.txt.   This is a list of immune-related genes.
# -dlbc*: TCGA data processed through Excel (check for date auto correct errors).
#       -dlbc_cna: Copy number abberation.
#       -dlbc_log2cna: Log2(cna + 1).
#       -dlbc_expr: median RNA-seq gene expression.
# - broad_info.csv: Broad Institute Cell line metadata.
# - broad_crispr.txt: Matrix of gene vs Cell Line.  Table values are the probability of negatively 
#       affecting proliferation if the gene is knocked out.
#       - Pre-filtered to genes in immune_genes.txt.
#       - A cell line is "dependent" on a gene if the gene negatively affects proliferation.

# Cleanse your palate.
rm(list = ls());

# Set global variables, read the data.
inputDir <- 'C:/Users/jparker/Code/Input/Other/';
outputDir <- 'C:/Users/jparker/Code/Output/';
dateStamp <- gsub("-","",as.Date(Sys.time()));

#####
# Question 1: Find the genes that are most commonly homozygous deleted.  
#   Name this set of genes homDel.
#   Statistically define "commonly" and justify the definition.
#   
# Thinks: 
# Use the dlbc_cna table to identify the homozygous deleted genes (-2s).
# Common is going to be a matter of homozygous deletion frequency across all samples, that is
# which genes are homozygous deleted in the most samples?.
setwd(inputDir);

# hm, R is throwing an error on reading this table that line 24960 doesn't have 50 elements.
# Not sure if this is intentional or just a me hiccup, but adding "fill = TRUE" to get around it.
cnaTable <- read.table("dlbc_cna.txt", header = TRUE, stringsAsFactors = FALSE, fill = TRUE);

# Did Excel eat the symbols?
head(cnaTable[order(cnaTable$Hugo_Symbol), "Hugo_Symbol"]);
# yes, yes it did.
# Use the Entre_Gene_Id column for reliable and unmunged variables.

# There are fifty columns, first two are gene identifiers (symbol and EntreID) so there are 48 samples.
# How many genes are homozygous deleted (-2) across samples?
cnaTable$homDelCount <- rowSums(cnaTable == -2, na.rm = TRUE);
table(cnaTable$homDelCount);
delHist <- hist(cnaTable$homDelCount);

# OK, the table shows that there are (... carry the one ...) 141 genes that are homozygous deleted 
# across  at least five samples but only 10 genes that are homozygous deleted across at least six
# samples.
# 
# What proportion is in the top?
sum(delHist$density[5:13]); # More than four samples? < 1%
sum(delHist$density[4:13]); # More than three samples? > 1%

# So the most commonly deleted homozygous deleted genes are those deleted in more than four samples,
# or those genes in the top 1% by frequency.
homDel <- cnaTable[which(cnaTable$homDelCount >= 4),];
homDel <- homDel[, which(names(homDel) != "homDelCount")];

# Other thoughts:
#   I had considered using Chi squared or ANOVA testing to determine the statistically significant 
#   set of most frequently homozygous deletion mutated genes, but kept running into walls 
#   remembering / figuring out how if that would work or how to set up the model.

# Cleanup
rm(delHist);

#####
# Question 2: Compare the calls of cna.txt with log2cna.txt and expr.txt for set homDel.
# Are they correlated?
log2Table <- read.table("dlbc_log2cna.txt", header = TRUE, stringsAsFactors = FALSE, fill = TRUE);
exprTable <- read.table("dlbc_expr.txt", header = TRUE, stringsAsFactors = FALSE, fill = TRUE);

# # Interestingly the numbers of genes after filtering each table to the list of Gene IDs in 
# # homDel is NOT the same.
# log2Table <- log2Table[which(log2Table$Entrez_Gene_Id %in% homDel$Entrez_Gene_Id), ];
# exprTable <- exprTable[which(exprTable$Entrez_Gene_Id %in% homDel$Entrez_Gene_Id), ];

# We need to merge the various tables to do a correlation test.
library(data.table);
homDel <- melt(data.table(homDel), id.vars = c("Hugo_Symbol", "Entrez_Gene_Id"), na.rm = TRUE);
log2Table <- melt(data.table(log2Table), id.vars = c("Hugo_Symbol", "Entrez_Gene_Id"), na.rm = TRUE);
exprTable <- melt(data.table(exprTable), id.vars = c("Hugo_Symbol", "Entrez_Gene_Id"), na.rm = TRUE);

colnames(homDel) <- c("Hugo_Symbol", "Entrez_Gene_Id", "Sample", "HomDel");
colnames(log2Table) <- c("Hugo_Symbol", "Entrez_Gene_Id", "Sample", "log2CNA");
colnames(exprTable) <- c("Hugo_Symbol", "Entrez_Gene_Id", "Sample", "meanExpr");

# The merge without using all.x or all.y will only keep values where both the sample
# and gene identifiers match from both tables.
delLog <- merge(homDel, log2Table);
delExpr <- merge(homDel, exprTable);
logExpr <- merge(log2Table, exprTable);

# Our homDel variable is effectively categorical with five levels - ANOVA
# The other two variables are continuous - Pearsons correlation

# plot(delLog$HomDel, delLog$log2CNA);
# plot(delExpr$HomDel, delExpr$meanExpr);
# plot(logExpr$log2CNA, logExpr$meanExpr);

### Test of Homozygous Deletion and Log2(CNA + 1) values by ANOVA.
aov.delLog <- aov(HomDel ~ log2CNA, delLog);
summary(aov.delLog);

# The summary output shows the two variables have a p-value <= 2E-16, so they are significantly
# correlated.

### Test of Homozygous Deletion and mean expression values by ANOVA.
aov.delExpr <- aov(HomDel ~ meanExpr, delExpr);
summary(aov.delExpr);

# The summary output shows the two variables have a p-value <= 7.88E-5, also significantly correlated.

### Test of Log2(CNA) and mean expression by Pearsons two sided correlation test.
cor.test(logExpr$log2CNA, logExpr$meanExpr, alternative = "two.sided", method = "pearson");

# Output shows p-value <= 2.2E-16, so these data are also significantly correlated.
#####
# Question 3
# For the set of patients with at least one gene in set homDel homozygous deleted (patientDel), AND
# for genes in the list immune_genes.txt, identify the genes that are very similarly expressed and 
# very differently expressed across the samples.
# - Filter dlbc_expr.txt to samples in homDel and genes in immune_genes.txt

immuneGenes <- read.table("immune_genes.txt", header = TRUE, stringsAsFactors = FALSE);
patientDel <- unique(homDel[which(homDel$HomDel == -2), "Sample"]);
patientDel <- as.character(patientDel$Sample);
exprTable <- read.table("dlbc_expr.txt", header = TRUE, stringsAsFactors = FALSE, fill = TRUE);

exprTable <- exprTable[which(exprTable$Hugo_Symbol %in% immuneGenes$GENE),
                       which(names(exprTable) %in% c("Hugo_Symbol", "Entrez_Gene_Id", patientDel))
            ];

# Assume FGF13 is expressed in all samples.
# Oddly there are two rows of FGF13...
exprTable <- exprTable[-which(exprTable$Hugo_Symbol == "FGF13" & exprTable$TCGA.FA.A7DS.01 == 0), ];

# To show the genes that are similary or differntly expressed in the filtered data set I am going to 
# go with hieararchical clustering useing PVClust.  This will group the genes by similar expression
# across samples.
# References: https://academic.oup.com/bioinformatics/article/22/12/1540/207339
#             https://github.com/shimo-lab/pvclust
library(pvclust);

# Clustering bootstrap analysis to identify biologically relevant groups of genes or outliers.
# Box groups with greater than 95% AU in red.
plotLabels <- apply(exprTable[, c("Hugo_Symbol", "Entrez_Gene_Id")], 1, paste, collapse = ".");
exprTable <- exprTable[, which(names(exprTable) != c("Hugo_Symbol", "Entrez_Gene_Id"))];
rownames(exprTable) <- plotLabels;
exprTable[is.na(exprTable)] <- 0;

# Remove zero expression rows.
exprTable <- exprTable[which(rowSums(exprTable) > 0), ];

# Cluster.  This takes a while.
pvclust_ma.sig <- pvclust(t(exprTable), method.dist="cor", method.hclust="ward.D", nboot=1000);

# Write clusters and lists to output. The genes that are grouped in similar clusters are those
# with similar expression across samples.  Genes in different groups have different epxression
# patterns.  The "AU" value is the approximately unbiased p-value that the clusters are 
# statistically significant with alpha = 0.95 (as indicated in the code).
# 
# Make the png file particularly large so the details are visible when zoomed in.
setwd(outputDir);
png(
    filename = paste0(dateStamp, " Clustering of HomDel Immune Genes.png"), 
    width = 15360, 
    height = 8640
);
plot(pvclust_ma.sig);
pvrect(pvclust_ma.sig, alpha=0.95, pv = "au", type = "geq");
dev.off();

# Write the significant cluster groups to files as well.
clusters <- pvpick(pvclust_ma.sig, alpha = 0.95);
sink(file = paste0(dateStamp, " HomDel Immune Gene Clusters.txt"));
print(clusters$clusters);
sink();

# Well, nothing was statistically significant at p <= 0.05, but they did group together so I am
# calling it a win.

###
# Question 4
# For the genes in homDel that are also present in broad_crisper.txt determine which cell lines,
# if any, are dependent on at least one gene in homDel.
# 
# I'm not going to do this one because I am well over the two hour suggested time frame and have to
# get back to my life.  I will desribe the process I would use to solve this.
# 
# - First, read the crisper data file.
# - Filter the genes by Entrez_Gene_ID in homDel.
# - From there, I'd probably try to find a reference that discusses the values associated with cell 
#   killing to identify a cutoff, then use that value to filter the genes across cell lines.
# - It would be useful to merge a human readable cell line name from the Broad_info.csv file.

# This was really fun.  I sank more hours into it than you wanted (probably four) but it stretched 
# my brain in ways I haven't used in a while and I really enjoyed it.  Thank you.
