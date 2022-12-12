library("phyloseq")
library("ggplot2")
library("vegan")
library("data.table")

setwd("path/to/directory/with/input/and/output/files")

#load your ASV/OTU count table, a taxonomy table and a table with metadata
otu<-as.matrix(read.table("ASV_counts.txt", header=TRUE, row.names = 1))
tax<-as.matrix(read.table("ASV_taxonomy.txt", header=TRUE,row.names = 1, na.strings=c(""," ","NA")))
meta_table<-read.table("metadata.csv",row.names = 1,sep=",",header=TRUE)
meta_table$date<-as.Date(meta_table$date)
meta_table$month<-month.abb[meta_table$month]

#create a phyloseq object from your input tables
OTU<-otu_table(otu, taxa_are_rows = TRUE)
TAX<-tax_table(tax)
SAM<-sample_data(meta_table)
ps <- merge_phyloseq(phyloseq(OTU, TAX), SAM)
rank_names(ps)
ps
View(tax_table(ps)@.Data)

#Remove ASVs/OTUs without taxonomic assignment
ps<-subset_taxa(ps, !is.na(Supergroup))
ps

#Remove metazoans and plants if you're only interested in protists/unicellular eukaryotes
filterPhylum<-c("Metazoa","Streptophyta")
ps<- subset_taxa(ps, !Phylum %in% filterPhylum)
ps

#show metadata affiliate with the samples
sample_variables(ps)
#show read number per sample
sample_sums(ps)

#Check occurrence across samples and taxonomic groups of your ASVs/OTUs
#Create table, number of ASV for each supergroup
table(tax_table(ps)[, "Supergroup"], exclude = NULL)

#Compute prevalence of each ASV (number of samples in which it appears at least once)
prevdf = apply(X = otu_table(ps),
               MARGIN = ifelse(taxa_are_rows(ps), yes = 1, no = 0),
               FUN = function(x){sum(x > 0)})
prevdf
#Add taxonomy and total read counts to this data.frame
prevdf = data.frame(Prevalence = prevdf,
                    TotalAbundance = taxa_sums(ps),
                    tax_table(ps))
prevdf

#Compute the  average and total prevalences of the ASVs in each supergroup across samples
plyr::ddply(prevdf, "Supergroup", function(df1){cbind(mean(df1$Prevalence),sum(df1$Prevalence))})

#Subset to the remaining phyla
prevdf1 = apply(X = otu_table(ps),
                MARGIN = ifelse(taxa_are_rows(ps), yes = 1, no = 0),
                FUN = function(x){sum(x > 0)})

#Add taxonomy and total read counts to this data.frame
prevdf1 = data.frame(Prevalence = prevdf1,
                     TotalAbundance = taxa_sums(ps),
                     tax_table(ps))
prevdf1 = subset(prevdf1, Supergroup %in% get_taxa_unique(ps, "Supergroup"))

#plot total abundance (read abundance) of ASVs in each supergroup vs. prevalence
ggplot(prevdf1, aes(TotalAbundance, Prevalence / nsamples(ps),color=Supergroup)) +
  geom_point(size = 0.8, alpha = 0.7) +
  scale_x_log10() + xlab("Total Abundance") + ylab("Prevalence [Frac. Samples]") +
  facet_wrap(~Supergroup) + theme(legend.position="none")

#calculate median/mean number of reads per ASV
median <- median(prevdf1$TotalAbundance)
median
mean <- mean(prevdf1$TotalAbundance)
mean

#It is advised to remove rare ASVs as they might only represent sequencing errors. There are different approaches depending on how strict you want to be.
#remove ASVs with prevalence of 1/occur only in one sample
ps1 <- filter_taxa(ps, function (x) {sum(x > 0) > 1}, prune=TRUE)
ps1
ps

#remove ASVs that have less reads than the median across all samples
ps1 <- prune_taxa(taxa_sums(ps) >= median, ps)
ps1
ps

#Plot distribution of the ASVs across samples
qplot(log10(colSums(otu_table(ps1)))) +
  xlab("Logged counts-per-sample")
boxplot(colSums(otu_table(ps1)))
stripchart(colSums(otu_table(ps1)), vertical = TRUE, 
           method = "jitter", add = TRUE, pch = 20, col = 'blue')

min(sample_sums(ps1))#
sample_sums(ps1)
which(colSums(otu_table(ps1)) < 5000)


#select only Erken samples
EK <- subset_samples(ps1, lake=="EK")

#remove ASV that don't occurr in any EK samples
EK <- prune_taxa(taxa_sums(EK) > 0 , EK)
EK

#Rarefy samples. This is only one approach to do this and there are many other (maybe better) ways.
EKrar <- rarefy_even_depth(EK, sample.size = min(sample_sums(EK)),
                           rngseed = TRUE, replace = TRUE, trimOTUs = TRUE, verbose = TRUE)
sample_sums(EKrar)

#Plot samples in ordination
NMDS <- ordinate(EKrar, "NMDS", "bray")
plot_ordination(EKrar, NMDS, type="sample", color="month", shape= "fraction", title="NMDS Erken")
plot_ordination(EKrar, NMDS, type="taxa", color="Phylum", shape= "Supergroup", title="NMDS Erken")
plot_ordination(EKrar, NMDS, type="taxa", color="Supergroup", title="NMDS Erken")

PCoA <- ordinate(EKrar, method = "PCoA", distance = "bray")
coords<-PCoA$vectors
evals<-PCoA$values$Eigenvalues
plot_ordination(EKrar, PCoA, color = "month", shape= "fraction", title="PCoA Erken") +
  coord_fixed(sqrt(evals[2] / evals[1])) + 
  theme_classic()

#Bray-Curtis distance matrix
#Sometimes it makes sense to merge the ASV at a specific taxonomic rank (e.g. Genus)
EKrarglomgen = tax_glom(EKrar, "Genus", NArm = FALSE)
BC.dist<-phyloseq::distance(EKrarglomgen, "bray")
BC.dist

#Compare agglomeration methods
par(mfcol = c(1,3)) ## To plot the three clustering trees side-by-side
plot(hclust(BC.dist, method = "average"))
plot(hclust(BC.dist, method = "ward.D2"))#My choice
plot(hclust(BC.dist, method = "single"))
par(mfrow=c(1,1))

#alpha-diversity
plot_richness(EKrar, x="date", color="fraction", measures = c("Observed", "Shannon", "Simpson"))

#barplot
#We merge ASVs at class level and subset the data set to only display the most abundant classes to improve readability
EK_glom = tax_glom(EK, "Class", NArm = FALSE)
EK_glom
SubSet = apply(X = otu_table(EK_glom),
               MARGIN = ifelse(taxa_are_rows(EK_glom), yes = 1, no = 2),
               FUN = function(x){sum(x > 0)})
SubSet2 = data.frame(Prevalence = SubSet,
                     TotalAbundance = taxa_sums(EK_glom))
SubSet2 = data.frame(Prevalence = SubSet,
                     TotalAbundance = taxa_sums(EK_glom),
                     Percent = 100 / sum(SubSet2$TotalAbundance)* SubSet2$TotalAbundance)
keepTaxa = rownames(SubSet2)[SubSet2$Prevalence > 0]
PrunedTaxaEKglom = prune_taxa(keepTaxa, EK_glom)
PrunedTaxaEKglom
keepTaxa2 = rownames(SubSet2)[SubSet2$Percent > 0.5]
PrunedTaxaEKglom2 = prune_taxa(keepTaxa2, PrunedTaxaEKglom)
PrunedTaxaEKglom2
EKglom_rar<-rarefy_even_depth(PrunedTaxaEKglom2, sample.size = min(sample_sums(PrunedTaxaEKglom2)),
                               rngseed = TRUE, replace = TRUE, trimOTUs = TRUE, verbose = TRUE)

phyloseq_barplot <- plot_bar(EKglom_rar, "date", fill="Class", facet_grid=~fraction)
print(phyloseq_barplot)

#heatmap 
plot_heatmap(EKglom_rar, "NMDS", "bray", "sample", "Class")
plot_heatmap(EKglom_rar, method = "NMDS", distance = "bray", sample.label = "sample", taxa.label = "Class", sample.order = "fraction")
