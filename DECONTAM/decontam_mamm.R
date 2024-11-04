#Importation des libraries 
library(phyloseq)
library(ggplot2)
library(decontam)
library(readxl)
library(dplyr)
library(reshape2)
library(tidyr)

#Import des données 
#Data externe 
extdata_mamm <- read.table(file = "extdata_mamm.tsv", sep = "\t", header=T)
#Table OTUs
asv_data_mamm_R1_bhum <- read.table(file = "mammal_run1_curated_ASVtable_full.tsv", header = TRUE, row.names = "ASV")
asv_data_mamm_R2_bhum <- read.table(file = "mammal_run2_curated_ASVtable_full.tsv", header = TRUE, row.names = "ASV")
#Assignation taxo
taxo_mamm <- read.table(file = "taxo_vsearch_mammal.tsv", sep=";", header = TRUE, row.names = "ASV")

#FORMATTAGE DONNÉES
taxo_mamm <- tax_table(as.matrix(taxo_mamm))
extdata_mamm <- na.omit(extdata_mamm)
row.names(extdata_mamm) <- extdata_mamm$Echantillons
extdata_mamm$Echantillons <- NULL
extdata_mamm <- sample_data(extdata_mamm)
colnames(asv_data_mamm_R1_bhum) <- toupper(colnames(asv_data_mamm_R1_bhum))
colnames(asv_data_mamm_R2_bhum) <- toupper(colnames(asv_data_mamm_R2_bhum))
asv_data_mamm_R1_bhum <- otu_table(asv_data_mamm_R1_bhum, taxa_are_rows = TRUE)
asv_data_mamm_R2_bhum <- otu_table(asv_data_mamm_R2_bhum, taxa_are_rows = TRUE)


ps_mamm_R1 <- phyloseq(asv_data_mamm_R1_bhum, taxo_mamm, extdata_mamm)
ps_mamm_R1
ps_mamm_R2 <- phyloseq(asv_data_mamm_R2_bhum, taxo_mamm, extdata_mamm)
ps_mamm_R2

#DECONTAM
#Taille des librairies 
head(sample_data(ps_mamm_R1))
head(sample_data(ps_mamm_R2))
df_mamm_R1 <- as.data.frame(sample_data(ps_mamm_R1)) # Put sample_data into a ggplot-friendly data.frame
df_mamm_R2 <- as.data.frame(sample_data(ps_mamm_R2)) # Put sample_data into a ggplot-friendly data.frame
df_mamm_R1$LibrarySize <- sample_sums(ps_mamm_R1)
df_mamm_R2$LibrarySize <- sample_sums(ps_mamm_R2)
df_mamm_R1 <- df_mamm_R1[order(df_mamm_R1$LibrarySize),]
df_mamm_R2 <- df_mamm_R2[order(df_mamm_R2$LibrarySize),]
df_mamm_R1$Index <- seq(nrow(df_mamm_R1))
df_mamm_R2$Index <- seq(nrow(df_mamm_R2))

ggplot(data=df_mamm_R1, aes(x=Index, y=LibrarySize, color=sample_or_control)) + 
  geom_point()
ggplot(data=df_mamm_R2, aes(x=Index, y=LibrarySize, color=sample_or_control)) + 
  geom_point()

#Identification des contaminants - fréquences
sample_data(ps_mamm_R1)$Qubit_ng_mL <- as.numeric(as.character(sample_data(ps_mamm_R1)$Qubit_ng_mL))
sample_data(ps_mamm_R2)$Qubit_ng_mL <- as.numeric(as.character(sample_data(ps_mamm_R2)$Qubit_ng_mL))
ps_mamm_R1 <- prune_samples(!is.na(sample_data(ps_mamm_R1)$Qubit_ng_mL) & sample_data(ps_mamm_R1)$Qubit_ng_mL > 0, ps_mamm_R1)
ps_mamm_R2 <- prune_samples(!is.na(sample_data(ps_mamm_R2)$Qubit_ng_mL) & sample_data(ps_mamm_R2)$Qubit_ng_mL > 0, ps_mamm_R2)
contamdf.freq_mamm_R1 <- isContaminant(ps_mamm_R1, method="frequency", conc="Qubit_ng_mL")
contamdf.freq_mamm_R2 <- isContaminant(ps_mamm_R2, method="frequency", conc="Qubit_ng_mL")
head(contamdf.freq_mamm_R1)
head(contamdf.freq_mamm_R2)
table(contamdf.freq_mamm_R1$contaminant)
table(contamdf.freq_mamm_R2$contaminant)
head(which(contamdf.freq_mamm_R1$contaminant))
head(which(contamdf.freq_mamm_R2$contaminant))

ps_filtered_mamm_R1 <- prune_samples(sample_sums(ps_mamm_R1) > 0, ps_mamm_R1)
ps_filtered_mamm_R2 <- prune_samples(sample_sums(ps_mamm_R2) > 0, ps_mamm_R2)
length(sample_names(ps_filtered_mamm_R1))
length(sample_names(ps_filtered_mamm_R2))

plot_frequency(ps_filtered_mamm_R1, taxa_names(ps_filtered_mamm_R1)[c(1,3)], conc="Qubit_ng_mL") + 
  xlab("DNA Concentration (Qubit)")
plot_frequency(ps_filtered_mamm_R2, taxa_names(ps_filtered_mamm_R2)[c(1,3)], conc="Qubit_ng_mL") + 
  xlab("DNA Concentration (Qubit)")

table(contamdf.freq_mamm_R1$contaminant) #que des FALSE = pas de contaminants détéctés 
set.seed(100)
plot_frequency(ps_filtered_mamm_R1, taxa_names(ps_filtered_mamm_R1)[sample(which(contamdf.freq_mamm_R1$contaminant),3)], conc="Qubit_ng_mL") +
  xlab("DNA Concentration (Qubit)") #plot impossible car pas de contaminants
table(contamdf.freq_mamm_R2$contaminant) #que des FALSE = pas de contaminants détéctés 
set.seed(100)
plot_frequency(ps_filtered_mamm_R2, taxa_names(ps_filtered_mamm_R2)[sample(which(contamdf.freq_mamm_R2$contaminant),3)], conc="Qubit_ng_mL") +
  xlab("DNA Concentration (Qubit)") #plot impossible car pas de contaminants 

#export tables filtrés par frequence 
ps_filtered_mamm_R1
ps.noncontam_mamm_R1 <- prune_taxa(!contamdf.freq_mamm_R1$contaminant, ps_filtered_mamm_R1)
ps.noncontam_mamm_R1
otu_df_R1_mamm <- as.data.frame(otu_table(ps.noncontam_mamm_R1))
write.csv(otu_df_R1_mamm, "otu_table_mamm_R1_freq.csv")
tax_df_R1_mamm <- as.data.frame(tax_table(ps.noncontam_mamm_R1))
write.csv(tax_df_R1_mamm, "taxonomy_table_mamm_R1_freq.csv")

ps_filtered_mamm_R2
ps.noncontam_mamm_R2 <- prune_taxa(!contamdf.freq_mamm_R2$contaminant, ps_filtered_mamm_R2)
ps.noncontam_mamm_R2
otu_df_R2 <- as.data.frame(otu_table(ps.noncontam_mamm_R2))
write.csv(otu_df_R2, "otu_table_mamm_R2_freq.csv")
tax_df_R2 <- as.data.frame(tax_table(ps.noncontam_mamm_R2))
write.csv(tax_df_R2, "taxonomy_table_mamm_R2_freq.csv")

#Joindre table OTU + TAXO   
merge_otu_taxo_mammR1freq <- cbind(otu_df_R1_mamm, tax_df_R1_mamm)
write.csv(merge_otu_taxo_mammR1freq, "merge_otu_taxo_mammR1freq.csv")
merge_otu_taxo_mammR2freq <- cbind(otu_df_R2, tax_df_R2)
write.csv(merge_otu_taxo_mammR2freq, "merge_otu_taxo_mammR2freq.csv")

#Identification des contaminants - prévalance 
sample_data(ps_filtered_mamm_R1)$is.neg <- sample_data(ps_filtered_mamm_R1)$sample_or_control == "Control"
sample_data(ps_filtered_mamm_R2)$is.neg <- sample_data(ps_filtered_mamm_R2)$sample_or_control == "Control"
contamdf.prev_mamm_R1 <- isContaminant(ps_filtered_mamm_R1, method="prevalence", neg="is.neg")
contamdf.prev_mamm_R2 <- isContaminant(ps_filtered_mamm_R2, method="prevalence", neg="is.neg")
table(contamdf.prev_mamm_R1$contaminant)
table(contamdf.prev_mamm_R2$contaminant)
head(which(contamdf.prev_mamm_R1$contaminant))
head(which(contamdf.prev_mamm_R2$contaminant))
contamdf.prev05_mamm_R1 <- isContaminant(ps_filtered_mamm_R1, method="prevalence", neg="is.neg", threshold=0.5)
contamdf.prev05_mamm_R2 <- isContaminant(ps_filtered_mamm_R2, method="prevalence", neg="is.neg", threshold=0.5)
table(contamdf.prev05_mamm_R1$contaminant)
table(contamdf.prev05_mamm_R2$contaminant)

#export tables filtrés par prevalance 
ps.noncontam_mamm_R1_prev <- prune_taxa(!contamdf.prev05_mamm_R1$contaminant, ps_filtered_mamm_R1)
ps.noncontam_mamm_R2_prev <- prune_taxa(!contamdf.prev05_mamm_R2$contaminant, ps_filtered_mamm_R2)

otu_df_mamm_R1_prev <- as.data.frame(otu_table(ps.noncontam_mamm_R1_prev))
write.csv(otu_df_mamm_R1_prev, "otu_table_mamm_R1_prevalence.csv")
tax_df_mamm_R1_prev <- as.data.frame(tax_table(ps.noncontam_mamm_R1_prev))
write.csv(tax_df_mamm_R1_prev, "taxonomy_table_mamm_R1_prevalence.csv")

otu_df_mamm_R2_prev <- as.data.frame(otu_table(ps.noncontam_mamm_R2_prev))
write.csv(otu_df_mamm_R2_prev, "otu_table_mamm_R2_prevalence.csv")
tax_df_mamm_R2_prev <- as.data.frame(tax_table(ps.noncontam_mamm_R2_prev))
write.csv(tax_df_mamm_R2_prev, "taxonomy_table_mamm_R2_prevalence.csv")

#Joindre table OTU + TAXO   
merge_otu_taxo_mammR1prev <- cbind(otu_df_mamm_R1_prev, tax_df_mamm_R1_prev)
write.csv(merge_otu_taxo_mammR1prev, "merge_otu_taxo_mammR1prev.csv")
merge_otu_taxo_mammR2prev <- cbind(otu_df_mamm_R1_prev, tax_df_mamm_R1_prev)
write.csv(merge_otu_taxo_mammR2prev, "merge_otu_taxo_mammR2prev.csv")

# Make phyloseq object of presence-absence in negative controls and true samples
ps.pa_mamm_R1 <- transform_sample_counts(ps_filtered_mamm_R1, function(abund) 1*(abund>0))
ps.pa_mamm_R2 <- transform_sample_counts(ps_filtered_mamm_R2, function(abund) 1*(abund>0))
ps.pa.neg_mamm_R1 <- prune_samples(sample_data(ps.pa_mamm_R1)$sample_or_control == "Control", ps.pa_mamm_R1)
ps.pa.neg_mamm_R2 <- prune_samples(sample_data(ps.pa_mamm_R2)$sample_or_control == "Control", ps.pa_mamm_R2)
ps.pa.pos_mamm_R1 <- prune_samples(sample_data(ps.pa_mamm_R1)$sample_or_control == "Sample", ps.pa_mamm_R1)
ps.pa.pos_mamm_R2 <- prune_samples(sample_data(ps.pa_mamm_R2)$sample_or_control == "Sample", ps.pa_mamm_R2)

# Make data.frame of prevalence in positive and negative samples
df.pa_mamm_R1 <- data.frame(pa.pos_R1=taxa_sums(ps.pa.pos_mamm_R1), pa.neg_R1=taxa_sums(ps.pa.neg_mamm_R1),
                    contaminant=contamdf.prev_mamm_R1$contaminant)
ggplot(data=df.pa_mamm_R1, aes(x=pa.neg_R1, y=pa.pos_R1, color=contaminant)) + geom_point() +
  xlab("Prevalence (Negative Controls)") + ylab("Prevalence (True Samples)")

df.pa_mamm_R2 <- data.frame(pa.pos_R2=taxa_sums(ps.pa.pos_mamm_R2), pa.neg_R2=taxa_sums(ps.pa.neg_mamm_R2),
                       contaminant=contamdf.prev_mamm_R2$contaminant)
ggplot(data=df.pa_mamm_R2, aes(x=pa.neg_R2, y=pa.pos_R2, color=contaminant)) + geom_point() +
  xlab("Prevalence (Negative Controls)") + ylab("Prevalence (True Samples)")

