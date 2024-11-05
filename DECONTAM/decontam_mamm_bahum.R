#Importation des libraries 
library(phyloseq)
library(ggplot2)
library(decontam)
library(readxl)
library(dplyr)
library(reshape2)
library(tidyr)

#Importation des données 
#Data externe 
extdata <- read.table(file = "extdata.tsv", sep = "\t", header=T)
#Tbale OTUs
asv_data_mamm_R1_bhum <- read.table(file = "mammal_run1_curated_ASVtable_full.tsv", header = TRUE, row.names = "ASV")
asv_data_mamm_R2_bhum <- read.table(file = "mammal_run2_curated_ASVtable_full.tsv", header = TRUE, row.names = "ASV")
#Assignation taxonomique 
taxo <- read.table(file = "taxo_vsearch_mammal.tsv", sep=";", header = TRUE, row.names = "ASV")

#FORMATTAGE DONNÉES
taxo <- taxo[apply(taxo, 1, function(row) any(row != "" & !is.na(row))), ]
taxo <- tax_table(as.matrix(taxo))
extdata <- na.omit(extdata)
row.names(extdata) <- extdata$Echantillons
extdata$Echantillons <- NULL
extdata <- sample_data(extdata)
colnames(asv_data_mamm_R1_bhum) <- toupper(colnames(asv_data_mamm_R1_bhum))
colnames(asv_data_mamm_R2_bhum) <- toupper(colnames(asv_data_mamm_R2_bhum))
asv_data_mamm_R1_bhum <- otu_table(asv_data_mamm_R1_bhum, taxa_are_rows = TRUE)
asv_data_mamm_R2_bhum <- otu_table(asv_data_mamm_R2_bhum, taxa_are_rows = TRUE)

taxo_R1 <- taxo[intersect(row.names(taxo), row.names(asv_data_mamm_R1_bhum)), ]
taxo_R2 <- taxo[intersect(row.names(taxo), row.names(asv_data_mamm_R2_bhum)), ]
asv_data_mamm_R1_bhum <- asv_data_mamm_R1_bhum[intersect(row.names(taxo_R1), row.names(asv_data_mamm_R1_bhum)), ]
asv_data_mamm_R2_bhum <- asv_data_mamm_R2_bhum[intersect(row.names(taxo_R2), row.names(asv_data_mamm_R1_bhum)), ]


ps_R1 <- phyloseq(asv_data_mamm_R1_bhum, taxo_R1, extdata)
ps_R1
ps_R2 <- phyloseq(asv_data_mamm_R2_bhum, taxo_R2, extdata)
ps_R2

#DECONTAM
head(sample_data(ps_R1))
head(sample_data(ps_R2))
df_R1 <- as.data.frame(sample_data(ps_R1)) # Put sample_data into a ggplot-friendly data.frame
df_R2 <- as.data.frame(sample_data(ps_R2)) # Put sample_data into a ggplot-friendly data.frame
df_R1$LibrarySize <- sample_sums(ps_R1)
df_R2$LibrarySize <- sample_sums(ps_R2)
df_R1 <- df_R1[order(df_R1$LibrarySize),]
df_R2 <- df_R2[order(df_R2$LibrarySize),]
df_R1$Index <- seq(nrow(df_R1))
df_R2$Index <- seq(nrow(df_R2))

ggplot(data=df_R1, aes(x=Index, y=LibrarySize, color=sample_or_control)) + 
  geom_point()
ggplot(data=df_R2, aes(x=Index, y=LibrarySize, color=sample_or_control)) + 
  geom_point()

sample_data(ps_R1)$Qubit_ng_mL <- as.numeric(as.character(sample_data(ps_R1)$Qubit_ng_mL))
sample_data(ps_R2)$Qubit_ng_mL <- as.numeric(as.character(sample_data(ps_R2)$Qubit_ng_mL))
ps_R1 <- prune_samples(!is.na(sample_data(ps_R1)$Qubit_ng_mL) & sample_data(ps_R1)$Qubit_ng_mL > 0, ps_R1)
ps_R2 <- prune_samples(!is.na(sample_data(ps_R2)$Qubit_ng_mL) & sample_data(ps_R2)$Qubit_ng_mL > 0, ps_R2)
contamdf.freq_R1 <- isContaminant(ps_R1, method="frequency", conc="Qubit_ng_mL")
contamdf.freq_R2 <- isContaminant(ps_R2, method="frequency", conc="Qubit_ng_mL")
head(contamdf.freq_R1)
head(contamdf.freq_R2)
table(contamdf.freq_R1$contaminant)
table(contamdf.freq_R2$contaminant)
head(which(contamdf.freq_R1$contaminant))
head(which(contamdf.freq_R2$contaminant))

ps_filtered_R1 <- prune_samples(sample_sums(ps_R1) > 0, ps_R1)
ps_filtered_R2 <- prune_samples(sample_sums(ps_R2) > 0, ps_R2)
length(sample_names(ps_filtered_R1))
length(sample_names(ps_filtered_R2))

plot_frequency(ps_filtered_R1, taxa_names(ps_filtered_R1), conc="Qubit_ng_mL") + 
  xlab("DNA Concentration (Qubit)")
plot_frequency(ps_filtered_R2, taxa_names(ps_filtered_R2), conc="Qubit_ng_mL") + 
  xlab("DNA Concentration (Qubit)")

  
set.seed(100)
plot_frequency(ps_filtered_R1, taxa_names(ps_filtered_R1)[sample(which(contamdf.freq_R1$contaminant),3)], conc="Qubit_ng_mL") +
  xlab("DNA Concentration (Qubit)")
set.seed(100)
plot_frequency(ps_filtered_R2, taxa_names(ps_filtered_R2)[sample(which(contamdf.freq_R2$contaminant),3)], conc="Qubit_ng_mL") +
  xlab("DNA Concentration (Qubit)")

ps_filtered_R1
ps.noncontam_R1 <- prune_taxa(!contamdf.freq_R1$contaminant, ps_filtered_R1)
ps.noncontam_R1
ps_filtered_R2
ps.noncontam_R2 <- prune_taxa(!contamdf.freq_R2$contaminant, ps_filtered_R2)
ps.noncontam_R2
sample_data(ps_filtered_R1)$is.neg <- sample_data(ps_filtered_R1)$sample_or_control == "Control"
sample_data(ps_filtered_R2)$is.neg <- sample_data(ps_filtered_R2)$sample_or_control == "Control"
contamdf.prev_R1 <- isContaminant(ps_filtered_R1, method="prevalence", neg="is.neg")
contamdf.prev_R2 <- isContaminant(ps_filtered_R2, method="prevalence", neg="is.neg")
table(contamdf.prev_R1$contaminant)
table(contamdf.prev_R2$contaminant)
head(which(contamdf.prev_R1$contaminant))
head(which(contamdf.prev_R2$contaminant))
contamdf.prev05_R1 <- isContaminant(ps_filtered_R1, method="prevalence", neg="is.neg", threshold=0.5)
contamdf.prev05_R2 <- isContaminant(ps_filtered_R2, method="prevalence", neg="is.neg", threshold=0.5)
table(contamdf.prev05_R1$contaminant)
table(contamdf.prev05_R2$contaminant)

# Make phyloseq object of presence-absence in negative controls and true samples
ps.pa_R1 <- transform_sample_counts(ps_filtered_R1, function(abund) 1*(abund>0))
ps.pa_R2 <- transform_sample_counts(ps_filtered_R2, function(abund) 1*(abund>0))
ps.pa.neg_R1 <- prune_samples(sample_data(ps.pa_R1)$sample_or_control == "Control", ps.pa_R1)
ps.pa.neg_R2 <- prune_samples(sample_data(ps.pa_R2)$sample_or_control == "Control", ps.pa_R2)
ps.pa.pos_R1 <- prune_samples(sample_data(ps.pa_R1)$sample_or_control == "Sample", ps.pa_R1)
ps.pa.pos_R2 <- prune_samples(sample_data(ps.pa_R2)$sample_or_control == "Sample", ps.pa_R2)

# Make data.frame of prevalence in positive and negative samples
df.pa_R1 <- data.frame(pa.pos_R1=taxa_sums(ps.pa.pos_R1), pa.neg_R1=taxa_sums(ps.pa.neg_R1),
                    contaminant=contamdf.prev_R1$contaminant)
ggplot(data=df.pa_R1, aes(x=pa.neg_R1, y=pa.pos_R1, color=contaminant)) + geom_point() +
  xlab("Prevalence (Negative Controls)") + ylab("Prevalence (True Samples)")
df.pa_R2 <- data.frame(pa.pos_R2=taxa_sums(ps.pa.pos_R2), pa.neg_R2=taxa_sums(ps.pa.neg_R2),
                       contaminant=contamdf.prev_R2$contaminant)
ggplot(data=df.pa_R2, aes(x=pa.neg_R2, y=pa.pos_R2, color=contaminant)) + geom_point() +
  xlab("Prevalence (Negative Controls)") + ylab("Prevalence (True Samples)")

