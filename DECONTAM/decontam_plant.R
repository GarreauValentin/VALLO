#LIBRARY
library(phyloseq)
library(ggplot2)
library(decontam)
library(readxl)
library(dplyr)
library(reshape2)
library(tidyr)

#Toutes les commandes pour fichier run 1 sont # car pas de fichier run 1

#Import des données 
#Data externe 
extdata_plant <- read.table(file = "extdata_plant.tsv", sep = "\t", header=T)
#Table OTUs
#asv_data_plant_R1 <- read.table(file = "", header = TRUE, row.names = "ASV")
asv_data_plant_R2 <- read.table(file = "plant_run2_curated_otutable.tsv", header = TRUE, row.names = "ASV")
#Assignation taxo_plantnomique 
taxo_plant_global <- read.table(file = "taxo_vsearch_plant_global.tsv", sep=";", header = TRUE, row.names = "ASV")
taxo_plant_local <- read.table(file = "taxo_vsearch_plant_local.tsv", sep=";", header = TRUE, row.names = "ASV")

#FORMATTAGE DONNÉES
taxo_plant_global <- tax_table(as.matrix(taxo_plant_global))
taxo_plant_local <- tax_table(as.matrix(taxo_plant_local))
extdata_plant <- na.omit(extdata_plant)
row.names(extdata_plant) <- extdata_plant$Echantillons
extdata_plant$Echantillons <- NULL
extdata_plant <- sample_data(extdata_plant)
#colnames(asv_data_plant_R1) <- toupper(colnames(asv_data_plant_R1))
colnames(asv_data_plant_R2) <- toupper(colnames(asv_data_plant_R2))
#asv_data_plant_R1 <- otu_table(asv_data_plant_R1, taxa_are_rows = TRUE)
asv_data_plant_R2 <- otu_table(asv_data_plant_R2, taxa_are_rows = TRUE)


#ps_plant_global_R1 <- phyloseq(asv_data_plant_R1, taxo_plant_global, extdata_plant)
#ps_plant_global_R1
#ps_plant_local_R2 <- phyloseq(asv_data_plant_R2, taxo_plant_local, extdata_plant)
#ps_plant_local_R2
ps_plant_global_R2 <- phyloseq(asv_data_plant_R2, taxo_plant_global, extdata_plant)
ps_plant_global_R2
ps_plant_local_R2 <- phyloseq(asv_data_plant_R2, taxo_plant_local, extdata_plant)
ps_plant_local_R2

#DECONTAM TAXO GLOBAL 
#Taille des librairies 
#head(sample_data(ps_plant_global_R1))
head(sample_data(ps_plant_global_R2))
#df_plant_global_R1 <- as.data.frame(sample_data(ps_plant_global_R1)) # Put sample_data into a ggplot-friendly data.frame
df_plant_global_R2 <- as.data.frame(sample_data(ps_plant_global_R2)) # Put sample_data into a ggplot-friendly data.frame
#df_plant_global_R1$LibrarySize <- sample_sums(ps_plant_global_R1)
df_plant_global_R2$LibrarySize <- sample_sums(ps_plant_global_R2)
#df_plant_global_R1 <- df_plant_global_R1[order(df_plant_global_R1$LibrarySize),]
df_plant_global_R2 <- df_plant_global_R2[order(df_plant_global_R2$LibrarySize),]
#df_plant_global_R1$Index <- seq(nrow(df_plant_global_R1))
df_plant_global_R2$Index <- seq(nrow(df_plant_global_R2))

#ggplot(data=df_plant_global_R1, aes(x=Index, y=LibrarySize, color=sample_or_control)) + 
#  geom_point()
ggplot(data=df_plant_global_R2, aes(x=Index, y=LibrarySize, color=sample_or_control)) + 
  geom_point()

#Identification des contaminants - fréquences
#sample_data(ps_plant_global_R1)$Qubit_ng_mL <- as.numeric(as.character(sample_data(ps_plant_global_R1)$Qubit_ng_mL))
sample_data(ps_plant_global_R2)$Qubit_ng_mL <- as.numeric(as.character(sample_data(ps_plant_global_R2)$Qubit_ng_mL))
#ps_plant_global_R1 <- prune_samples(!is.na(sample_data(ps_plant_global_R1)$Qubit_ng_mL) & sample_data(ps_plant_global_R1)$Qubit_ng_mL > 0, ps_plant_global_R1)
ps_plant_global_R2 <- prune_samples(!is.na(sample_data(ps_plant_global_R2)$Qubit_ng_mL) & sample_data(ps_plant_global_R2)$Qubit_ng_mL > 0, ps_plant_global_R2)
#contamdf.freq_plant_global_R1 <- isContaminant(ps_plant_global_R1, method="frequency", conc="Qubit_ng_mL")
contamdf.freq_plant_global_R2 <- isContaminant(ps_plant_global_R2, method="frequency", conc="Qubit_ng_mL")
#head(contamdf.freq_plant_global_R1)
head(contamdf.freq_plant_global_R2)
#table(contamdf.freq_plant_global_R1$contaminant)
table(contamdf.freq_plant_global_R2$contaminant)
#head(which(contamdf.freq_plant_global_R1$contaminant))
head(which(contamdf.freq_plant_global_R2$contaminant))

#ps_filtered_plant_global_R1 <- prune_samples(sample_sums(ps_plant_global_R1) > 0, ps_plant_global_R1)
ps_filtered_plant_global_R2 <- prune_samples(sample_sums(ps_plant_global_R2) > 0, ps_plant_global_R2)
#length(sample_names(ps_filtered_plant_global_R1))
length(sample_names(ps_filtered_plant_global_R2))

#plot_frequency(ps_filtered_plant_global_R1, taxa_names(ps_filtered_plant_global_R1), conc="Qubit_ng_mL") + 
#  xlab("DNA Concentration (Qubit)")
plot_frequency(ps_filtered_plant_global_R2, taxa_names(ps_filtered_plant_global_R2), conc="Qubit_ng_mL") + 
  xlab("DNA Concentration (Qubit)")

#table(contamdf.freq_plant_global_R1$contaminant) #que des FALSE = pas de contaminants détéctés 
#set.seed(100)
#plot_frequency(ps_filtered_plant_global_R1, taxa_names(ps_filtered_plant_global_R1)[sample(which(contamdf.freq_plant_global_R1$contaminant),3)], conc="Qubit_ng_mL") +
#  xlab("DNA Concentration (Qubit)")
table(contamdf.freq_plant_global_R2$contaminant) #que des FALSE = pas de contaminants détéctés 
set.seed(100)
plot_frequency(ps_filtered_plant_global_R2, taxa_names(ps_filtered_plant_global_R2)[sample(which(contamdf.freq_plant_global_R2$contaminant),3)], conc="Qubit_ng_mL") +
  xlab("DNA Concentration (Qubit)")

#export tables filtrés par frequence 
#ps_filtered_plant_global_R1
#ps.noncontam_plant_global_R1 <- prune_taxa(!contamdf.freq_plant_global_R1$contaminant, ps_filtered_plant_global_R1)
#ps.noncontam_plant_global_R1
#otu_df_R1_plant <- as.data.frame(otu_table(ps.noncontam_plant_global_R1))
#write.csv(otu_df_R1_plant, "otu_table_plant_global_R1_freq")
#tax_df_R1_plant <- as.data.frame(tax_table(ps.noncontam_plant_global_R1))
#write.csv(tax_df_R1_plant, "taxonomy_table_plant_global_R1_freq")
ps_filtered_plant_global_R2
ps.noncontam_plant_global_R2 <- prune_taxa(!contamdf.freq_plant_global_R2$contaminant, ps_filtered_plant_global_R2)
ps.noncontam_plant_global_R2
otu_df_R2_plant <- as.data.frame(otu_table(ps.noncontam_plant_global_R2))
write.csv(otu_df_R2_plant, "otu_table_plant_global_R2_freq.csv")
tax_df_R2_plant <- as.data.frame(tax_table(ps.noncontam_plant_global_R2))
write.csv(tax_df_R2_plant, "taxonomy_table_plant_global_R2_freq.csv")

#Joindre table OTU + TAXO   
merge_otu_taxo_plantglobalR2freq <- cbind(otu_df_R2_plant, tax_df_R2_plant)
write.csv(merge_otu_taxo_plantglobalR2freq, "merge_otu_taxo_plantglobalR2freq.csv")

#Identification des contaminants - prévalance 
#sample_data(ps_filtered_plant_global_R1)$is.neg <- sample_data(ps_filtered_plant_global_R1)$sample_or_control == "Control"
sample_data(ps_filtered_plant_global_R2)$is.neg <- sample_data(ps_filtered_plant_global_R2)$sample_or_control == "Control"
#ps.contamdf.prev_plant_global_R1 <- isContaminant(ps_filtered_plant_global_R1, method="prevalence", neg="is.neg")
contamdf.prev_plant_global_R2 <- isContaminant(ps_filtered_plant_global_R2, method="prevalence", neg="is.neg")
#table(ps.contamdf.prev_plant_global_R1$contaminant)
table(contamdf.prev_plant_global_R2$contaminant)
#head(which(ps.contamdf.prev_plant_global_R1$contaminant))
head(which(contamdf.prev_plant_global_R2$contaminant))
#contamdf.prev05_plant_global_R1 <- isContaminant(ps_filtered_plant_global_R1, method="prevalence", neg="is.neg", threshold=0.5)
contamdf.prev05_plant_global_R2 <- isContaminant(ps_filtered_plant_global_R2, method="prevalence", neg="is.neg", threshold=0.5)
#table(contamdf.prev05_plant_global_R1$contaminant)
table(contamdf.prev05_plant_global_R2$contaminant)

#export tables filtrés par prevalance 
#ps.noncontam_plant_glob_R1_prev <- prune_taxa(!contamdf.prev05_plant_global_R1$contaminant, ps_filtered_plant_global_R1)
ps.noncontam_plant_glob_R2_prev <- prune_taxa(!contamdf.prev05_plant_global_R2$contaminant, ps_filtered_plant_global_R2)

#otu_df_plant_glob_R1_prev <- as.data.frame(otu_table(ps.noncontam_plant_glob_R1_prev))
#write.csv(otu_df_plant_glob_R1_prev, "otu_table_plant_global_R1_prevalence.csv")
#tax_df_plant_glob_R1_prev <- as.data.frame(tax_table(ps.noncontam_plant_glob_R1_prev))
#write.csv(tax_df_plant_glob_R1_prev, "taxonomy_table_plant_global_R1_prevalence.csv")

otu_df_plant_glob_R2_prev <- as.data.frame(otu_table(ps.noncontam_plant_glob_R2_prev))
write.csv(otu_df_plant_glob_R2_prev, "otu_table_plant_global_R2_prevalence.csv")
tax_df_plant_glob_R2_prev <- as.data.frame(tax_table(ps.noncontam_plant_glob_R2_prev))
write.csv(tax_df_plant_glob_R2_prev, "taxonomy_table_plant_global_R2_prevalence.csv")

#Joindre table OTU + TAXO   
merge_otu_taxo_plantglobalR2prev <- cbind(otu_df_plant_glob_R2_prev, tax_df_plant_glob_R2_prev)
write.csv(merge_otu_taxo_plantglobalR2prev, "merge_otu_taxo_plantglobalR2prev.csv")

# Make phyloseq object of presence-absence in negative controls and true samples
#ps.pa_plant_global_R1 <- transform_sample_counts(ps_filtered_plant_global_R1, function(abund) 1*(abund>0))
ps.pa_plant_global_R2 <- transform_sample_counts(ps_filtered_plant_global_R2, function(abund) 1*(abund>0))
#ps.pa.neg_plant_global_R1 <- prune_samples(sample_data(ps.pa_plant_global_R1)$sample_or_control == "Control", ps.pa_plant_global_R1)
ps.pa.neg_plant_global_R2 <- prune_samples(sample_data(ps.pa_plant_global_R2)$sample_or_control == "Control", ps.pa_plant_global_R2)
#ps.pa.pos_plant_global_R1 <- prune_samples(sample_data(ps.pa_plant_global_R1)$sample_or_control == "Sample", ps.pa_plant_global_R1)
ps.pa.pos_plant_global_R2 <- prune_samples(sample_data(ps.pa_plant_global_R2)$sample_or_control == "Sample", ps.pa_plant_global_R2)

# Make data.frame of prevalence in positive and negative samples
#df.pa_plant_global_R1 <- data.frame(pa.pos_R1=taxa_sums(ps.pa.pos_plant_global_R1), pa.neg_R1=taxa_sums(ps.pa.neg_plant_global_R1),
#                       contaminant=contamdf.prev_plant_global_R1$contaminant)
#ggplot(data=df.pa_plant_global_R1, aes(x=pa.neg_R1, y=pa.pos_R1, color=contaminant)) + 
#  geom_point() +
#  xlab("Prevalence (Negative Controls)") + ylab("Prevalence (True Samples)")

df.pa_plant_global_R2 <- data.frame(pa.pos_R2=taxa_sums(ps.pa.pos_plant_global_R2), pa.neg_R2=taxa_sums(ps.pa.neg_plant_global_R2),
                       contaminant=contamdf.prev_plant_global_R2$contaminant)
ggplot(data=df.pa_plant_global_R2, aes(x=pa.neg_R2, y=pa.pos_R2, color=contaminant)) + 
  geom_point() +
  xlab("Prevalence (Negative Controls)") + ylab("Prevalence (True Samples)")


#DECONTAM TAXO LOCAL 
#Taille des librairies 
#head(sample_data(ps_plant_local_R1))
head(sample_data(ps_plant_local_R2))
#df_plant_local_R1 <- as.data.frame(sample_data(ps_plant_local_R1)) # Put sample_data into a ggplot-friendly data.frame
df_plant_local_R2 <- as.data.frame(sample_data(ps_plant_local_R2)) # Put sample_data into a ggplot-friendly data.frame
#df_plant_local_R1$LibrarySize <- sample_sums(ps_plant_local_R1)
df_plant_local_R2$LibrarySize <- sample_sums(ps_plant_local_R2)
#df_plant_local_R1 <- df_plant_local_R1[order(df_plant_local_R1$LibrarySize),]
df_plant_local_R2 <- df_plant_local_R2[order(df_plant_local_R2$LibrarySize),]
#df_plant_local_R1$Index <- seq(nrow(df_plant_local_R1))
df_plant_local_R2$Index <- seq(nrow(df_plant_local_R2))

#ggplot(data=df_plant_local_R1, aes(x=Index, y=LibrarySize, color=sample_or_control)) + 
#  geom_point()
ggplot(data=df_plant_local_R2, aes(x=Index, y=LibrarySize, color=sample_or_control)) + 
  geom_point()

#Identification des contaminants - fréquences
#sample_data(ps_plant_local_R1)$Qubit_ng_mL <- as.numeric(as.character(sample_data(ps_plant_local_R1)$Qubit_ng_mL))
sample_data(ps_plant_local_R2)$Qubit_ng_mL <- as.numeric(as.character(sample_data(ps_plant_local_R2)$Qubit_ng_mL))
#ps_plant_local_R1 <- prune_samples(!is.na(sample_data(ps_plant_local_R1)$Qubit_ng_mL) & sample_data(ps_plant_local_R1)$Qubit_ng_mL > 0, ps_plant_local_R1)
ps_plant_local_R2 <- prune_samples(!is.na(sample_data(ps_plant_local_R2)$Qubit_ng_mL) & sample_data(ps_plant_local_R2)$Qubit_ng_mL > 0, ps_plant_local_R2)
#contamdf.freq_plant_local_R1 <- isContaminant(ps_plant_local_R1, method="frequency", conc="Qubit_ng_mL")
contamdf.freq_plant_local_R2 <- isContaminant(ps_plant_local_R2, method="frequency", conc="Qubit_ng_mL")
#head(contamdf.freq_plant_local_R1)
head(contamdf.freq_plant_local_R2)
#table(contamdf.freq_plant_local_R1$contaminant)
table(contamdf.freq_plant_local_R2$contaminant)
#head(which(contamdf.freq_plant_local_R1$contaminant))
head(which(contamdf.freq_plant_local_R2$contaminant))

#ps_filtered_plant_local_R1 <- prune_samples(sample_sums(ps_plant_local_R1) > 0, ps_plant_local_R1)
ps_filtered_plant_local_R2 <- prune_samples(sample_sums(ps_plant_local_R2) > 0, ps_plant_local_R2)
#length(sample_names(ps_filtered_plant_local_R1))
length(sample_names(ps_filtered_plant_local_R2))

#plot_frequency(ps_filtered_plant_local_R1, taxa_names(ps_filtered_plant_local_R1), conc="Qubit_ng_mL") + 
#  xlab("DNA Concentration (Qubit)")
plot_frequency(ps_filtered_plant_local_R2, taxa_names(ps_filtered_plant_local_R2), conc="Qubit_ng_mL") + 
  xlab("DNA Concentration (Qubit)")

#table(contamdf.freq_plant_local_R1$contaminant) #que des FALSE = pas de contaminants détéctés 
#set.seed(100)
#plot_frequency(ps_filtered_plant_local_R1, taxa_names(ps_filtered_plant_local_R1)[sample(which(contamdf.freq_plant_local_R1$contaminant),3)], conc="Qubit_ng_mL") +
#  xlab("DNA Concentration (Qubit)")
table(contamdf.freq_plant_local_R2$contaminant) #que des FALSE = pas de contaminants détéctés 
set.seed(100)
plot_frequency(ps_filtered_plant_local_R2, taxa_names(ps_filtered_plant_local_R2)[sample(which(contamdf.freq_plant_local_R2$contaminant),3)], conc="Qubit_ng_mL") +
  xlab("DNA Concentration (Qubit)")

#export tables filtrés par frequence 
#ps_filtered_plant_local_R1
#ps.noncontam_plant_local_R1 <- prune_taxa(!contamdf.freq_plant_local_R1$contaminant, ps_filtered_plant_local_R1)
#ps.noncontam_plant_local_R1
#otu_df_R1_plant_loc <- as.data.frame(otu_table(ps.noncontam_plant_local_R1))
#write.csv(otu_df_R1_plant_loc, "otu_table_plant_local_R1_freq")
#tax_df_R1_plant_loc <- as.data.frame(tax_table(ps.noncontam_plant_local_R1))
#write.csv(tax_df_R1_plant_loc, "taxonomy_table_plant_local_R1_freq")
ps_filtered_plant_local_R2
ps.noncontam_plant_local_R2 <- prune_taxa(!contamdf.freq_plant_local_R2$contaminant, ps_filtered_plant_local_R2)
ps.noncontam_plant_local_R2
otu_df_R2_plant_loc <- as.data.frame(otu_table(ps.noncontam_plant_local_R2))
write.csv(otu_df_R2_plant_loc, "otu_table_plant_local_R2_freq.csv")
tax_df_R2_plant_loc <- as.data.frame(tax_table(ps.noncontam_plant_local_R2))
write.csv(tax_df_R2_plant_loc, "taxonomy_table_plant_local_R2_freq.csv")

#Joindre table OTU + TAXO   
merge_otu_taxo_plantlocaalR2freq <- cbind(otu_df_R2_plant_loc, tax_df_R2_plant_loc)
write.csv(merge_otu_taxo_plantlocaalR2freq, "merge_otu_taxo_plantlocalR2freq.csv")

#ps_filtered_plant_local_R1
#ps.noncontam_plant_gloabal_R1 <- prune_taxa(!contamdf.freq_plant_local_R1$contaminant, ps_filtered_plant_local_R1)
#ps.noncontam_plant_gloabal_R1
ps_filtered_plant_local_R2
ps.noncontam_plant_gloabal_R2 <- prune_taxa(!contamdf.freq_plant_local_R2$contaminant, ps_filtered_plant_local_R2)
ps.noncontam_plant_gloabal_R2

#Identification des contaminants - prévalance 
#sample_data(ps_filtered_plant_local_R1)$is.neg <- sample_data(ps_filtered_plant_local_R1)$sample_or_control == "Control"
sample_data(ps_filtered_plant_local_R2)$is.neg <- sample_data(ps_filtered_plant_local_R2)$sample_or_control == "Control"
#ps.contamdf.prev_plant_local_R1 <- isContaminant(ps_filtered_plant_local_R1, method="prevalence", neg="is.neg")
contamdf.prev_plant_local_R2 <- isContaminant(ps_filtered_plant_local_R2, method="prevalence", neg="is.neg")
#table(ps.contamdf.prev_plant_local_R1$contaminant)
table(contamdf.prev_plant_local_R2$contaminant)
#head(which(ps.contamdf.prev_plant_local_R1$contaminant))
head(which(contamdf.prev_plant_local_R2$contaminant))
#contamdf.prev05_plant_local_R1 <- isContaminant(ps_filtered_plant_local_R1, method="prevalence", neg="is.neg", threshold=0.5)
contamdf.prev05_plant_local_R2 <- isContaminant(ps_filtered_plant_local_R2, method="prevalence", neg="is.neg", threshold=0.5)
#table(contamdf.prev05_plant_local_R1$contaminant)
table(contamdf.prev05_plant_local_R2$contaminant)

#export tables filtrés par prevalance 
#ps.noncontam_plant_loc_R1_prev <- prune_taxa(!contamdf.prev05_plant_local_R1$contaminant, ps_filtered_plant_local_R1)
ps.noncontam_plant_loc_R2_prev <- prune_taxa(!contamdf.prev05_plant_local_R2$contaminant, ps_filtered_plant_local_R2)

#otu_df_plant_loc_R1_prev <- as.data.frame(otu_table(ps.noncontam_plant_loc_R1_prev))
#write.csv(otu_df_plant_loc_R1_prev, "otu_table_plant_local_R1_prevalence.csv")
#tax_df_plant_loc_R1_prev <- as.data.frame(tax_table(ps.noncontam_plant_loc_R1_prev))
#write.csv(tax_df_plant_loc_R1_prev, "taxonomy_table_plant_local_R1_prevalence.csv")

otu_df_plant_loc_R2_prev <- as.data.frame(otu_table(ps.noncontam_plant_loc_R2_prev))
write.csv(otu_df_plant_loc_R2_prev, "otu_table_plant_local_R2_prevalence.csv")
tax_df_plant_loc_R2_prev <- as.data.frame(tax_table(ps.noncontam_plant_loc_R2_prev))
write.csv(tax_df_plant_loc_R2_prev, "taxonomy_table_plant_local_R2_prevalence.csv")

#Joindre table OTU + TAXO   
merge_otu_taxo_plantlocaalR2prev <- cbind(otu_df_plant_loc_R2_prev, tax_df_plant_loc_R2_prev)
write.csv(merge_otu_taxo_plantlocaalR2prev, "merge_otu_taxo_plantlocalR2prev.csv")

# Make phyloseq object of presence-absence in negative controls and true samples
#ps.pa_plant_local_R1 <- transform_sample_counts(ps_filtered_plant_local_R1, function(abund) 1*(abund>0))
ps.pa_plant_local_R2 <- transform_sample_counts(ps_filtered_plant_local_R2, function(abund) 1*(abund>0))
#ps.pa.neg_plant_local_R1 <- prune_samples(sample_data(ps.pa_plant_local_R1)$sample_or_control == "Control", ps.pa_plant_local_R1)
ps.pa.neg_plant_local_R2 <- prune_samples(sample_data(ps.pa_plant_local_R2)$sample_or_control == "Control", ps.pa_plant_local_R2)
#ps.pa.pos_plant_local_R1 <- prune_samples(sample_data(ps.pa_plant_local_R1)$sample_or_control == "Sample", ps.pa_plant_local_R1)
ps.pa.pos_plant_local_R2 <- prune_samples(sample_data(ps.pa_plant_local_R2)$sample_or_control == "Sample", ps.pa_plant_local_R2)

# Make data.frame of prevalence in positive and negative samples
#df.pa_plant_local_R1 <- data.frame(pa.pos_R1=taxa_sums(ps.pa.pos_plant_local_R1), pa.neg_R1=taxa_sums(ps.pa.neg_plant_local_R1),
#                       contaminant=contamdf.prev_plant_local_R1$contaminant)
#ggplot(data=df.pa_plant_local_R1, aes(x=pa.neg_R1, y=pa.pos_R1, color=contaminant)) + 
#  geom_point() +
#  xlab("Prevalence (Negative Controls)") + ylab("Prevalence (True Samples)")

df.pa_plant_local_R2 <- data.frame(pa.pos_R2=taxa_sums(ps.pa.pos_plant_local_R2), pa.neg_R2=taxa_sums(ps.pa.neg_plant_local_R2),
                                    contaminant=contamdf.prev_plant_local_R2$contaminant)
ggplot(data=df.pa_plant_local_R2, aes(x=pa.neg_R2, y=pa.pos_R2, color=contaminant)) + 
  geom_point() +
  xlab("Prevalence (Negative Controls)") + ylab("Prevalence (True Samples)")

