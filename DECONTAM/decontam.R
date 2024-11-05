#Importation des libraries 
library(phyloseq)
library(ggplot2)
library(decontam)
library(readxl)
library(dplyr)
library(reshape2)
library(tidyr)


qubit_data <- read_excel("~/Documents/MASTER/M2/ATELIERS/R/DECONTAM/Qubit.xlsx")
#Tbale OTUs
#Attention si # devant le titre des colonnes pas pris en compte 
asv_data_mamm_R1 <- read.csv("~/Documents/MASTER/M2/ATELIERS/R/DECONTAM/mammal_run1_curated_ASVtable.csv", header = TRUE, sep = ",")
asv_data_mamm_R2 <- read.table("~/Documents/MASTER/M2/ATELIERS/R/DECONTAM/mammal_run2_curated_ASVtable.csv", header = TRUE, sep = ",")
str(asv_data_mamm_R1)

# Créer un tableau de caractéristiques (otu_table)
# Assurez-vous que la première colonne (X) est l'index des ASVs et que les autres colonnes sont les échantillons
otu_table <- asv_data_mamm_R1 %>%
  column_to_rownames(var = "X") %>%  # Définir la colonne X comme noms de lignes
  as.matrix() %>%  # Convertir en matrice
  otu_table(taxa_are_rows = TRUE)  # Créer l'OTU table, les taxa sont en lignes

# Créer un tableau de métadonnées (sample_data)
# Il faut que la colonne des échantillons dans le tableau de métadonnées soit alignée avec l'OTU table
# Assurez-vous que les échantillons de qubit_data sont en accord avec ceux d'asv_data_mamm_R1
qubit_data_sample_data <- qubit_data %>%
  select(Echantillons, Qubit_ng_mL, Concentration_dans_le_tube_(ng/µL)) %>%
  column_to_rownames(var = "Echantillons") %>%
  sample_data()  # Créer le sample data

# Créer l'objet phyloseq
ps <- phyloseq(otu_table, qubit_data_sample_data)

# Vérifiez l'objet phyloseq
print(ps)

# Visualiser les données de l'objet phyloseq
head(sample_data(ps))  # Afficher les premières lignes des métadonnées
