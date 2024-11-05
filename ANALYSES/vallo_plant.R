#Importation des libraries 
library(phyloseq)
library(ggplot2)
library(decontam)
library(readxl)
library(dplyr)
library(reshape2)
library(tidyr)

#Import des données GLOBAL
#Data externe 
extdata <- read.table(file = "extdata_plant.tsv", sep = "\t", header=T)
#Table OTUs
data_plant_R2_freq <- read.table(file = "merge_otu_taxo_plantglobalR2freq.csv", header = TRUE, sep = ",", row.names = "X")

#formarrage donnéés
data_plant_R2_freq <- subset(data_plant_R2_freq, Subgroups != "")


#plant_R2_freq
# Initialiser un tableau vide pour stocker le résultat
result_plant_R2_freq <- data.frame()

# Boucle sur chaque colonne de données pour extraire les informations taxonomiques
for (col in colnames(data_plant_R2_freq)[!colnames(data_plant_R2_freq) %in% c("Subgroups", "Infragroups", "class", "Order", "Family", "Genus", "Species")]) {
  
  # Filtrer les lignes où la valeur est > 0 dans la colonne actuelle
  rows_with_data <- data_plant_R2_freq[data_plant_R2_freq[[col]] > 0, ]
  
  # Si des lignes existent, ajouter leurs informations taxonomiques au tableau résultat
  if (nrow(rows_with_data) > 0) {
    tax_info <- rows_with_data[, c("Subgroups", "Infragroups", "Class", "Order", "Family", "Genus", "Species")]
    
    # Générer des noms de ligne uniques en combinant le nom de la colonne et l’indice de la ligne
    rownames(tax_info) <- paste(col, seq_len(nrow(tax_info)), sep = "_")
    
    # Ajouter les informations au tableau résultat
    result_plant_R2_freq <- rbind(result_plant_R2_freq, tax_info)
  }
}

#noms lignes 
result_plant_R2_freq$Sample <- sub("_.*", "", rownames(result_plant_R2_freq))
#merge données extérieurs et resultats
result_plant_R2_freq <- merge(result_plant_R2_freq, extdata, by.x = "Sample", by.y = "Echantillons", all.x = TRUE)

#ANALYSES DONNEES
ggplot(data = result_plant_R2_freq, aes(x = Altitude_.m., y = Species)) +
  geom_point() + # You can add a layer, for example, points
  labs(x = "Altitude (m)", y = "Species", title = "Species vs Altitude") +
  theme_minimal()

