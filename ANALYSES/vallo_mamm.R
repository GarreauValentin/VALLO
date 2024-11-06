#Importation des libraries 
library(phyloseq)
library(ggplot2)
library(decontam)
library(readxl)
library(dplyr)
library(reshape2)
library(tidyr)
library(vegan)

#Import des données 
#Data externe 
extdata <- read.table(file = "extdata_mamm.tsv", sep = "\t", header=T)
#Table OTUs
data_mamm_R2_freq <- read.table(file = "merge_otu_taxo_mammR2freq.csv", header = TRUE, sep = ",", row.names = "X")

#formattage donnéés
data_mamm_R2_freq <- subset(data_mamm_R2_freq, Phylum != "")
data_mamm_div <- data_mamm_R2_freq
data_mamm_div <- data_mamm_div[, !names(data_mamm_div) %in% c("ETEX_BHUM", "ETF_BHUM","STEX1_BHUM")]
data_mamm_div$sum_1 <- rowSums(data_mamm_div[, grep("1", names(data_mamm_div))])
data_mamm_div$sum_2 <- rowSums(data_mamm_div[, grep("2", names(data_mamm_div))])
data_mamm_div$sum_3 <- rowSums(data_mamm_div[, grep("3", names(data_mamm_div))])
data_mamm_div$sum_4 <- rowSums(data_mamm_div[, grep("4", names(data_mamm_div))])
data_mamm_div <- data_mamm_div[, 33:ncol(data_mamm_div)]
# Calculer la somme des colonnes sum_1, sum_2, sum_3, sum_4 pour chaque ligne
#sum_values <- rowSums(data_mamm_div[, c("sum_1", "sum_2", "sum_3", "sum_4")], na.rm = TRUE)
#new_row <- colSums(data_mamm_div[, c("sum_1", "sum_2", "sum_3", "sum_4")], na.rm = TRUE)
#data_mamm_div <- rbind(data_mamm_div, c(NA, NA, NA, NA, NA, NA, new_row))
#rownames(data_mamm_div)[nrow(data_mamm_div)] <- "sum_asv"

#AUTRES
colnames(data_mamm_R2_freq)
echeau <- c(1:12,28:31)
echsed <- c(15:26,32:35)
tabeau <- colSums(data_mamm_R2_freq[,echeau])
tabsed <- colSums(data_mamm_R2_freq[,echsed])
summary(tabeau)
summary(tabsed)
boxplot(tabeau)
boxplot(tabsed)
shapiro.test(tabeau)
shapiro.test(tabsed)
var.test(tabeau, tabsed)
t.test(tabeau, tabsed)
wilcox.test(tabeau, tabsed)

#RAREFACTION 
rarefied_table_mamm <- rrarefy(data_mamm_div[,7:10], 94939)

#MAMM_R2_freq
# Initialiser un tableau vide pour stocker le résultat
result_mamm_R2_freq <- data.frame()

# Boucle sur chaque colonne de données pour extraire les informations taxonomiques
for (col in colnames(data_mamm_R2_freq)[!colnames(data_mamm_R2_freq) %in% c("Phylum", "Class", "Order", "Family", "Genus", "Species")]) {
  
  # Filtrer les lignes où la valeur est > 0 dans la colonne actuelle
  rows_with_data <- data_mamm_R2_freq[data_mamm_R2_freq[[col]] > 0, ]
  
  # Si des lignes existent, ajouter leurs informations taxonomiques au tableau résultat
  if (nrow(rows_with_data) > 0) {
    tax_info <- rows_with_data[, c("Phylum", "Class", "Order", "Family", "Genus", "Species")]
    
    # Générer des noms de ligne uniques en combinant le nom de la colonne et l’indice de la ligne
    rownames(tax_info) <- paste(col, seq_len(nrow(tax_info)), sep = "_")
    
    # Ajouter les informations au tableau résultat
    result_mamm_R2_freq <- rbind(result_mamm_R2_freq, tax_info)
  }
}

#noms lignes 
result_mamm_R2_freq$Sample <- sub("_.*", "", rownames(result_mamm_R2_freq))
#merge données extérieurs et resultats
result_mamm_R2_freq <- merge(result_mamm_R2_freq, extdata, by.x = "Sample", by.y = "Echantillons", all.x = TRUE)

#ANALYSES DONNEES
ggplot(data = result_mamm_R2_freq, aes(x = Altitude_.m., y = Species)) +
  geom_point() + # You can add a layer, for example, points
  labs(x = "Altitude (m)", y = "Species", title = "Species vs Altitude") +
  theme_minimal()
