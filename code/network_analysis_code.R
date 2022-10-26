###############################NETWORK ANALYSIS
#creating lists for network analysis 

#if (!requireNamespace("BiocManager", quietly = TRUE))
# install.packages("BiocManager")
#BiocManager::install("OmnipathR")

library(readxl)
library(dplyr)
library(OmnipathR)
library(tidyverse)
library(dplyr) # data manipulation
library(ggforce) # extension of ggplot. 
library(ggplot2)

setwd("C:/Users/czagami/NanostringData")

zonality <- read_excel("zonality.xlsx")
zonalitydf <- as.data.frame(zonality, col.names = T, row.names = T)
rownames(zonalitydf) <- zonalitydf[,1]
zonalitydf[,1] <- NULL


#import interactions from omnipath
all_interactions <-  OmnipathR::import_all_interactions(
  organism = 9606,
  directed = 'no'
)

#filter the interactions to only contain protein-protein and TF-TG interactions

#get_interaction_resources()
#interactions = import_omnipath_interactions() %>% as_tibble()
all_interactions %>% distinct(type) #checking all the interactions present 


interactions <- dplyr::filter(all_interactions, type == "post_translational" | type == "transcriptional")

unique(all_interactions$dorothea_level)
#Filter for dorotheal levels a, b, c and NAs.

interactions <- dplyr::filter(all_interactions, dorothea_level != "D" | dorothea_level != "A;D")

#Foveola epithelium 
fov_epi <- select(zonalitydf,'FOV EPI')
fov_epi_filtr <- filter(fov_epi, `FOV EPI` > 11.0)
fov_epi_express <- c(rownames(fov_epi_filtr))

interactions_fov_epi <- filter(interactions, source_genesymbol %in% fov_epi_express & 
                                 target_genesymbol %in% fov_epi_express) 

#Isthmus epithelium 
ist_epi <- select(zonalitydf,'IST EPI')
ist_epi_filtr <- filter(ist_epi, `IST EPI` > 11.0)
ist_epi_express <- c(rownames(ist_epi_filtr))
interactions_ist_epi <- filter(interactions, source_genesymbol %in% ist_epi_express & 
                                 target_genesymbol %in% ist_epi_express) 


#try to make a loop

zones_express = list()
for (i in colnames(zonalitydf)){
  zones_express[[i]] <- filter(zonalitydf, zonalitydf[[i]] > 11.0)
}

#neck epithelium
neck_epi_express <- (rownames(zones_express[['NECK EPI']]))
interactions_neck_epi <- filter(interactions, source_genesymbol %in% neck_epi_express & 
                                  target_genesymbol %in% neck_epi_express)

#base epithelium
base_epi_express <- (rownames(zones_express[['BASE EPI']]))
interactions_base_epi <- filter(interactions, source_genesymbol %in% base_epi_express & 
                                  target_genesymbol %in% base_epi_express)

#Fov stroma

fov_str_express <- (rownames(zones_express[['FOV STR']]))
interactions_fov_str <- filter(interactions, source_genesymbol %in% fov_str_express & 
                                 target_genesymbol %in% fov_str_express)

#isthmus stroma

ist_str_express <- (rownames(zones_express[['IST STR']]))
interactions_ist_str <- filter(interactions, source_genesymbol %in% ist_str_express & 
                                 target_genesymbol %in% ist_str_express)

#neck stroma

neck_str_express <- (rownames(zones_express[['NECK STR']]))
interactions_neck_str <- filter(interactions, source_genesymbol %in% neck_str_express & 
                                  target_genesymbol %in% neck_str_express)

#Base stroma 

base_str_express <- (rownames(zones_express[['BASE STR']]))
interactions_base_str <- filter(interactions, source_genesymbol %in% base_str_express & 
                                  target_genesymbol %in% base_str_express)

#muscularis mucosa

musc_str_express <- (rownames(zones_express[['MUSC STR']]))
interactions_musc_str <- filter(interactions, source_genesymbol %in% musc_str_express & 
                                  target_genesymbol %in% musc_str_express)

#create the networks
write.table(interactions_fov_epi, file="interactions_fov_epi.csv", sep = ",")
write.table(interactions_ist_epi, file="interactions_ist_epi.csv", sep = ",")
write.table(interactions_neck_epi, file="interactions_neck_epi.csv", sep = ",")
write.table(interactions_base_epi, file="interactions_base_epi.csv", sep = ",")
write.table(interactions_fov_str, file="interactions_fov_str.csv", sep = ",")
write.table(interactions_ist_str, file="interactions_ist_str.csv", sep = ",")
write.table(interactions_neck_str, file="interactions_neck_str.csv", sep = ",")
write.table(interactions_base_str, file="interactions_base_str.csv", sep = ",")
write.table(interactions_musc_str, file="interactions_musc_str.csv", sep = ",")

write.table(gradients_df, file = 'list_gradients.csv', sep = ",")
