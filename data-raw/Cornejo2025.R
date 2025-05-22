## code to prepare `Cornejo2025` dataset goes here

library(tidyverse)
library(parafac4microbiome)

# Tongue
tongue = read.csv("./data-raw/Cornejo2025/20240503_UNOISE_new/tongueCounts.csv", header=FALSE) %>% as_tibble()
tongue_taxonomy = read.csv("./data-raw/Cornejo2025/20240503_UNOISE_new/taxonomyTongue_fixed.csv", sep=" ") %>% as_tibble()
tongue_sampleMeta = read.csv("./data-raw/Cornejo2025/20240503_UNOISE_new/tongueSampleMeta.csv",header=FALSE) %>% as_tibble()

colnames(tongue) = tongue_taxonomy$zOTU
temp = read.csv("./data-raw/Cornejo2025/sampleInfo_fixed.csv", sep=" ") %>% as_tibble()
colnames(tongue_sampleMeta) = c(temp %>% select(-Description,-subject,-newTimepoint) %>% colnames, "Description", "subject", "newTimepoint")

# Load other metadata and ph
ph_BOMP = read_delim("./data-raw/Cornejo2025/GOH-TRANS_csv_export_20240205114955/GOH-TRANS_export_20240205.csv",
                     delim = ";", escape_double = FALSE, trim_ws = TRUE) %>% as_tibble()

df1 = ph_BOMP %>% select(`Participant Id`, starts_with("5.")) %>% mutate(subject = 1:42, numTeeth = `5.1|Number of teeth`, DMFT = `5.2|DMFT`, numBleedingSites = `5.3|Bleeding sites`, boppercent = `5.4|BOP%`, DPSI = `5.5|DPSI`, pH = `5.8|pH`) %>% select(subject, numTeeth, DMFT, numBleedingSites, boppercent, DPSI, pH)
df2 = ph_BOMP %>% select(`Participant Id`, starts_with("12.")) %>% mutate(subject = 1:42, numTeeth = `12.1|Number of teeth`, DMFT = `12.2|DMFT`, numBleedingSites = `12.3|Bleeding sites`, boppercent = `12.4|BOP%`, DPSI = `12.5|DPSI`, pH = `12.8|pH`) %>% select(subject, numTeeth, DMFT, numBleedingSites, boppercent, DPSI, pH)
df3 = ph_BOMP %>% select(`Participant Id`, starts_with("19.")) %>% mutate(subject = 1:42, numTeeth = `19.1|Number of teeth`, DMFT = `19.2|DMFT`, numBleedingSites = `19.3|Bleeding sites`, boppercent = `19.4|BOP%`, DPSI = `19.5|DPSI`, pH = `19.8|pH`) %>% select(subject, numTeeth, DMFT, numBleedingSites, boppercent, DPSI, pH)
df4 = ph_BOMP %>% select(`Participant Id`, starts_with("26.")) %>% mutate(subject = 1:42, numTeeth = `26.1|Number of teeth`, DMFT = `26.2|DMFT`, numBleedingSites = `26.3|Bleeding sites`, boppercent = `26.4|BOP%`, DPSI = `26.5|DPSI`, pH = `26.8|pH`) %>% select(subject, numTeeth, DMFT, numBleedingSites, boppercent, DPSI, pH)

otherMeta = rbind(df1, df2, df3, df4) %>% as_tibble() %>% mutate(newTimepoint = rep(c(0,3,6,12), each=42))

# Remove pH measurements of lower than 0
otherMeta = otherMeta[!otherMeta$pH < 1,] %>% as_tibble()

# New approach per 20250402:

# Remove timepoint 9
mask = tongue_sampleMeta$newTimepoint != 9
tongue = tongue[mask,]
tongue_sampleMeta = tongue_sampleMeta[mask,]

# Put into cube
tongue = parafac4microbiome::reshapeData(tongue, tongue_sampleMeta$subject, tongue_taxonomy, tongue_sampleMeta$newTimepoint)

# Repair modes to avoid breaking changes
colnames(tongue$mode1) = c("subject", "index")
tongue$mode1 = tongue$mode1 %>% left_join(tongue_sampleMeta %>% mutate(subject=as.character(subject)) %>% select(subject,GenderID) %>% unique()) %>% select(-index)
tongue$mode2 = tongue$mode2 %>% select(-index)
tongue$mode3 = tongue$mode3 %>% mutate(newTimepoint=timepointMetadata) %>% select(-timepointMetadata,-index)

# Old approach
#
# Put into cube
# I = tongue_sampleMeta$subject %>% unique() %>% length()
# J = ncol(tongue)
# K = 4
# timepoints = c(0, 3, 6, 12)
# tongueCube = array(0L, dim=c(I,J,K))
#
# for(k in 1:K){
#   temp = cbind(tongue, tongue_sampleMeta) %>% as_tibble()
#   tongueCube[,,k] = temp %>%
#     filter(newTimepoint == timepoints[k]) %>%
#     right_join(tongue_sampleMeta %>% select(subject) %>% unique()) %>%
#     arrange(subject) %>%
#     select(-all_of(colnames(tongue_sampleMeta))) %>%
#     as.matrix()
# }
#
# tongueCube_mode1 = tongue_sampleMeta %>% select(subject, GenderID) %>% arrange(subject) %>% unique()
# tongueCube_mode2 = tongue_taxonomy
# tongueCube_mode3 = tongue_sampleMeta %>% filter(newTimepoint %in% timepoints) %>% select(newTimepoint) %>% unique()
#
# tongueData = list("data"=tongueCube, "mode1"=tongueCube_mode1, "mode2"=tongueCube_mode2, "mode3"=tongueCube_mode3)

# Saliva
saliva = read.csv("./data-raw/Cornejo2025/20240503_UNOISE_new/salivaCounts.csv", header=FALSE) %>% as_tibble()
saliva_taxonomy = read.csv("./data-raw/Cornejo2025/20240503_UNOISE_new/taxonomysaliva_fixed.csv", sep=" ") %>% as_tibble()
saliva_sampleMeta = read.csv("./data-raw/Cornejo2025/20240503_UNOISE_new/salivaSampleMeta.csv",header=FALSE) %>% as_tibble()

colnames(saliva) = saliva_taxonomy$zOTU
temp = read.csv("./data-raw/Cornejo2025/sampleInfo_fixed.csv", sep=" ") %>% as_tibble()
colnames(saliva_sampleMeta) = c(temp %>% select(-Description,-subject,-newTimepoint) %>% colnames, "Description", "subject", "newTimepoint")

# New approach per 20250402:

# Remove timepoint 9
mask = saliva_sampleMeta$newTimepoint != 9
saliva = saliva[mask,]
saliva_sampleMeta = saliva_sampleMeta[mask,]

# Put into cube
saliva = parafac4microbiome::reshapeData(saliva, saliva_sampleMeta$subject, saliva_taxonomy, saliva_sampleMeta$newTimepoint)

# Repair modes to avoid breaking changes
colnames(saliva$mode1) = c("subject", "index")
saliva$mode1 = saliva$mode1 %>% left_join(saliva_sampleMeta %>% mutate(subject=as.character(subject)) %>% select(subject,GenderID) %>% unique()) %>% select(-index)
saliva$mode2 = saliva$mode2 %>% select(-index)
saliva$mode3 = saliva$mode3 %>% mutate(newTimepoint=timepointMetadata) %>% select(-timepointMetadata,-index)

# Old approach:
#
# Put into cube
# I = saliva_sampleMeta$subject %>% unique() %>% length()
# J = ncol(saliva)
# K = 4
# timepoints = c(0, 3, 6, 12)
# salivaCube = array(0L, dim=c(I,J,K))
#
# for(k in 1:K){
#   temp = cbind(saliva, saliva_sampleMeta) %>% as_tibble()
#   salivaCube[,,k] = temp %>%
#     filter(newTimepoint == timepoints[k]) %>%
#     right_join(saliva_sampleMeta %>% select(subject) %>% unique()) %>%
#     arrange(subject) %>%
#     select(-all_of(colnames(saliva_sampleMeta))) %>%
#     as.matrix()
# }
#
# salivaCube_mode1 = saliva_sampleMeta %>% select(subject, GenderID) %>% arrange(subject) %>% unique()
# salivaCube_mode2 = saliva_taxonomy
# salivaCube_mode3 = saliva_sampleMeta %>% filter(newTimepoint %in% timepoints) %>% select(newTimepoint) %>% unique()
#
# salivaData = list("data"=salivaCube, "mode1"=salivaCube_mode1, "mode2"=salivaCube_mode2, "mode3"=salivaCube_mode3)

# Process microbiome data
# processedTongue = processDataCube(tongueData, sparsityThreshold = 0.5, considerGroups=TRUE, groupVariable="GenderID", CLR=TRUE, centerMode=1, scaleMode=2)
# processedSaliva = processDataCube(salivaData, sparsityThreshold = 0.5, considerGroups=TRUE, groupVariable="GenderID", CLR=TRUE, centerMode=1, scaleMode=2)

# Cytokines
df = read.csv("./data-raw/Cornejo2025/20241209_cytokines.csv", header=FALSE, sep=" ") %>% as_tibble()
featureMeta = read.csv("./data-raw/Cornejo2025/20241209_cytokines_featureMeta.csv", header=FALSE) %>% as_tibble()
sampleInfo = read.csv("./data-raw/Cornejo2025/20241209_cytokines_sampleMeta.csv", header=FALSE, sep=" ") %>% as_tibble()
colnames(sampleInfo) = c("subject", "GenderID", "newTimepoint", "unknown", "unknown2")

# Remove dilution factor using total protein
dilutionFactor = sampleInfo %>% left_join(saliva_sampleMeta) %>% select(Total_protein_ug_ml) %>% pull()
df = sweep(df, 1, dilutionFactor, FUN="/")

# New approach per 20250402
cytokine = parafac4microbiome::reshapeData(df, sampleInfo$subject, featureMeta, sampleInfo$newTimepoint)

# Repair modes to avoid breaking changes
colnames(cytokine$mode1) = c("subject", "index")
cytokine$mode1 = cytokine$mode1 %>% left_join(saliva_sampleMeta %>% mutate(subject=as.character(subject)) %>% select(subject,GenderID) %>% unique()) %>% select(-index)
cytokine$mode2 = cytokine$mode2 %>% select(-index)
cytokine$mode3 = cytokine$mode3 %>% mutate(newTimepoint=timepointMetadata) %>% select(-timepointMetadata,-index)

# Remove subject 9 and 19 due to being outliers
cytokine$data = cytokine$data[-c(9,19),,]
cytokine$mode1 = cytokine$mode1[-c(9,19),]

# Process
cytokine$data = log(cytokine$data + 0.1300)
cytokine$data = multiwayCenter(cytokine$data, 1)
cytokine$data = multiwayScale(cytokine$data, 2)

# Old approach
#
# Put into cube
# I = 27
# J = 22
# K = 3
# timepoints = c(0, 3, 12)
# cytokineCube = array(0L, dim=c(I,J,K))
#
# for(k in 1:K){
#   temp = cbind(df, sampleInfo) %>% as_tibble()
#   cytokineCube[,,k] = temp %>%
#     filter(newTimepoint == timepoints[k]) %>%
#     right_join(sampleInfo %>% select(subject) %>% unique()) %>%
#     arrange(subject) %>%
#     select(-all_of(colnames(sampleInfo))) %>%
#     as.matrix()
# }
#
# # Log transform
# pseudocount = 0.1300
# cytokineCube_log = log(cytokineCube + pseudocount)
#
# # Center and scale
# cytokineCube_cnt = multiwayCenter(cytokineCube_log, mode=1)
# cytokineCube_cnt_scl = multiwayScale(cytokineCube_cnt, mode=2)
#
# # Prep metadata
# cytokineCube_mode1 = sampleInfo %>% select(subject, GenderID) %>% arrange(subject) %>% unique()
# cytokineCube_mode2 = featureMeta
# cytokineCube_mode3 = sampleInfo %>% filter(newTimepoint %in% timepoints) %>% select(newTimepoint) %>% unique()
#
# cytokineData = list("data"=cytokineCube_cnt_scl, "mode1"=cytokineCube_mode1, "mode2"=cytokineCube_mode2, "mode3"=cytokineCube_mode3)

Cornejo2025 = list()
Cornejo2025$Tongue = tongue
Cornejo2025$Saliva = saliva
Cornejo2025$Cytokines = cytokine
usethis::use_data(Cornejo2025, overwrite = TRUE)
