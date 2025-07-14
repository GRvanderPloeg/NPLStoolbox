## code to prepare `Jakobsen2025` dataset goes here
library(tidyverse)
library(ggplot2)
library(ggpubr)
library(stringr)
library(parafac4microbiome)
library(CMTFtoolbox)

# Load subject info
infant_anthropometrics = read.csv("./data-raw/Jakobsen2025/infant_anthropometrics.csv") %>% as_tibble()

# Faeces subject info
sampleInfo = read.csv("./data-raw/Jakobsen2025/faeces_sampleMeta.csv", header=FALSE, sep=" ") %>% as_tibble()
colnames(sampleInfo) = c("Sample", "RCID", "BMI", "BMI.group", "Days", "Gestational.age", "C.section", "AB.infant", "AB.mother", "Secretor", "Lewis", "subject")
sampleInfo = sampleInfo %>% left_join(infant_anthropometrics %>% select(RCID, whz.6m))
subjectMeta_faeces = sampleInfo %>% select(subject, BMI, BMI.group, C.section, Secretor, Lewis, whz.6m) %>% unique() %>% mutate(subject=as.character(subject)) %>% arrange(subject)
subjectMeta_faeces = subjectMeta_faeces[-90,]

# Milk subject info
sampleInfo = read.csv("./data-raw/Jakobsen2025/milk_sampleMeta.csv", header=FALSE, sep=" ") %>% as_tibble()
colnames(sampleInfo) = c("Sample", "RCID", "BMI", "BMI.group", "Days", "Gestational.age", "C.section", "AB.infant", "AB.mother", "Secretor", "Lewis", "subject")
subjectMeta = sampleInfo %>% select(subject, BMI, BMI.group, C.section, Secretor, Lewis) %>% unique() %>% mutate(subject=as.character(subject)) %>% arrange(subject)
subjectMeta_milk = subjectMeta %>% left_join(subjectMeta_faeces)

# Milk metab subject info
sampleInfo = read.csv("./data-raw/Jakobsen2025/milkMetab_sampleMeta.csv", header=FALSE, sep=" ") %>% as_tibble()
colnames(sampleInfo) = c("RCID", "BMI", "BMI.group", "Days", "Gestational.age", "C.section", "AB.infant", "AB.mother", "Secretor", "Lewis", "subject")
subjectMeta = sampleInfo %>% select(subject, BMI, BMI.group, C.section, Secretor, Lewis) %>% unique() %>% mutate(subject=as.character(subject)) %>% arrange(subject)
subjectMeta_milkMetab = subjectMeta %>% left_join(subjectMeta_faeces %>% select(subject, whz.6m))

# Identify shared subjects
sharedSubjects = intersect(intersect(subjectMeta_faeces$subject, subjectMeta_milk$subject), subjectMeta_milkMetab$subject)
homogenized_subjectMeta_bmi = subjectMeta_faeces %>% filter(subject %in% sharedSubjects) %>% arrange(subject)
homogenized_subjectMeta_whz = homogenized_subjectMeta_bmi %>% filter(whz.6m != "NA")

# Faecal microbiome
df = read.csv("./data-raw/Jakobsen2025/faecesCounts.csv", header=FALSE, sep=" ") %>% as_tibble()
taxonomy = read.csv("./data-raw/Jakobsen2025/newTaxonomy_faeces.csv", header=FALSE, sep=" ") %>% as_tibble()
sampleInfo = read.csv("./data-raw/Jakobsen2025/faeces_sampleMeta.csv", header=FALSE, sep=" ") %>% as_tibble()
colnames(sampleInfo) = c("Sample", "RCID", "BMI", "BMI.group", "Days", "Gestational.age", "C.section", "AB.infant", "AB.mother", "Secretor", "Lewis", "subject")

subjectMeta = sampleInfo %>% select(subject, BMI, BMI.group, C.section, Secretor, Lewis) %>% unique() %>% mutate(subject=as.character(subject)) %>% arrange(subject)

# Filter taxa to taxa with at least 1 non-zero value
featureMask = colSums(df) > 0
df = df[,featureMask]
taxonomy = taxonomy[featureMask,]

# New approach per 20250402:
faeces = parafac4microbiome::reshapeData(df, sampleInfo$subject, taxonomy, sampleInfo$Days)

# Repair modes
colnames(faeces$mode1) = c("subject", "index")
faeces$mode1 = faeces$mode1 %>% left_join(subjectMeta) %>% select(-index) %>% left_join(subjectMeta_faeces)
faeces$mode2 = faeces$mode2 %>% select(-index)
faeces$mode3 = faeces$mode3 %>% mutate(visit=1:3, timepoint=timepointMetadata) %>% select(-timepointMetadata,-index)

# Old approach:
#
# Filter based on sparsity
# threshold = 0.75
# sparsity = colSums(df==0) / nrow(df)
# featureSelection = sparsity <= threshold

# # CLR
# df_clr = t(apply(df+1, 1, function(x){log(x / compositions::geometricmean(x))})) %>% as_tibble()
#
# # Apply feature filter
# df_clr = df_clr[,featureSelection]
# taxonomy_filtered = taxonomy[featureSelection,]

# Make into cube
# I = length(unique(sampleInfo$subject))
# J = ncol(df)
# K = length(unique(sampleInfo$Days))
# X = array(0L, c(I,J,K))
# timepoints = sampleInfo %>% arrange(Days) %>% select(Days) %>% unique() %>% pull()
#
# for(k in 1:K){
#   Day = timepoints[k]
#   X[,,k] = cbind(df, sampleInfo) %>%
#     as_tibble() %>%
#     mutate(subject=as.character(subject)) %>%
#     filter(Days == Day) %>%
#     select(c(colnames(df),subject)) %>%
#     right_join(subjectMeta) %>%
#     arrange(subject) %>%
#     select(-colnames(subjectMeta)) %>%
#     as.matrix()
# }
#
#
# Mask based on shared subjects for BMI and WHZ
# X_bmi = X[subjectMeta$subject %in% homogenized_subjectMeta_bmi$subject,,]
# X_whz = X[subjectMeta$subject %in% homogenized_subjectMeta_whz$subject,,]
# X_bmi = X
# X_whz = X
#
# # Center and scale
# X_bmi_cnt = parafac4microbiome::multiwayCenter(X_bmi, mode=1)
# X_bmi_cnt_scl = parafac4microbiome::multiwayScale(X_bmi_cnt, mode=2)
#
# X_whz_cnt = parafac4microbiome::multiwayCenter(X_whz, mode=1)
# X_whz_cnt_scl = parafac4microbiome::multiwayScale(X_whz_cnt, mode=2)
#
# # Save
# faeces_df_bmi = X_bmi_cnt_scl
# faeces_df_whz = X_whz_cnt_scl
# faeces_subjectMeta = subjectMeta
# faeces_taxonomy = taxonomy_filtered
# faeces_timepoints = matrix(NA, nrow=3, ncol=2) %>% as_tibble() %>% mutate(visit = 1:3, timepoint=timepoints) %>% select(-V1,-V2)

# Milk microbiome
df = read.csv("./data-raw/Jakobsen2025/milkCounts.csv", header=FALSE, sep=" ") %>% as_tibble()
taxonomy = read.csv("./data-raw/Jakobsen2025/newTaxonomy_milk.csv", header=FALSE, sep=" ") %>% as_tibble()
sampleInfo = read.csv("./data-raw/Jakobsen2025/milk_sampleMeta.csv", header=FALSE, sep=" ") %>% as_tibble()
colnames(sampleInfo) = c("Sample", "RCID", "BMI", "BMI.group", "Days", "Gestational.age", "C.section", "AB.infant", "AB.mother", "Secretor", "Lewis", "subject")

# Make subject metadata
subjectMeta = sampleInfo %>% select(subject, BMI, BMI.group, C.section, Secretor, Lewis) %>% unique() %>% mutate(subject=as.character(subject)) %>% arrange(subject)

# Filter taxa to taxa with at least 1 non-zero value
featureMask = colSums(df) > 0
df = df[,featureMask]
taxonomy = taxonomy[featureMask,]

# New approach per 20250402
milk = parafac4microbiome::reshapeData(df, sampleInfo$subject, taxonomy, sampleInfo$Days)

# Repair modes
colnames(milk$mode1) = c("subject", "index")
milk$mode1 = milk$mode1 %>% left_join(subjectMeta) %>% select(-index) %>% left_join(subjectMeta_milk)
milk$mode2 = milk$mode2 %>% select(-index)
milk$mode3 = milk$mode3 %>% mutate(visit=1:4, timepoint=timepointMetadata) %>% select(-timepointMetadata,-index)

# Old approach:
#
# # Filter based on sparsity
# threshold = 0.85
# sparsity = colSums(df==0) / nrow(df)
# featureSelection = sparsity <= threshold
#
# # CLR
# df_clr = t(apply(df+1, 1, function(x){log(x / compositions::geometricmean(x))})) %>% as_tibble()
#
# # Apply feature filter
# df_clr = df_clr[,featureSelection]
# taxonomy_filtered = taxonomy[featureSelection,]
#
# # Make into cube
# I = length(unique(sampleInfo$subject))
# J = ncol(df_clr)
# K = length(unique(sampleInfo$Days))
# X = array(0L, c(I,J,K))
# timepoints = sampleInfo %>% arrange(Days) %>% select(Days) %>% unique() %>% pull()
#
# for(k in 1:K){
#   Day = timepoints[k]
#   X[,,k] = cbind(df_clr, sampleInfo) %>% as_tibble() %>% mutate(subject=as.character(subject)) %>% filter(Days == Day) %>% select(c(colnames(df_clr),subject)) %>% right_join(subjectMeta) %>% arrange(subject) %>% select(-colnames(subjectMeta)) %>% as.matrix()
# }

# Mask based on shared subjects for BMI and WHZ
# X_bmi = X[subjectMeta$subject %in% homogenized_subjectMeta_bmi$subject,,]
# X_whz = X[subjectMeta$subject %in% homogenized_subjectMeta_whz$subject,,]
# X_bmi = X
# X_whz = X
#
# # Center and scale
# X_bmi_cnt = parafac4microbiome::multiwayCenter(X_bmi, mode=1)
# X_bmi_cnt_scl = parafac4microbiome::multiwayScale(X_bmi_cnt, mode=2)
#
# X_whz_cnt = parafac4microbiome::multiwayCenter(X_whz, mode=1)
# X_whz_cnt_scl = parafac4microbiome::multiwayScale(X_whz_cnt, mode=2)
#
# milk_df_bmi = X_bmi_cnt_scl
# milk_df_whz = X_whz_cnt_scl
# milk_subjectMeta = subjectMeta
# milk_taxonomy = taxonomy_filtered
# milk_timepoints = matrix(NA, nrow=4, ncol=2) %>% as_tibble() %>% mutate(visit = 1:4, timepoint=timepoints) %>% select(-V1,-V2)

# Milk metabolomics
df = read.csv("./data-raw/Jakobsen2025/milkMetabNumeric.csv", header=FALSE, sep=" ") %>% as_tibble()
taxonomy = read.csv("./data-raw/Jakobsen2025/milk_metab_CAS_numbers.csv", header=TRUE, sep=",") %>% as_tibble()
sampleInfo = read.csv("./data-raw/Jakobsen2025/milkMetab_sampleMeta.csv", header=FALSE, sep=" ") %>% as_tibble()
colnames(sampleInfo) = c("RCID", "BMI", "BMI.group", "Days", "Gestational.age", "C.section", "AB.infant", "AB.mother", "Secretor", "Lewis", "subject")
metaboliteCategories = read.csv("./data-raw/Jakobsen2025/Metabolite_categories.txt", sep="\t") %>% as_tibble()
metaboliteCategories = metaboliteCategories %>% mutate(Metabolite = vctrs::vec_as_names(Metabolite, repair="universal_quiet"))

# Make subject metadata
subjectMeta = sampleInfo %>% select(subject, BMI, BMI.group, C.section, Secretor, Lewis) %>% unique() %>% mutate(subject=as.character(subject)) %>% arrange(subject)

# Remove duplicates
drop = c(61,81,146)
subjectMeta = subjectMeta %>% mutate(index=1:nrow(.)) %>% filter(!index %in% drop) %>% select(-index) %>% arrange(subject)

# Log transform
df_log = log(df)

# Make into cube - new approach per 20250402
milkMetab = parafac4microbiome::reshapeData(df_log, sampleInfo$subject, taxonomy, sampleInfo$Days)

# Edit metadata per mode to avoid breaking changes
colnames(milkMetab$mode1) = c("subject", "index")
milkMetab$mode1 = milkMetab$mode1 %>% left_join(subjectMeta) %>% select(-index) %>% left_join(subjectMeta_milkMetab)

milkMetab$mode3 = milkMetab$mode3 %>% mutate(visit = index, timepoint = timepointMetadata) %>% select(-timepointMetadata,-index)

# Fix milkMetab feature annotations
df = milkMetab$mode2 %>% mutate(Metabolite = make.names(X), CAS.Registry = make.names(CAS.Registry)) %>% left_join(metaboliteCategories)

# Fix mismatches by hand
df[df$Metabolite == "X2.Aminobutyrate","Class"] = "Amino acids and derivatives"
df[df$Metabolite == "X2.Fucosyllactose","Class"] = "Oligosaccharides"
df[df$Metabolite == "X2.Hydroxybutyrate","Class"] = "Amino acids and derivatives"
df[df$Metabolite == "X2.Oxoglutarate","Class"] = "Energy related"
df[df$Metabolite == "X3.Fucosyllactose","Class"] = "Oligosaccharides"
df[df$Metabolite == "X3SL.partial","Class"] = "Oligosaccharides"
df[df$Metabolite == "X6SL.partial","Class"] = "Oligosaccharides"
df[df$Metabolite == "Methionine","Class"] = "Amino acids and derivatives"
df[70,"Metabolite"] = "tau.Methylhistidine" # Fix non-ascii tau character

# Repair metabolite names
df[df$Metabolite == "X2.Aminobutyrate","Metabolite"] = "2-Aminobutyrate"
df[df$Metabolite == "X2.Fucosyllactose","Metabolite"] = "2-Fucosyllactose"
df[df$Metabolite == "X2.Hydroxybutyrate","Metabolite"] = "2-Hydroxybutyrate"
df[df$Metabolite == "X2.Oxoglutarate","Metabolite"] = "2-Oxoglutarate"
df[df$Metabolite == "X3.Fucosyllactose","Metabolite"] = "3-Fucosyllactose"
df[df$Metabolite == "X3SL.partial","Metabolite"] = "3SL partial"
df[df$Metabolite == "X6SL.partial","Metabolite"] = "6SL partial"
df[df$Metabolite == "Unk..HMO.1.70.Sialylated","Metabolite"] = "Unk. HMO 1.70 Sialylated"
df[df$Metabolite == "Unk..HMO.5.076","Metabolite"] = "Unk. HMO 5.076"
df[df$Metabolite == "Unk..HMO.5.09","Metabolite"] = "Unk. HMO 5.09"
df[df$Metabolite == "Unk..HMO.5.10","Metabolite"] = "Unk. HMO 5.10"
df[df$Metabolite == "Unk..HMO.5.11","Metabolite"] = "Unk. HMO 5.11"
df[df$Metabolite == "Unk..HMO.5.29","Metabolite"] = "Unk. HMO 5.29"
df[df$Metabolite == "Unk..HMO.5.359","Metabolite"] = "Unk. HMO 5.359"
df[df$Metabolite == "sn.Glycero.3.phosphocholine","Metabolite"] = "sn-Glycero-3-phosphocholine"
df[df$Metabolite == "myo.Inositol","Metabolite"] = "myo-Inositol"
df[df$Metabolite == "cis.Aconitate","Metabolite"] = "cis-Aconitate"
df[df$Metabolite == "DSLNT.partial","Metabolite"] = "DSLNT partial"
df[df$Metabolite == "O.Phosphocholine","Metabolite"] = "O-Phosphocholine"
df[df$Metabolite == "O.Acetylcarnitine","Metabolite"] = "O-Acetylcarnitine"

milkMetab$mode2 = df %>% select(-X,-index)
# Old approach:
#
# I = length(unique(sampleInfo$subject))
# J = ncol(df_log)
# K = length(unique(sampleInfo$Days))
# X = array(0L, c(I,J,K))
# timepoints = sampleInfo %>% arrange(Days) %>% select(Days) %>% unique() %>% pull()
#
# for(k in 1:K){
#   Day = timepoints[k]
#   X[,,k] = cbind(df_log, sampleInfo) %>% as_tibble() %>% mutate(subject=as.character(subject)) %>% filter(Days == Day) %>% select(c(colnames(df_log),subject)) %>% right_join(subjectMeta) %>% arrange(subject) %>% select(-colnames(subjectMeta)) %>% as.matrix()
# }
#
# Remove problematic sample
# mask = !is.na(subjectMeta$BMI)
# X = X[mask,,]
# subjectMeta = subjectMeta[mask,]
# #
# # Mask based on shared subjects for BMI and WHZ
# # X_bmi = X[subjectMeta$subject %in% homogenized_subjectMeta_bmi$subject,,]
# # X_whz = X[subjectMeta$subject %in% homogenized_subjectMeta_whz$subject,,]
# X_bmi = X
# X_whz = X
#
# # Center and scale
# X_bmi_cnt = parafac4microbiome::multiwayCenter(X_bmi, mode=1)
# X_bmi_cnt_scl = parafac4microbiome::multiwayScale(X_bmi_cnt, mode=2)
#
# X_whz_cnt = parafac4microbiome::multiwayCenter(X_whz, mode=1)
# X_whz_cnt_scl = parafac4microbiome::multiwayScale(X_whz_cnt, mode=2)
#
# milkMetab_df_bmi = X_bmi_cnt_scl
# milkMetab_df_whz = X_whz_cnt_scl
# milkMetab_subjectMeta = subjectMeta %>% left_join(subjectMeta_milkMetab)
# milkMetab_taxonomy = taxonomy
# milkMetab_timepoints = matrix(NA, nrow=4, ncol=2) %>% as_tibble() %>% mutate(visit = 1:4, timepoint=sampleInfo %>% arrange(Days) %>% select(Days) %>% unique() %>% pull()) %>% select(-V1,-V2)
#
# # Fix milkMetab feature annotations
# df = milkMetab_taxonomy %>% mutate(Metabolite = make.names(X), CAS.Registry = make.names(CAS.Registry)) %>% left_join(metaboliteCategories)
#
# # Fix mismatches by hand
# df[df$Metabolite == "X2.Aminobutyrate","Class"] = "Amino acids and derivatives"
# df[df$Metabolite == "X2.Fucosyllactose","Class"] = "Oligosaccharides"
# df[df$Metabolite == "X2.Hydroxybutyrate","Class"] = "Amino acids and derivatives"
# df[df$Metabolite == "X2.Oxoglutarate","Class"] = "Energy related"
# df[df$Metabolite == "X3.Fucosyllactose","Class"] = "Oligosaccharides"
# df[df$Metabolite == "X3SL.partial","Class"] = "Oligosaccharides"
# df[df$Metabolite == "X6SL.partial","Class"] = "Oligosaccharides"
# df[df$Metabolite == "Methionine","Class"] = "Amino acids and derivatives"
# df[70,"Metabolite"] = "tau.Methylhistidine" # Fix non-ascii tau character
#
# milkMetab_featureMeta = df %>% select(-X)
#
# # Save
# milkMetab = list("data"=milkMetab_df_bmi, "mode1"=milkMetab_subjectMeta, "mode2"=milkMetab_featureMeta, "mode3"=milkMetab_timepoints)

Jakobsen2025 = list()
Jakobsen2025$faeces = faeces
Jakobsen2025$milkMicrobiome = milk
Jakobsen2025$milkMetabolomics = milkMetab

usethis::use_data(Jakobsen2025, overwrite = TRUE)
