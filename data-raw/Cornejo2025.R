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

# Saliva
saliva = read.csv("./data-raw/Cornejo2025/20240503_UNOISE_new/salivaCounts.csv", header=FALSE) %>% as_tibble()
saliva_taxonomy = read.csv("./data-raw/Cornejo2025/20240503_UNOISE_new/taxonomysaliva_fixed.csv", sep=" ") %>% as_tibble()
saliva_sampleMeta = read.csv("./data-raw/Cornejo2025/20240503_UNOISE_new/salivaSampleMeta.csv",header=FALSE) %>% as_tibble()

colnames(saliva) = saliva_taxonomy$zOTU
temp = read.csv("./data-raw/Cornejo2025/sampleInfo_fixed.csv", sep=" ") %>% as_tibble()
colnames(saliva_sampleMeta) = c(temp %>% select(-Description,-subject,-newTimepoint) %>% colnames, "Description", "subject", "newTimepoint")

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

# Cytokines
df = read.csv("./data-raw/Cornejo2025/20241209_cytokines.csv", header=FALSE, sep=" ") %>% as_tibble()
featureMeta = read.csv("./data-raw/Cornejo2025/20241209_cytokines_featureMeta.csv", header=FALSE) %>% as_tibble()
sampleInfo = read.csv("./data-raw/Cornejo2025/20241209_cytokines_sampleMeta.csv", header=FALSE, sep=" ") %>% as_tibble()
colnames(sampleInfo) = c("subject", "GenderID", "newTimepoint", "unknown", "unknown2")

# Put into cube
cytokine = parafac4microbiome::reshapeData(df, sampleInfo$subject, featureMeta, sampleInfo$newTimepoint)

# Repair modes to avoid breaking changes
colnames(cytokine$mode1) = c("subject", "index")
cytokine$mode1 = cytokine$mode1 %>% left_join(saliva_sampleMeta %>% mutate(subject=as.character(subject)) %>% select(subject,GenderID) %>% unique()) %>% select(-index)
cytokine$mode2 = cytokine$mode2 %>% select(-index)
cytokine$mode3 = cytokine$mode3 %>% mutate(newTimepoint=timepointMetadata) %>% select(-timepointMetadata,-index)

# Salivary biochemistry
biochemistry = saliva_sampleMeta %>% left_join(otherMeta %>% select(subject,newTimepoint,pH)) %>% select(-value, -BarcodeSequence, -LinkerPrimerSequence, -Niche, -Timepoint, -Age, -qPCR_16S_ng_ul, -FQ_ng_ul, -Participant_code, -Testosterone_nmol_L, -Free_testosterone_Vermeulen_pmol_L, -Estradiol_pmol_ml, -Description,-GenderID, -LH_U_L, -SHBG_nmol_L)
Xlong = biochemistry %>% select(-subject,-newTimepoint) %>% as.matrix()
biochemistry = parafac4microbiome::reshapeData(Xlong, biochemistry$subject, colnames(Xlong), biochemistry$newTimepoint)

# Change modes for consistency
colnames(biochemistry$mode1) = c("subject", "index")
biochemistry$mode1 = biochemistry$mode1 %>% left_join(saliva_sampleMeta %>% mutate(subject=as.character(subject)) %>% select(subject,GenderID) %>% unique()) %>% select(-index)
biochemistry$mode2 = biochemistry$mode2 %>% select(-index)
biochemistry$mode3 = biochemistry$mode3 %>% mutate(newTimepoint=timepointMetadata) %>% select(-timepointMetadata,-index)

# Blood hormone levels
df = saliva_sampleMeta %>% filter(newTimepoint != 6) %>% select(subject, newTimepoint, Free_testosterone_Vermeulen_pmol_L, Estradiol_pmol_ml, LH_U_L, SHBG_nmol_L)
blood_hormones = parafac4microbiome::reshapeData(df %>% select(-subject,-newTimepoint),
                                                 subjectMetadata = df$subject,
                                                 timepointMetadata = df$newTimepoint,
                                                 featureMetadata = c("Free_testosterone_Vermeulen_pmol_L", "Estradiol_pmol_ml", "LH_U_L", "SHBG_nmol_L"))
colnames(blood_hormones$mode1) = c("subject", "index")
blood_hormones$mode1 = blood_hormones$mode1 %>% left_join(saliva_sampleMeta %>% mutate(subject=as.character(subject)) %>% select(subject,GenderID) %>% unique()) %>% select(-index)
blood_hormones$mode2 = blood_hormones$mode2 %>% select(-index)
blood_hormones$mode3 = blood_hormones$mode3 %>% mutate(newTimepoint=timepointMetadata) %>% select(-timepointMetadata,-index)

# Clinical measurements
df = read.csv("./data-raw/Cornejo2025/CODS_XI.csv")
df = df[1:40,1:10] %>% as_tibble() %>% filter(Participant_code %in% saliva_sampleMeta$Participant_code)

CODS = df %>% select(Participant_code, CODS_BASELINE, CODS_3_MONTHS, CODS_6_MONTHS, CODS_12_MONTHS) %>% pivot_longer(-Participant_code) %>% mutate(timepoint = NA)
CODS[CODS$name == "CODS_BASELINE", "timepoint"] = 0
CODS[CODS$name == "CODS_3_MONTHS", "timepoint"] = 3
CODS[CODS$name == "CODS_6_MONTHS", "timepoint"] = 6
CODS[CODS$name == "CODS_12_MONTHS", "timepoint"] = 12
CODS[CODS$value == -99, "value"] = NA
CODS = CODS %>% mutate(subject = Participant_code, newTimepoint = timepoint, CODS = value) %>% select(-Participant_code, -name, -value, -timepoint)

XI = df %>% select(Participant_code, XI_BASELINE, XI_3_MONTHS, XI_6_MONTHS, XI_12_MONTHS) %>% pivot_longer(-Participant_code) %>% mutate(timepoint = NA)
XI[XI$name == "XI_BASELINE", "timepoint"] = 0
XI[XI$name == "XI_3_MONTHS", "timepoint"] = 3
XI[XI$name == "XI_6_MONTHS", "timepoint"] = 6
XI[XI$name == "XI_12_MONTHS", "timepoint"] = 12
XI[XI$value == -99, "value"] = NA
XI = XI %>% mutate(subject = Participant_code, newTimepoint = timepoint, XI = value) %>% select(-Participant_code, -name, -value, -timepoint)

df = otherMeta %>% left_join(CODS) %>% left_join(XI) %>% select(subject, newTimepoint, boppercent, CODS, XI) %>% filter(newTimepoint != "NA")
clinical = parafac4microbiome::reshapeData(df %>% select(-subject,-newTimepoint),
                                           subjectMeta=df$subject,
                                           featureMetadata = c("BOP", "CODS", "XI"),
                                           timepointMetadata = df$newTimepoint)
colnames(clinical$mode1) = c("subject", "index")
clinical$mode1 = clinical$mode1 %>% left_join(saliva_sampleMeta %>% mutate(subject=as.character(subject)) %>% select(subject,GenderID) %>% unique()) %>% select(-index)
clinical$mode2 = clinical$mode2 %>% select(-index)
clinical$mode3 = clinical$mode3 %>% mutate(newTimepoint=timepointMetadata) %>% select(-timepointMetadata,-index)

# Subject metadata
metadata = saliva_sampleMeta %>% select(subject, GenderID, newTimepoint, Age) %>% left_join(otherMeta) %>% select(subject, GenderID, Age) %>% unique()

# Store and save
Cornejo2025 = list()
Cornejo2025$Tongue_microbiome = tongue
Cornejo2025$Salivary_microbiome = saliva
Cornejo2025$Salivary_biochemistry = biochemistry
Cornejo2025$Salivary_cytokines = cytokine
Cornejo2025$Circulatory_hormones = blood_hormones
Cornejo2025$Clinical_measurements = clinical
Cornejo2025$Subject_metadata = metadata
usethis::use_data(Cornejo2025, overwrite = TRUE)
