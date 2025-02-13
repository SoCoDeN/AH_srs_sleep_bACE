
# CODE TO PREPARE VCU AND ABCD DATA FOR TWIN ANALYSIS
# Last updated 2/13/2025 by Amanda Halliday

rm(list = ls())
library(tidyverse)
library(OpenMx)

# ------------------------------------------------------------------------------------------------------------------
# READ IN DATA
vcu_long <- read_csv("/Volumes/SoCoDeN/SSCDN_Private/data/VCU_MATR/Amanda_analysis/github/vcu_long.csv")
  # Data previously organized to include variables of interest, recode race, score SRS and SDS, and create "long" format
  # Data dictionary:
    # FID = Family ID = unique identifier for each family group
    # TID = Twin ID = unique identifier for each twin pair within a family group
    # zygo = zygosity as determined by MATR
    # age = age in years
    # sex = 1=males, 2=females
    # race is categorized to match ABCD
      # isCA = caucasian
      # isAA = african american
      # isHI = hispanic
      # isAS = asian
      # isMI = mixed
      # isOT = other/missing
    # srs = 65 item Social Responsiveness Scale (SRS) raw total score
    # ssrs = 11 item Social Responsiveness Scale (SRS) raw score
    # sds = 26 item Sleep Disturbance Scale for Children (SDS) raw total score
abcd_long <- read_csv("/Volumes/SoCoDeN/SSCDN_Private/data/VCU_MATR/Amanda_analysis/github/abcd_long.csv")
  # Data previously organized to separate twins, include variables of interest, recode race, and recode site
  # Data dictionary:
    # IID = individual ID (src_subject_id) = unique indentifier for each individual participant
    # FID = Family ID = unique identifier for each family group
    # TID = Twin ID = unique identifier for each twin pair within a family group
    # zygo = zygosity as determined by ABCD, not genetic
    # ssrs = 11 item Social Responsiveness Scale (SRS) raw score
    # sds = 26 item Sleep Disturbance Scale for Children (SDS) raw total score
    # age = age in years
    # sex = 1=males, 2=females
    # race is categorized to match ABCD
      # isCA = caucasian
      # isAA = african american
      # isHI = hispanic
      # isAS = asian
      # isMI = mixed
      # isOT = other/missing
    # site
      # isCUB = Colorado Boulder
      # isUMN = Minnesota
      # isVCU = Virginia Commonwealth University *ABCD participants were EXCLUDED from participating in SoCoDeN data
      # isWSL = Washington University St. Louis
      # isOTH = other ABCD site


# ------------------------------------------------------------------------------------------------------------------
# REMOVING NAs in variables of interest
# Change "VCU" to "ABCD" in data frames to prep ABCD data. Code is otherwise the same for both cohorts.

cols_to_check <- c("srs", "ssrs", "sds")

# Step 1: Count the number of NAs per variable before removal
na_counts_before <- sapply(vcu_long[cols_to_check], function(x) sum(is.na(x)))

# Step 2: Remove rows with NAs in any of the specified columns
vcu_long_clean <- vcu_long %>%
  filter(complete.cases(vcu_long[cols_to_check])) # Removes rows where any NA exists in the specified columns

# Step 3: Count the number of NAs per variable after removal
na_counts_after <- sapply(vcu_long_clean[cols_to_check], function(x) sum(is.na(x)))

# Step 4: Count the number of rows removed
rows_removed <- nrow(vcu_long) - nrow(vcu_long_clean)

# Step 5: Calculate how many NAs were removed per variable
na_removed_per_variable <- na_counts_before - na_counts_after

# Display results
cat("Number of rows removed:", rows_removed, "\n")
cat("NA counts removed per variable:\n")
print(na_removed_per_variable)

# Removing solo pairs now that NAs have been removed
# Step 1: Count the number of rows per twin pair
pair_counts <- vcu_long_clean %>%
  group_by(TID) %>%
  summarise(pair_count = n()) # Count rows per pair

# Step 2: Identify complete pairs (pair_count == 2)
complete_pairs <- pair_counts %>%
  filter(pair_count == 2) %>%
  pull(TID) # Extract pair IDs with 2 twins

# Step 3: Filter the dataset to retain only complete pairs
vcu_long_complete <- vcu_long_clean %>%
  filter(TID %in% complete_pairs)

# Check the number of rows removed
rows_removed <- nrow(vcu_long_clean) - nrow(vcu_long_complete)
cat("Number of rows removed:", rows_removed)

vcu25l <- vcu_long_complete


# ------------------------------------------------------------------------------------------------------------------
# RESIDUALIZE DATA

vcu25l <- vcu25l %>%
  mutate(IID = row_number()) # give every twin an individual ID; only necessary for VCU

# Identify variables of interest
vcuVars <- names(vcu25l[, c(12:14)])  # Variables for regression

# Perform regression and store residuals
for (i in vcuVars) {
  # Fit linear model
  temp <- lm(vcu25l[[i]] ~ age + isCA + isAA + isHI + isAS + isMI + isOT, 
             data = vcu25l, 
             na.action = na.exclude)
  
  # Store residuals in a new column
  vcu25l[[paste0('r', i)]] <- residuals(temp)
}


# ------------------------------------------------------------------------------------------------------------------
# REMOVE OUTLIERS

# Function to standardize variables
myStan <- function(x) {
  (x - mean(x, na.rm = TRUE)) / sd(x, na.rm = TRUE)
}

# Select the variables to standardize
rVars <- names(vcu25l)[16:18] # Standardizing residualized data

# Standardize and handle outliers
for (var in rVars) {
  # Standardize the column
  standardized <- myStan(vcu25l[[var]])
  
  # Replace outliers (values beyond ±4 SD) with NA
  standardized[standardized < -4 | standardized > 4] <- NA
  
  # Assign back to the corresponding column in vcu25l
  vcu25l[[var]] <- standardized
}

# Count the number of NAs per column
na_counts_per_column <- colSums(is.na(vcu25l))
print(na_counts_per_column)

# Remove rows with any NA values (outliers)
vcu25lrs <- na.omit(vcu25l) #long, residualized, standardized


# ------------------------------------------------------------------------------------------------------------------
# RANDOMIZE TWINS AND MAKE WIDE DATA FRAME (each family has unique row)

# Identify unmatched TIDs (those that appear only once). These are now unmatched because their twin had an outlier score
unmatched_TIDs <- vcu25lrs %>%
  group_by(TID) %>%
  filter(n() == 1) %>%
  ungroup() %>%
  select(TID) %>%
  distinct()

# Remove rows where TID has no pair (i.e., those in unmatched_TIDs)
vcu25lrso <- vcu25lrs %>% # NAs and outliers removed
  filter(!TID %in% unmatched_TIDs$TID)

# Check that every twin has a pair (the result should be empty df)
vcu25lrso %>%
  group_by(TID) %>%
  summarise(n = n()) %>%
  filter(n != 2)

set.seed(123)  # Set the seed to 123
vcu25lrso <- vcu25lrso %>%
  group_by(TID) %>%
  mutate(ah_id_new = sample(c("1", "2"))) %>%
  ungroup()

# Check that each pair has a T1 and a T2 (the result should be empty df)
vcu25lrso %>%
  group_by(TID, ah_id_new) %>%
  summarise(n = n(), .groups = "drop") %>%
  filter(n != 1)

# Make separate dfs for T1 and T2
vcu25_T1 <- vcu25lrso %>%
  filter(ah_id_new == "1") %>%
  select(TID:zygo, sex, srs:ah_id_new) %>%
  rename(
    IID_T1 = IID,
    sex_T1 = sex,
    srs_T1 = srs,
    sds_T1 = sds,
    ssrs_T1 = ssrs,
    rsrs_T1 = rsrs,
    rssrs_T1 = rssrs,
    rsds_T1 = rsds)

vcu25_T2 <- vcu25lrso %>%
  filter(ah_id_new == "2") %>%
  select(TID:zygo, sex, srs:ah_id_new) %>%
  rename(
    IID_T2 = IID,
    sex_T2 = sex,
    srs_T2 = srs,
    sds_T2 = sds,
    ssrs_T2 = ssrs,
    rsrs_T2 = rsrs,
    rssrs_T2 = rssrs,
    rsds_T2 = rsds)

# Putting dfs together
vcu25w <- merge(vcu25_T1, vcu25_T2, by = "TID") %>%
  rename(zygo = zygo.x) %>%
  select(-zygo.y, -ah_id_new.x, -ah_id_new.y)