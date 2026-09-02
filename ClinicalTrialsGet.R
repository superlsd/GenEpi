library(httr)
library(jsonlite)
library(dplyr)
library(purrr)
library(stringr)
library(tibble)
library(readr)
library(tidyr)

setwd("/Users/elmadervic/Documents/Github/GenEpi")

`%||%` <- function(a, b) {
  if (is.null(a)) return(b)
  if (length(a) == 0) return(b)
  a
}

safe_phrase <- function(x) {
  x <- as.character(x)
  x[is.na(x)] <- ""
  x <- trimws(x)
  x
}

# ------------------------------------------------------------------
# ClinicalTrials.gov API v2 query function
# Returns total number of registered trials matching disease + drug.
# ------------------------------------------------------------------
get_trials_count <- function(disease, drug, sleep_sec = 1.2) {
  disease <- safe_phrase(disease)
  drug    <- safe_phrase(drug)
  
  if (disease == "" || drug == "") return(NA_integer_)
  
  url <- "https://clinicaltrials.gov/api/v2/studies"
  
  res <- tryCatch({
    GET(
      url,
      query = list(
        "query.cond" = disease,
        "query.intr" = drug,
        "countTotal" = "true",
        "pageSize"   = 1
      ),
      timeout(6)
    )
  }, error = function(e) NULL)
  
  if (sleep_sec > 0) Sys.sleep(sleep_sec)
  
  if (is.null(res) || status_code(res) != 200) return(NA_integer_)
  
  parsed <- tryCatch(
    fromJSON(content(res, "text", encoding = "UTF-8")),
    error = function(e) NULL
  )
  if (is.null(parsed)) return(NA_integer_)
  
  as.integer(parsed$totalCount %||% NA_integer_)
}

# ------------------------------------------------------------------
# Optional: pull a status breakdown (RECRUITING / COMPLETED / etc.)
# for a disease-drug pair, in addition to the raw count.
# Set INCLUDE_STATUS_BREAKDOWN <- TRUE below to use it.
# ------------------------------------------------------------------
get_trials_status_breakdown <- function(disease, drug, sleep_sec = 1.2) {
  disease <- safe_phrase(disease)
  drug    <- safe_phrase(drug)
  
  empty_result <- tibble(
    n_recruiting = NA_integer_,
    n_completed  = NA_integer_,
    n_terminated = NA_integer_,
    n_withdrawn  = NA_integer_,
    n_active     = NA_integer_,
    n_other      = NA_integer_
  )
  
  if (disease == "" || drug == "") return(empty_result)
  
  url <- "https://clinicaltrials.gov/api/v2/studies"
  
  res <- tryCatch({
    GET(
      url,
      query = list(
        "query.cond" = disease,
        "query.intr" = drug,
        "fields"     = "OverallStatus",
        "pageSize"   = 1000
      ),
      timeout(6)
    )
  }, error = function(e) NULL)
  
  if (sleep_sec > 0) Sys.sleep(sleep_sec)
  
  if (is.null(res) || status_code(res) != 200) return(empty_result)
  
  parsed <- tryCatch(
    fromJSON(content(res, "text", encoding = "UTF-8"), simplifyVector = TRUE),
    error = function(e) NULL
  )
  if (is.null(parsed) || is.null(parsed$studies) || length(parsed$studies) == 0) {
    return(empty_result)
  }
  
  statuses <- tryCatch(
    parsed$studies$protocolSection$statusModule$overallStatus,
    error = function(e) NULL
  )
  if (is.null(statuses)) return(empty_result)
  
  tibble(
    n_recruiting = sum(statuses == "RECRUITING", na.rm = TRUE),
    n_completed  = sum(statuses == "COMPLETED", na.rm = TRUE),
    n_terminated = sum(statuses == "TERMINATED", na.rm = TRUE),
    n_withdrawn  = sum(statuses == "WITHDRAWN", na.rm = TRUE),
    n_active     = sum(statuses %in% c("ACTIVE_NOT_RECRUITING", "ENROLLING_BY_INVITATION"), na.rm = TRUE),
    n_other      = sum(!statuses %in% c("RECRUITING", "COMPLETED", "TERMINATED",
                                        "WITHDRAWN", "ACTIVE_NOT_RECRUITING",
                                        "ENROLLING_BY_INVITATION"), na.rm = TRUE)
  )
}

# Set TRUE if you also want the status breakdown columns (slower: 1 extra call/row)
INCLUDE_STATUS_BREAKDOWN <- FALSE

############################################################
#
########## MALE ##########
#
############################################################

female_data <- read.csv("NewFilesV2/random_male_trj_mol_drug_v2026_official.csv")
str(female_data)

female_long <- female_data %>%
  mutate(
    disease_member = str_split(
      str_replace_all(`Disease.Consensus.Name.List`, "\\[|\\]|'|\"", ""),
      "\\s*,\\s*"
    )
  ) %>%
  unnest(disease_member) %>%
  mutate(disease_member = str_trim(disease_member)) %>%
  filter(disease_member != "")

str(female_long)
head(female_long)

# ---- Dedupe before querying: many rows repeat the same disease-drug pair
# across trajectories, so query unique pairs only, then join back. ----
female_unique_pairs <- female_long %>%
  distinct(disease_member, Drug)

cat("Female: ", nrow(female_long), " total rows -> ",
    nrow(female_unique_pairs), " unique disease-drug pairs to query\n", sep = "")

female_unique_pairs$result_trials_count <- NA_integer_
female_unique_pairs$query_status <- NA_character_  # "ok" / "failed" / not yet attempted
for (i in seq_len(nrow(female_unique_pairs))) {
  if (is.na(female_unique_pairs$result_trials_count[i])) {
    # tryCatch here (on top of get_trials_count's own internal tryCatch)
    # guarantees a network hiccup, timeout, or malformed response never
    # kills the whole loop -- it just records NA + "failed" and moves on.
    val <- tryCatch(
      get_trials_count(
        disease = female_unique_pairs$disease_member[i],
        drug    = female_unique_pairs$Drug[i]
      ),
      error = function(e) NA_integer_
    )
    female_unique_pairs$result_trials_count[i] <- val
    female_unique_pairs$query_status[i] <- if (is.na(val)) "failed" else "ok"
    # print(val)
  }
  # checkpoint every 10000 pairs in case of a crash or rate-limit interruption
  if (i %% 10000 == 0) {
    write.csv(female_unique_pairs, "NewFilesV2/male_unique_pairs_Trials_checkpoint.csv", row.names = F)
    n_failed <- sum(female_unique_pairs$query_status[seq_len(i)] == "failed", na.rm = TRUE)
    cat("Checkpoint saved at pair ", i, " / ", nrow(female_unique_pairs),
        " (", n_failed, " failed so far)\n", sep = "")
  }
}

if (INCLUDE_STATUS_BREAKDOWN) {
  status_list <- vector("list", nrow(female_unique_pairs))
  for (i in seq_len(nrow(female_unique_pairs))) {
    status_list[[i]] <- get_trials_status_breakdown(
      disease = female_unique_pairs$disease_member[i],
      drug    = female_unique_pairs$Drug[i]
    )
  }
  female_unique_pairs <- bind_cols(female_unique_pairs, bind_rows(status_list))
}

female_long <- female_long %>%
  left_join(female_unique_pairs, by = c("disease_member", "Drug"))


tail(female_long$result_trials_count)
summary(female_long)
sum(!is.na(female_long$result_trials_count) & female_long$result_trials_count != 0)

write.csv(female_long, "NewFilesV2/random_male_trj_mol_drug_v2026_official_Trials.csv", row.names = F)

############################################################
#
########## MALE ##########
#
############################################################

male_data <- read.csv("NewFilesV2/male_trj_mol_drug_v2026_official.csv")
str(male_data)

male_long <- male_data %>%
  mutate(
    disease_member = str_split(
      str_replace_all(`Disease.Consensus.Name.List`, "\\[|\\]|'|\"", ""),
      "\\s*,\\s*"
    )
  ) %>%
  unnest(disease_member) %>%
  mutate(disease_member = str_trim(disease_member)) %>%
  filter(disease_member != "")

str(male_long)
head(male_long)

# ---- Dedupe before querying: many rows repeat the same disease-drug pair
# across trajectories, so query unique pairs only, then join back. ----
male_unique_pairs <- male_long %>%
  distinct(disease_member, Drug)

cat("Male: ", nrow(male_long), " total rows -> ",
    nrow(male_unique_pairs), " unique disease-drug pairs to query\n", sep = "")

male_unique_pairs$result_trials_count <- NA_integer_
male_unique_pairs$query_status <- NA_character_  # "ok" / "failed" / not yet attempted
for (i in seq_len(nrow(male_unique_pairs))) {
  if (is.na(male_unique_pairs$result_trials_count[i])) {
    val <- tryCatch(
      get_trials_count(
        disease = male_unique_pairs$disease_member[i],
        drug    = male_unique_pairs$Drug[i]
      ),
      error = function(e) NA_integer_
    )
    male_unique_pairs$result_trials_count[i] <- val
    male_unique_pairs$query_status[i] <- if (is.na(val)) "failed" else "ok"
    # print(val)
  }
  # checkpoint every 10000 pairs in case of a crash or rate-limit interruption
  if (i %% 10000 == 0) {
    write.csv(male_unique_pairs, "NewFilesV2/male_unique_pairs_Trials_checkpoint.csv", row.names = F)
    n_failed <- sum(male_unique_pairs$query_status[seq_len(i)] == "failed", na.rm = TRUE)
    cat("Checkpoint saved at pair ", i, " / ", nrow(male_unique_pairs),
        " (", n_failed, " failed so far)\n", sep = "")
  }
}

if (INCLUDE_STATUS_BREAKDOWN) {
  status_list <- vector("list", nrow(male_unique_pairs))
  for (i in seq_len(nrow(male_unique_pairs))) {
    status_list[[i]] <- get_trials_status_breakdown(
      disease = male_unique_pairs$disease_member[i],
      drug    = male_unique_pairs$Drug[i]
    )
  }
  male_unique_pairs <- bind_cols(male_unique_pairs, bind_rows(status_list))
}

male_long <- male_long %>%
  left_join(male_unique_pairs, by = c("disease_member", "Drug"))


tail(male_long$result_trials_count)
summary(male_long)
sum(!is.na(male_long$result_trials_count) & male_long$result_trials_count != 0)

write.csv(male_long, "NewFilesV2/male_trj_mol_drug_v2026_official_Trials.csv", row.names = F)