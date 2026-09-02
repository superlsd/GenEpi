library(httr)
library(jsonlite)
library(dplyr)
library(readr)

setwd("/Users/elmadervic/Documents/Github/GenEpi")

library(dplyr)
library(readr)

setwd("/Users/elmadervic/Documents/Github/GenEpi")

# ---- Read the saved output back in ----
male_long <- read.csv("NewFilesV2/male_trj_mol_drug_v2026_official_Trials.csv")

str(male_long)
summary(male_long$result_trials_count)

# ---- How many NAs, and where ----
n_na <- sum(is.na(male_long$result_trials_count))
cat("NA rows:", n_na, "/", nrow(male_long), "\n")

na_rows <- male_long %>% filter(is.na(result_trials_count))
print(na_rows)

# ---- Unique disease-drug pairs behind those NA rows ----
# (this is what you actually want to retry -- not every NA *row*,
# since many rows share the same pair across trajectories)
na_pairs <- na_rows %>% distinct(disease_member, Drug)
cat("Unique disease-drug pairs to retry:", nrow(na_pairs), "\n")
print(na_pairs)

# If you saved query_status earlier, this tells you how many were
# genuine "failed" queries vs anything else
if ("query_status" %in% names(male_long)) {
  print(table(male_long$query_status, useNA = "ifany"))
}

write.csv(na_pairs, "NewFilesV2/male_na_pairs_to_retry.csv", row.names = FALSE)


# ============================================================
# 0. get_trials_count() -- same function as before, in case this
#    is a fresh R session and it's not already defined
# ============================================================
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

# ============================================================
# 1. LOAD THE NA PAIRS TO RETRY
# ============================================================
# From the previous step (find_na_trials.R), which wrote out
# just the unique disease-drug pairs behind the NA rows.
na_pairs <- read.csv("NewFilesV2/male_na_pairs_to_retry.csv")
cat("Pairs to retry:", nrow(na_pairs), "\n")

# ============================================================
# 2. RE-QUERY JUST THESE PAIRS
# ============================================================
na_pairs$result_trials_count <- NA_integer_
na_pairs$query_status <- NA_character_

for (i in seq_len(nrow(na_pairs))) {
  val <- tryCatch(
    get_trials_count(
      disease = na_pairs$disease_member[i],
      drug    = na_pairs$Drug[i]
    ),
    error = function(e) NA_integer_
  )
  na_pairs$result_trials_count[i] <- val
  na_pairs$query_status[i] <- if (is.na(val)) "failed" else "ok"
  
  if (i %% 25 == 0) {
    cat("Retried", i, "/", nrow(na_pairs), "\n")
    write.csv(na_pairs, "NewFilesV2/male_na_pairs_retry_checkpoint.csv", row.names = FALSE)
  }
}

cat("Still NA after retry:", sum(is.na(na_pairs$result_trials_count)), "\n")
print(na_pairs %>% filter(is.na(result_trials_count)))

write.csv(na_pairs, "NewFilesV2/male_na_pairs_retried_final.csv", row.names = FALSE)

# ============================================================
# 3. MERGE RETRIED VALUES BACK INTO THE FULL DATASET
# ============================================================
male_long <- read.csv("NewFilesV2/male_trj_mol_drug_v2026_official_Trials.csv")

# Use the retried value only where the original was NA -- leave every
# already-successful row untouched.
male_long <- male_long %>%
  left_join(
    na_pairs %>% select(disease_member, Drug, retried_count = result_trials_count),
    by = c("disease_member", "Drug")
  ) %>%
  mutate(
    result_trials_count = ifelse(is.na(result_trials_count), retried_count, result_trials_count)
  ) %>%
  select(-retried_count)

cat("Remaining NAs after merge:", sum(is.na(male_long$result_trials_count)), "\n")

write.csv(male_long, "NewFilesV2/male_trj_mol_drug_v2026_official_Trials.csv", row.names = FALSE)