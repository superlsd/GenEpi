library(httr)
library(jsonlite)
library(dplyr)
library(readr)

setwd("/Users/elmadervic/Documents/Github/GenEpi")

# ============================================================
# 0. Helpers
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

# ---- Query all matching studies' phases for one disease-drug pair ----
# Returns a single string like "PHASE1; PHASE2; PHASE3" (unique, sorted),
# or NA_character_ if nothing usable was found.
get_trial_stages <- function(disease, drug, sleep_sec = 1.2,
                             page_size = 100, max_pages = 20) {
  disease <- safe_phrase(disease)
  drug    <- safe_phrase(drug)
  
  if (disease == "" || drug == "") return(NA_character_)
  
  url <- "https://clinicaltrials.gov/api/v2/studies"
  
  all_phases <- character(0)
  page_token <- NULL
  page <- 0
  
  repeat {
    page <- page + 1
    if (page > max_pages) break  # safety cap
    
    q <- list(
      "query.cond" = disease,
      "query.intr" = drug,
      "countTotal" = "true",
      "pageSize"   = page_size,
      "fields"     = "protocolSection.designModule.phases"
    )
    if (!is.null(page_token)) q[["pageToken"]] <- page_token
    
    res <- tryCatch({
      GET(url, query = q, timeout(6))
    }, error = function(e) NULL)
    
    if (sleep_sec > 0) Sys.sleep(sleep_sec)
    
    if (is.null(res) || status_code(res) != 200) break
    
    parsed <- tryCatch(
      fromJSON(content(res, "text", encoding = "UTF-8"), flatten = TRUE),
      error = function(e) NULL
    )
    if (is.null(parsed)) break
    
    studies <- parsed$studies
    if (is.null(studies) || length(studies) == 0) break
    
    # Depending on flattening, the column may be named exactly this:
    phase_col <- "protocolSection.designModule.phases"
    
    if (phase_col %in% names(studies)) {
      phases_this_page <- studies[[phase_col]]
      # phases_this_page is a list-column (one list per study, since a
      # study can have multiple phases, e.g. "PHASE2","PHASE3")
      phases_flat <- unlist(phases_this_page, use.names = FALSE)
      all_phases <- c(all_phases, phases_flat)
    }
    
    page_token <- parsed$nextPageToken %||% NULL
    if (is.null(page_token)) break
  }
  
  all_phases <- all_phases[all_phases != "" & !is.na(all_phases)]
  
  if (length(all_phases) == 0) return(NA_character_)
  
  paste(all_phases, collapse = "; ")
}

# ============================================================
# 1. Read the resulting file
# ============================================================
male_long <- read.csv("NewFilesV2/male_trj_mol_drug_v2026_official_Trials.csv")

str(male_long)

male_long$stages.x <- NULL
male_long$stages.y <- NULL
cat("Rows with >0 trials:", sum(male_long$result_trials_count > 0, na.rm = TRUE), "\n")

# ============================================================
# 2. Unique disease-drug pairs with >0 trials
#    (query once per pair, not once per row, to save API calls)
# ============================================================
pairs_to_query <- male_long %>%
  filter(!is.na(result_trials_count), result_trials_count > 0) %>%
  distinct(disease_member, Drug)

cat("Unique disease-drug pairs to query for stages:", nrow(pairs_to_query), "\n")

pairs_to_query$stages <- NA_character_

for (i in seq_len(nrow(pairs_to_query))) {
  pairs_to_query$stages[i] <- tryCatch(
    get_trial_stages(
      disease = pairs_to_query$disease_member[i],
      drug    = pairs_to_query$Drug[i]
    ),
    error = function(e) NA_character_
  )
  
  if (i %% 250 == 0) {
    cat("Queried", i, "/", nrow(pairs_to_query), "\n")
    write.csv(pairs_to_query, "NewFilesV2/male_stages_checkpoint.csv", row.names = FALSE)
  }
}

cat("Pairs with no stage info found:", sum(is.na(pairs_to_query$stages)), "\n")

write.csv(pairs_to_query, "NewFilesV2/male_pairs_stages_final.csv", row.names = FALSE)

# ============================================================
# 3. Merge stages back into the full dataset
# ============================================================
male_long <- male_long %>%
  left_join(pairs_to_query %>% select(disease_member, Drug, stages),
            by = c("disease_member", "Drug"))

write.csv(male_long, "NewFilesV2/male_trj_mol_drug_v2026_official_Trials.csv", row.names = FALSE)



summary(male_long)

retry_pairs <- pairs_to_query %>%
  filter(is.na(stages)) %>%
  left_join(male_long %>% distinct(disease_member, Drug, result_trials_count),
            by = c("disease_member", "Drug")) %>%
  filter(result_trials_count >= 2)

nrow(retry_pairs)  # should be ~282

for (i in seq_len(nrow(retry_pairs))) {
  new_stage <- tryCatch(
    get_trial_stages(disease = retry_pairs$disease_member[i], drug = retry_pairs$Drug[i]),
    error = function(e) NA_character_
  )
  if (!is.na(new_stage)) {
    pairs_to_query$stages[pairs_to_query$disease_member == retry_pairs$disease_member[i] &
                            pairs_to_query$Drug == retry_pairs$Drug[i]] <- new_stage
  }
  if (i %% 25 == 0) cat("Retried", i, "/", nrow(retry_pairs), "\n")
}

write.csv(pairs_to_query, "NewFilesV2/male_pairs_stages_final.csv", row.names = FALSE)


summary(pairs_to_query)
sum(is.na(pairs_to_query$stages))