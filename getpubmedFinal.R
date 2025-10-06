# --- deps ---
library(dplyr)
library(purrr)
library(stringr)
library(readr)   # only if you want to write CSV
library(rentrez)

library(rentrez)
library(dplyr)
library(purrr)
library(stringr)
library(tibble)
library(readr)


setwd("/Users/elmadervic/Documents/Github/GenEpi")

icd <- readRDS("ICD_eng.RDS")
icd <- icd[icd$Code<"N99",]

icd$Id <- NULL
icd$CodeWithSeparator <- NULL
icd$HippaCovered <- NULL
icd$Deleted <- NULL
names(icd) <- c( "icd3",     "name",     "synonyms")


icd_patch <- tribble(
  ~icd3, ~name,
  "A16", "Respiratory tuberculosis, not confirmed bacteriologically or histologically",
  "A97", "Dengue",
  "B21", "Human immunodeficiency virus [HIV] disease resulting in malignant neoplasms",
  "B22", "Human immunodeficiency virus [HIV] disease resulting in other specified diseases",
  "B23", "Human immunodeficiency virus [HIV] disease resulting in other conditions",
  "B24", "Unspecified human immunodeficiency virus [HIV] disease",
  "B98", "Other specified infectious agents as the cause of diseases classified to other chapters",
  "C97", "Malignant neoplasms of independent (primary) multiple sites",
  "E12", "Malnutrition-related diabetes mellitus",
  "E14", "Unspecified diabetes mellitus",
  "E90", "Nutritional and metabolic disorders in diseases classified elsewhere",
  "F00", "Dementia in Alzheimer disease",
  "F38", "Other mood [affective] disorders",
  "F61", "Mixed and other personality disorders",
  "F62", "Enduring personality changes, not attributable to brain damage and disease",
  "F83", "Mixed specific developmental disorders",
  "F92", "Mixed disorders of conduct and emotions",
  "G22", "Parkinsonism in diseases classified elsewhere",
  "G41", "Status epilepticus",
  "H03", "Disorders of eyelid in diseases classified elsewhere",
  "H06", "Disorders of lacrimal system and orbit in diseases classified elsewhere",
  "H13", "Disorders of conjunctiva in diseases classified elsewhere",
  "H19", "Disorders of sclera and cornea in diseases classified elsewhere",
  "H45", "Disorders of vitreous body and globe in diseases classified elsewhere",
  "H48", "Disorders of optic [2nd] nerve and visual pathways in diseases classified elsewhere",
  "H58", "Other disorders of eye and adnexa in diseases classified elsewhere",
  "I64", "Stroke, not specified as haemorrhage or infarction",
  "I84", "Haemorrhoids",
  "I98", "Other disorders of circulatory system in diseases classified elsewhere",
  "J46", "Status asthmaticus",
  "K07", "Dentofacial anomalies [including malocclusion]",
  "K10", "Other diseases of jaws",
  "K93", "Disorders of other digestive organs in diseases classified elsewhere",
  "M03", "Postinfective and reactive arthropathies in diseases classified elsewhere",
  "M09", "Juvenile arthritis in diseases classified elsewhere",
  "M68", "Disorders of synovium and tendon in diseases classified elsewhere",
  "M73", "Soft tissue disorders in diseases classified elsewhere",
  "M82", "Osteoporosis in diseases classified elsewhere",
  "N99", "Postprocedural disorders of genitourinary system, not elsewhere classified"
)

icd_patch$synonyms <- NA


icd <- rbind(icd, icd_patch)

icd <- icd[order(icd$icd3),]


prevalence <- read.csv("Prevalence_ICD.csv")
hospitalisations <- read.csv("Hospitalisation_pri_ICD.csv")
gene <- read.table("icd_characteristics_13082025.tsv", sep = "\t", header = T)


str(prevalence)
str(hospitalisations)
str(gene)

icd <- icd[icd$icd3 %in% prevalence$icd_code,]




`%||%` <- function(a,b) if (is.null(a)) b else a

# Heuristic expansions from a single name (handles diabetes variants out of the box)
expand_terms <- function(name) {
  nm <- str_squish(name)
  out <- c(nm)
  # Drop 'mellitus' (helps T2DM recall)
  out <- c(out, str_squish(gsub("\\bmellitus\\b", "", nm, ignore.case = TRUE)))
  # Roman numerals -> Arabic for "type"
  repl <- nm |>
    str_to_lower() |>
    str_replace_all("\\btype\\s*i\\b",  "type 1") |>
    str_replace_all("\\btype\\s*ii\\b", "type 2") |>
    str_replace_all("\\btype\\s*iii\\b","type 3") |>
    str_replace_all("\\btype\\s*iv\\b", "type 4")
  out <- unique(str_squish(c(out, repl)))
  out[nzchar(out)]
}

# Optional: add a MeSH anchor for very common cases (diabetes)
guess_mesh <- function(icd3, name) {
  icd3 <- toupper(icd3)
  if (icd3 == "E11") return("Diabetes Mellitus, Type 2")
  if (icd3 == "E10") return("Diabetes Mellitus, Type 1")
  if (icd3 %in% c("E13","E14")) return("Diabetes Mellitus")
  # fallback: none
  NA_character_
}

make_query <- function(icd3, name, include_bare_code = FALSE, mesh_term = NULL) {
  terms <- expand_terms(name)
  
  # For each term: unquoted TIAB, quoted TIAB, and All Fields (ATM)
  terms_name <- map_chr(
    terms,
    ~ sprintf("(%s[TIAB] OR \"%s\"[TIAB] OR %s[All Fields])", .x, .x, .x)
  )
  
  icd_terms <- c(
    sprintf("\"ICD-10 %s\"[All Fields]", icd3),
    sprintf("\"ICD10 %s\"[All Fields]" , icd3),
    sprintf("\"ICD 10 %s\"[All Fields]", icd3),
    sprintf("(%s[All Fields] AND (ICD-10[All Fields] OR ICD10[All Fields]))", icd3)
  )
  
  bare_terms <- character()
  if (isTRUE(include_bare_code)) {
    bare_terms <- c(sprintf("\"%s\"[All Fields]", icd3))
  }
  
  mesh_clause <- if (!is.null(mesh_term) && nzchar(mesh_term)) {
    sprintf("(\"%s\"[MeSH Terms])", mesh_term)
  } else NULL
  
  all_terms <- c(terms_name, icd_terms, bare_terms, mesh_clause)
  paste0("(", paste(all_terms, collapse = " OR "), ")")
}

pubmed_count <- function(term, email = NULL, retries = 3L, pause = 0.35) {
  if (!is.null(email) && nzchar(email)) options(entrez.email = email)
  key <- Sys.getenv("RENTREZ_KEY", "")
  if (nzchar(key) && exists("set_entrez_key", mode = "function")) {
    try(set_entrez_key(key), silent = TRUE)
  }
  last_err <- NULL
  for (i in seq_len(retries)) {
    Sys.sleep(pause)
    out <- try(rentrez::entrez_search(db = "pubmed", term = term, retmax = 0L),
               silent = TRUE)
    if (!inherits(out, "try-error")) return(out$count %||% NA_integer_)
    last_err <- out
    Sys.sleep(pause * i)
  }
  warning(sprintf("PubMed query failed after %d retries: %s", retries, last_err))
  NA_integer_
}

# ---- MAIN: vectors in, tibble out ----
# icd_vec : character vector of ICD-10 3-char codes (e.g., "E11", "K20", ...)
# name_vec: character vector of corresponding diagnosis names (same length)
count_icd_pubmed_from_vectors <- function(icd_vec, name_vec,
                                          include_bare_code = FALSE,
                                          email = NULL,
                                          out_csv = NULL) {
  stopifnot(length(icd_vec) == length(name_vec))
  df <- tibble(icd3 = toupper(icd_vec), name = name_vec)
  
  res <- df |>
    mutate(
      mesh_term = map_chr(icd3, ~ guess_mesh(.x, name)),
      query     = pmap_chr(list(icd3, name, include_bare_code, mesh_term),
                           ~ make_query(..1, ..2, ..3, ..4)),
      pubmed_n  = map_int(query, ~ pubmed_count(.x, email = email))
    ) |>
    select(icd3, name, pubmed_n, query)
  
  if (!is.null(out_csv)) readr::write_csv(res, out_csv)
  res
}

# ---- EXAMPLE USAGE ----
icd_vec  <- c("E11","K20")
name_vec <- c("Type 2 diabetes mellitus","Oesophagitis")
result <- count_icd_pubmed_from_vectors(icd_vec, name_vec, include_bare_code = FALSE)
View(result)






result <- count_icd_pubmed_from_vectors(icd$icd3, icd$name, include_bare_code = FALSE)





str(result)


names(prevalence)[4] <- "prevalence"
prevalence$n <- NULL
prevalence$n_all <- NULL

names(hospitalisations)[4] <- "hospitalisation"
hospitalisations$n <- NULL

gene$X <- NULL
gene$Connectivity.z.score <- NULL
names(gene) <- c("icd_code", "NumberOfGenes")


result$query <- NULL
result$synonyms  <- NULL
names(result)[1] <- "icd_code"

df_final <- prevalence
df_final <- left_join(df_final,hospitalisations )
df_final <- left_join(df_final,gene )
df_final <- left_join(df_final, result)

str(df_final)
summary(df_final)


saveRDS(df_final, "ICDs_Parameters.rds")

write.csv(df_final, "ICDs_Parameters.csv", row.names = F)

