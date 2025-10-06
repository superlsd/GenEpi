
setwd("/Users/elmadervic/Documents/Github/GenEpi")

# # 1) Expand to long format
# df_long <- db_icd_list %>%
#   unnest(codes) %>%
#   mutate(member = 1)
# 
# # 2) Pivot to wide binary matrix
# df_wide <- df_long %>%
#   pivot_wider(
#     names_from = database,
#     values_from = member,
#     values_fill = list(member = 0)
#   )
# 
# # 3) Remove codes column for UpSetR
# df_wide_mat <- df_wide %>% select(-codes) %>% as.data.frame()
# 
# # 4) UpSetR plot with blue bars
# UpSetR::upset(
#   df_wide_mat,
#   sets = unique(db_icd_list$database),
#   keep.order = TRUE,
#   main.bar.color = BLUE,
#   matrix.color   = BLUE,
#   sets.bar.color = BLUE
# )
# 



library(jsonlite)
library(tidyverse)
library(UpSetR)
all_dicts <- fromJSON("all_disease_gene_dicts.json")


# Make tibble: one row per database, shorten names
db_icd_list <- tibble(
  database = sub("_.*", "", names(all_dicts)),  # keep only text before first "_"
  codes    = lapply(all_dicts, function(db_dict) names(db_dict))
)

db_icd_list
library(tidyverse)
library(UpSetR)

BLUE <- "#1f77b4"

# 1) Expand to long format
df_long <- db_icd_list %>%
  unnest(codes) %>%
  mutate(member = 1)

# 2) Pivot to wide binary matrix
df_wide <- df_long %>%
  pivot_wider(
    names_from = database,
    values_from = member,
    values_fill = list(member = 0)
  )

# 3) Remove codes column for UpSetR
df_wide_mat <- df_wide %>% select(-codes) %>% as.data.frame()

# 4) UpSetR plot with blue bars
UpSetR::upset(
  df_wide_mat,
  sets = unique(db_icd_list$database),
  keep.order = TRUE,
  main.bar.color = BLUE,
  matrix.color   = BLUE,
  sets.bar.color = BLUE
)


library(tidyverse)
library(UpSetR)
library(jsonlite)

# Load data
all_dicts <- fromJSON("all_disease_gene_dicts.json")

str(all_dicts )

# Make tibble: one row per database, shorten names
db_icd_list <- tibble(
  database = sub("_.*", "", names(all_dicts)),  # keep only text before first "_"
  codes    = lapply(all_dicts, function(db_dict) names(db_dict))
)

db_icd_list
library(tidyverse)
library(UpSetR)

BLUE <- "#1f77b4"



######################################



BLUE <- "#1f77b4"

# 1) Expand to long format & capitalize database names
df_long <- db_icd_list %>%
  unnest(codes) %>%
  mutate(
    database = str_to_title(database),  # first letter uppercase
    member   = 1
  )

# 2) Pivot to wide binary matrix
df_wide <- df_long %>%
  pivot_wider(
    names_from = database,
    values_from = member,
    values_fill = list(member = 0)
  )

# 3) Remove codes column for UpSetR
df_wide_mat <- df_wide %>% select(-codes) %>% as.data.frame()

# 4) Order sets by size (largest first)
set_order <- df_long %>%
  count(database, name = "size") %>%
  arrange(desc(size)) %>%
  pull(database)
# --- Set font to Helvetica Neue (macOS) ---
# quartzFonts(
#   helvNeue = quartzFont(c("Helvetica Neue", "Helvetica Neue", "Helvetica Neue", "Helvetica Neue"))
# )
# quartz(family = "helvNeue")
# 

# 5) UpSetR plot with ordering and blue theme
UpSetR::upset(
  df_wide_mat,
  sets           = set_order,
  keep.order     = TRUE,
  main.bar.color = BLUE,
  matrix.color   = BLUE,
  sets.bar.color = BLUE
)



# text.scale     = 1.5  # roughly font size 15 in UpSetR



library(tidyverse)
library(tidyverse)

str(df_long )
# A–Z palette
letter_cols <- c(
  "#1AF239","#58F21A","#961D1A","#B41AF2","#FFCA01","#581AF2","#1AF295","#1A95F2",
  "#F2761A","#1A39F2","#F21A1A","#F2D31A","#B41AF2","#1AF2F2","#FF6F61","#6B5B95",
  "#88B04B","#F7CAA0","#92A8D1","#955251","#B5A642","#009B77","#DD41AC","#45B8AC",
  "#EFF05A","#9B2323"
)
names(letter_cols) <- LETTERS

# counts by database × first letter
set_size_by_letter <- df_long %>%
  mutate(first_letter = toupper(substr(codes, 1, 1))) %>%
  filter(first_letter %in% LETTERS) %>%
  count(database, first_letter, name = "n") %>%
  mutate(database = factor(database, levels = set_order))



# Expect df_long with column `codes` like "K21", "D50", "H65", etc.

library(dplyr)
library(stringr)
library(forcats)

# Order (adjust if you have a preferred order)
set_order <- if (exists("set_order")) set_order else NULL

# ---- Map ICD-10 code -> Chapter (I–XXII) ----
icd_chapter <- function(code) {
  code <- toupper(code)
  L <- substr(code, 1, 1)
  N <- suppressWarnings(as.numeric(str_extract(code, "(?<=^[A-Z])\\d{2}")))
  if (is.na(N)) return(NA_character_)
  
  # Chapters per WHO ICD-10
  if (L %in% c("A","B")) return("I")                               # A00–B99
  if (L == "C" || (L == "D" && N <= 48)) return("II")              # C00–D48
  if (L == "D" && N >= 50) return("III")                           # D50–D89
  if (L == "E") return("IV")                                       # E00–E90
  if (L == "F") return("V")                                        # F00–F99
  if (L == "G") return("VI")                                       # G00–G99
  if (L == "H" && N <= 59) return("VII")                           # H00–H59
  if (L == "H" && N >= 60) return("VIII")                          # H60–H95
  if (L == "I") return("IX")                                       # I00–I99
  if (L == "J") return("X")                                        # J00–J99
  if (L == "K") return("XI")                                       # K00–K93
  if (L == "L") return("XII")                                      # L00–L99
  if (L == "M") return("XIII")                                     # M00–M99
  if (L == "N") return("XIV")                                      # N00–N99
  if (L == "O") return("XV")                                       # O00–O99
  if (L == "P") return("XVI")                                      # P00–P96
  if (L == "Q") return("XVII")                                     # Q00–Q99
  if (L == "R") return("XVIII")                                    # R00–R99
  if (L %in% c("S","T")) return("XIX")                             # S00–T98
  if (L %in% c("V","W","X","Y")) return("XX")                      # V01–Y98
  if (L == "Z") return("XXI")                                      # Z00–Z99
  if (L == "U") return("XXII")                                     # U00–U99
  return(NA_character_)
}

# ---- Your requested 14-color palette (I–XIV) ----
chapter_cols <- c(
  I   = "#4589FF", # infectious
  II  = "#00A3FF", # neoplasms
  III = "#00B5F1", # blood/immune
  IV  = "#00C2CD", # endocrine
  V   = "#00CAA5", # mental
  VI  = "#53CD85", # nervous
  VII = "#7CCD66", # eye
  VIII= "#A3CA47", # ear
  IX  = "#CAC42A", # circulatory
  X   = "#EBB71E", # respiratory
  XI  = "#F9A22C", # digestive
  XII = "#FF8C3F", # skin
  XIII= "#FF7852", # musculoskeletal
  XIV = "#FF6666", # genitourinary
  XV  = "#D65DB1", # pregnancy (new)
  XVI = "#845EC2", # perinatal   (new)
  XVII= "#0081CF", # congenital  (new)
  XVIII="#2C73D2", # symptoms    (new)
  XIX = "#008E9B", # injury      (new)
  XX  = "#008F7A", # external    (new)
  XXI = "#FFC75F", # health services (new)
  XXII= "#F9F871"  # special purposes (new)
)
library(dplyr)
library(stringr)
library(ggplot2)

# ---------------------------
# ICD-10 chapter mapping
# ---------------------------
icd_chapter <- function(code) {
  code <- toupper(code)
  L <- substr(code, 1, 1)
  N <- suppressWarnings(as.numeric(str_extract(code, "(?<=^[A-Z])\\d{2}")))
  if (is.na(N)) return(NA_character_)
  
  if (L %in% c("A","B")) return("I")                               # A00–B99
  if (L == "C" || (L == "D" && N <= 48)) return("II")              # C00–D48
  if (L == "D" && N >= 50) return("III")                           # D50–D89
  if (L == "E") return("IV")                                       # E00–E90
  if (L == "F") return("V")                                        # F00–F99
  if (L == "G") return("VI")                                       # G00–G99
  if (L == "H" && N <= 59) return("VII")                           # H00–H59
  if (L == "H" && N >= 60) return("VIII")                          # H60–H95
  if (L == "I") return("IX")                                       # I00–I99
  if (L == "J") return("X")                                        # J00–J99
  if (L == "K") return("XI")                                       # K00–K93
  if (L == "L") return("XII")                                      # L00–L99
  if (L == "M") return("XIII")                                     # M00–M99
  if (L == "N") return("XIV")                                      # N00–N99
  if (L == "O") return("XV")                                       # O00–O99
  if (L == "P") return("XVI")                                      # P00–P96
  if (L == "Q") return("XVII")                                     # Q00–Q99
  if (L == "R") return("XVIII")                                    # R00–R99
  if (L %in% c("S","T")) return("XIX")                             # S00–T98
  if (L %in% c("V","W","X","Y")) return("XX")                      # V01–Y98
  if (L == "Z") return("XXI")                                      # Z00–Z99
  if (L == "U") return("XXII")                                     # U00–U99
  return(NA_character_)
}

# ---------------------------
# Palette for 22 chapters
# ---------------------------
chapter_cols <- c(
  I   = "#4589FF",
  II  = "#00A3FF",
  III = "#00B5F1",
  IV  = "#00C2CD",
  V   = "#00CAA5",
  VI  = "#53CD85",
  VII = "#7CCD66",
  VIII= "#A3CA47",
  IX  = "#CAC42A",
  X   = "#EBB71E",
  XI  = "#F9A22C",
  XII = "#FF8C3F",
  XIII= "#FF7852",
  XIV = "#FF6666",
  XV  = "#D65DB1",
  XVI = "#845EC2",
  XVII= "#0081CF",
  XVIII="#2C73D2",
  XIX = "#008E9B",
  XX  = "#008F7A",
  XXI = "#FFC75F",
  XXII= "#F9F871"
)

# ---------------------------
# Aggregate to chapters
# ---------------------------
set_size_by_chapter <- df_long %>%
  mutate(
    icd_chapter = vapply(codes, icd_chapter, FUN.VALUE = character(1)),
    icd_chapter = factor(icd_chapter,
                         levels = names(chapter_cols))
  ) %>%
  filter(!is.na(icd_chapter)) %>%
  count(database, icd_chapter, name = "n") %>%
  mutate(database = if (exists("set_order")) factor(database, levels = set_order) else database)

# ---------------------------
# Plot mirrored horizontal bars (chapter colors)
# ---------------------------
p_sets <- ggplot(set_size_by_chapter,
                 aes(x = n, y = database, fill = icd_chapter)) +
  geom_col() +
  scale_fill_manual(values = chapter_cols, drop = FALSE, guide = "none") +
  scale_x_reverse(expand = expansion(mult = c(0, 0.05))) +
  labs(x = NULL, y = NULL) +
  theme_minimal(base_size = 16) +
  theme(
    panel.background = element_rect(fill = "white", colour = NA),
    plot.background  = element_rect(fill = "white", colour = NA),
    panel.grid.major = element_line(colour = "grey80", size = 1),
    panel.grid.minor = element_blank(),
    axis.text.y      = element_blank(),
    axis.title       = element_blank(),
    axis.text.x      = element_blank(),
    axis.ticks.x     = element_blank(),
    plot.margin      = margin(2, 2, 2, 2)
  )

p_sets

p_sets <- ggplot(set_size_by_chapter,
                 aes(x = n, y = database, fill = icd_chapter)) +
  geom_col(colour = "white", size = 0.2) +   # <- thin white border
  scale_fill_manual(values = chapter_cols, drop = FALSE, guide = "none") +
  scale_x_reverse(expand = expansion(mult = c(0, 0.05))) +
  labs(x = NULL, y = NULL) +
  theme_minimal(base_size = 16) +
  theme(
    panel.background = element_rect(fill = "white", colour = NA),
    plot.background  = element_rect(fill = "white", colour = NA),
    panel.grid.major = element_line(colour = "grey80", size = 1),
    panel.grid.minor = element_blank(),
    axis.text.y      = element_blank(),
    axis.title       = element_blank(),
    axis.text.x      = element_blank(),
    axis.ticks.x     = element_blank(),
    plot.margin      = margin(2, 2, 2, 2)
  )

p_sets



# save at width 4, height 1
ggsave("set_size_by_letter_mirrored.png", p_sets, width = 31, height = 8.5, dpi = 300)
ggsave("set_size_by_letter_mirrored.pdf", p_sets, width = 31, height = 8.5, device = cairo_pdf)





library(tidyverse)
library(UpSetR)
library(jsonlite)

BLUE <- "#155a8a"  # darker blue


# Load combined JSON
all_dicts <- fromJSON("all_disease_gene_dicts.json")

# Build per-database GENE lists (unique, drop NA/empty)
db_gene_list <- tibble(
  database = sub("_.*", "", names(all_dicts)),
  genes = lapply(all_dicts, function(db_dict) {
    g <- unique(unlist(db_dict, use.names = FALSE))
    g <- g[!is.na(g) & nzchar(g)]
    sort(g)
  })
)

# Long -> Wide for UpSetR
df_long_g <- db_gene_list %>%
  unnest(genes, keep_empty = FALSE) %>%
  mutate(
    database = str_to_title(database),
    member   = 1L
  )

# Order databases by gene count (desc)
set_order_g <- df_long_g %>%
  count(database, name = "size") %>%
  arrange(desc(size)) %>%
  pull(database)

# Binary matrix: rows = genes, cols = databases
df_wide_g <- df_long_g %>%
  pivot_wider(
    names_from  = database,
    values_from = member,
    values_fill = list(member = 0L)
  )

# UpSetR expects only the set columns (drop the gene id column)
df_wide_mat_g <- df_wide_g %>% select(-genes) %>% as.data.frame()

# OPTIONAL (macOS) — Helvetica Neue
# quartzFonts(helvNeue = quartzFont(rep("Helvetica Neue", 4)))
# quartz(family = "helvNeue")

# Plot
UpSetR::upset(
  df_wide_mat_g,
  sets           = set_order_g,
  keep.order     = TRUE,
  main.bar.color = BLUE,
  matrix.color   = BLUE,
  sets.bar.color = BLUE
)


