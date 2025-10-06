# ---- Packages ----
library(dplyr)
library(tidyr)
library(ggplot2)
library(ggbump)    # install.packages("ggbump")
library(ggrepel)

# ---- Params ----
TOP_HILIGHT <- 30   # how many ICDs (by highest prevalence) to color
LABEL_N     <- 20   # how many to label on ends
OUT_PNG     <- "bump_icd2.png"
OUT_PDF     <- "bump_icd1.pdf"

# ---- Data prep ----
metrics_order <- c("prevalence", "hospitalisation", "NumberOfGenes", "pubmed_count")

names(df_final)[7] <- "pubmed_count"

df_clean <- df_final %>%
  mutate(
    prevalence      = coalesce(prevalence, 0),
    hospitalisation = coalesce(hospitalisation, 0),
    pubmed_count    = coalesce(pubmed_count, 0L),
    NumberOfGenes   = coalesce(NumberOfGenes, 0L)
  )

# top ICDs to highlight by prevalence
hilight_icd <- df_clean %>%
  arrange(desc(prevalence)) %>%
  slice_head(n = min(TOP_HILIGHT, nrow(.))) %>%
  pull(icd_code)

# long format + ranks
long <- df_clean %>%
  select(icd_code, all_of(metrics_order)) %>%
  pivot_longer(cols = all_of(metrics_order), names_to = "metric", values_to = "value") %>%
  mutate(metric = factor(metric, levels = metrics_order)) %>%
  group_by(metric) %>%
  mutate(rank = min_rank(desc(value))) %>%
  ungroup() %>%
  mutate(
    is_hilight = icd_code %in% hilight_icd,
    label_left  = ifelse(metric == first(metrics_order) & is_hilight, icd_code, NA),
    label_right = ifelse(metric == last(metrics_order)  & is_hilight, icd_code, NA)
  )

# limit labels to top LABEL_N
label_keep <- df_clean %>%
  arrange(desc(prevalence)) %>%
  slice_head(n = min(LABEL_N, nrow(.))) %>%
  pull(icd_code)

long <- long %>%
  mutate(
    label_left  = ifelse(icd_code %in% label_keep, label_left, NA),
    label_right = ifelse(icd_code %in% label_keep, label_right, NA)
  )

# ---- Plot ----
p <- ggplot(long, aes(x = metric, y = rank, group = icd_code)) +
  # grey background lines
  geom_bump(data = subset(long, !is_hilight), linewidth = 0.4, alpha = 0.25, colour = "grey60") +
  geom_point(data = subset(long, !is_hilight), size = 0.8, alpha = 0.25, colour = "grey60") +
  
  # highlighted lines
  geom_bump(data = subset(long, is_hilight), aes(color = icd_code), linewidth = 1.1, alpha = 0.9) +
  geom_point(data = subset(long, is_hilight), aes(color = icd_code), size = 1.6, alpha = 0.9) +
  
  # labels on left & right
  geom_text_repel(
    data = subset(long, !is.na(label_left)),
    aes(label = label_left),
    nudge_x = -0.25, direction = "y", hjust = 1,
    size = 3.2, min.segment.length = 0.05, segment.alpha = 0.4, seed = 1
  ) +
  geom_text_repel(
    data = subset(long, !is.na(label_right)),
    aes(label = label_right),
    nudge_x = 0.25, direction = "y", hjust = 0,
    size = 3.2, min.segment.length = 0.05, segment.alpha = 0.4, seed = 1
  ) +
  
  scale_y_reverse(expand = expansion(mult = c(0.03, 0.08))) +
  scale_x_discrete(expand = expansion(add = 0.4)) +
  scale_color_discrete(guide = "none") +
  labs(
    x = NULL, y = "Rank (1 = highest)",
    title = "ICD rank flow across metrics",
    subtitle = "prevalence → hospitalisation → NumberOfGenes → pubmed_count"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    panel.grid.minor = element_blank(),
    axis.text.x = element_text(face = "bold")
  )

print(p)

# ---- Save ----
ggsave(OUT_PNG, p, width = 12, height = 7, dpi = 300)
ggsave(OUT_PDF, p, width = 12, height = 7)




library(dplyr)
library(tidyr)
library(ggplot2)
library(ggbump)
library(ggrepel)

# ---- Params ----
TOP_VISIBLE   <- 30   # top ICDs (fully colored)
LABEL_N       <- 20   # labels on ends
metrics_order <- c("prevalence","hospitalisation","NumberOfGenes","pubmed_count")

# ---- A–Z color palette ----
letter_cols <- c(
  "#1AF239","#58F21A","#961D1A","#B41AF2","#FFCA01","#581AF2",
  "#1AF295","#1A95F2","#F2761A","#1A39F2","#F21A1A","#F2D31A",
  "#B41AF2","#1AF2F2","#FF6F61","#6B5B95","#88B04B","#F7CAA0",
  "#92A8D1","#955251","#B5A642","#009B77","#DD41AC","#45B8AC",
  "#EFF05A","#9B2323"
)
color_by_letter <- setNames(letter_cols, LETTERS)

# ---- Data prep ----
df_clean <- df %>%
  mutate(
    prevalence      = coalesce(prevalence, 0),
    hospitalisation = coalesce(hospitalisation, 0),
    pubmed_count    = coalesce(pubmed_count, 0L),
    NumberOfGenes   = coalesce(NumberOfGenes, 0L),
    first_letter    = substr(icd_code, 1, 1)
  )

# Top ICDs to highlight
top_icd <- df_clean %>%
  arrange(desc(prevalence)) %>%
  slice_head(n = min(TOP_VISIBLE, nrow(.))) %>%
  pull(icd_code)

# Long format with ranks
long <- df_clean %>%
  select(icd_code, first_letter, all_of(metrics_order)) %>%
  pivot_longer(cols = all_of(metrics_order), names_to = "metric", values_to = "value") %>%
  mutate(metric = factor(metric, levels = metrics_order)) %>%
  group_by(metric) %>%
  mutate(rank = min_rank(desc(value))) %>%
  ungroup() %>%
  mutate(
    is_top     = icd_code %in% top_icd,
    alpha_val  = ifelse(is_top, 1, 0.5),
    label_left  = ifelse(metric == first(metrics_order) & is_top, icd_code, NA),
    label_right = ifelse(metric == last(metrics_order)  & is_top, icd_code, NA)
  )

# Limit labels to top LABEL_N
label_keep <- df_clean %>%
  arrange(desc(prevalence)) %>%
  slice_head(n = min(LABEL_N, nrow(.))) %>%
  pull(icd_code)

long <- long %>%
  mutate(
    label_left  = ifelse(icd_code %in% label_keep, label_left, NA),
    label_right = ifelse(icd_code %in% label_keep, label_right, NA)
  )

# ---- Plot ----
p <- ggplot(long, aes(x = metric, y = rank, group = icd_code)) +
  # Non-top ICDs: grey links
  geom_bump(data = subset(long, !is_top),
            colour = "grey70", linewidth = 0.5, alpha = 0.5) +
  # Top ICDs: colored links
  geom_bump(data = subset(long, is_top),
            aes(color = first_letter), linewidth = 1.1, alpha = 1) +
  
  # Points: all colored by first letter (opacity depends on top/not)
  geom_point(aes(color = first_letter, alpha = alpha_val), size = 1.5) +
  
  # Labels on ends
  geom_text_repel(
    data = subset(long, !is.na(label_left)),
    aes(label = label_left), nudge_x = -0.25, direction = "y", hjust = 1,
    size = 3.1, min.segment.length = 0.05, segment.alpha = 0.4, seed = 1
  ) +
  geom_text_repel(
    data = subset(long, !is.na(label_right)),
    aes(label = label_right), nudge_x = 0.25, direction = "y", hjust = 0,
    size = 3.1, min.segment.length = 0.05, segment.alpha = 0.4, seed = 1
  ) +
  
  scale_y_reverse(expand = expansion(mult = c(0.03, 0.08))) +
  scale_x_discrete(expand = expansion(add = 0.4)) +
  scale_color_manual(values = color_by_letter, drop = FALSE, guide = "legend") +
  scale_alpha_identity() +
  labs(
    x = NULL, y = "Rank (1 = highest)",
     ) +
  theme_minimal(base_size = 12) +
  theme(
    panel.grid.minor = element_blank(),
    axis.text.x = element_text(face = "bold")
  )

print(p)


metrics_order    <- c("prevalence","hospitalisation","NumberOfGenes","pubmed_count")
display_levels   <- c("Prevalence","Hospitalisation Rate","Number of Genes","PubMed Count")
metric_labels    <- setNames(display_levels, metrics_order)

# keep these from before
# top_icd, label_keep already computed

long <- df_clean %>%
  select(icd_code, first_letter, all_of(metrics_order)) %>%
  tidyr::pivot_longer(cols = tidyselect::all_of(metrics_order),
                      names_to = "metric", values_to = "value") %>%
  mutate(
    metric = factor(metric, levels = metrics_order, labels = display_levels)
  ) %>%
  group_by(metric) %>%
  mutate(rank = dplyr::min_rank(dplyr::desc(value))) %>%
  ungroup() %>%
  mutate(
    is_top   = icd_code %in% top_icd,
    alpha_bg = ifelse(is_top, 1, 0.5),
    
    # label only ends, only for top-N you chose to label
    label_left  = ifelse(metric == display_levels[1] & icd_code %in% label_keep, icd_code, NA),
    label_right = ifelse(metric == display_levels[length(display_levels)] & icd_code %in% label_keep, icd_code, NA)
  )



# ---- Plot ----
p <- ggplot(long, aes(x = metric, y = rank, group = icd_code)) +
  # Non-top ICDs: grey links
  geom_bump(data = subset(long, !is_top),
            colour = "grey70", linewidth = 0.5, alpha = 0.5) +
  
  # Top ICDs: colored links
  geom_bump(data = subset(long, is_top),
            aes(color = first_letter), linewidth = 1.1, alpha = 1) +
  
  # Points: all colored by first letter (alpha only for drawing, not legend)
  geom_point(aes(color = first_letter),
             size = 1.5, alpha = 0.9, show.legend = TRUE) +
  
  # Labels on ends
  geom_text_repel(
    data = subset(long, !is.na(label_left)),
    aes(label = label_left), nudge_x = -0.25, direction = "y", hjust = 1,
    size = 3.1, min.segment.length = 0.05, segment.alpha = 0.4, seed = 1
  ) +
  geom_text_repel(
    data = subset(long, !is.na(label_right)),
    aes(label = label_right), nudge_x = 0.25, direction = "y", hjust = 0,
    size = 3.1, min.segment.length = 0.05, segment.alpha = 0.4, seed = 1
  ) +
  
  scale_y_reverse(expand = expansion(mult = c(0.03, 0.08))) +
  scale_x_discrete(expand = expansion(add = 0.4)) +
  scale_color_manual(values = color_by_letter, drop = FALSE, guide = guide_legend(override.aes = list(size = 3, alpha = 1))) +
  labs(x = NULL, y = "Rank (1 = highest)", color = NULL) +   # remove legend title
  theme_minimal(base_size = 12) +
  theme(
    panel.grid.minor = element_blank(),
    axis.text.x = element_text(face = "bold"),
    legend.title = element_blank()  # extra safeguard to drop title
  ) 

print(p)

# ---- Save ----
ggsave(OUT_PNG, p, width = 12, height = 7, dpi = 300)
ggsave(OUT_PDF, p, width = 12, height = 7)


long$metric <- factor(long$metric,
                      levels =  c("Prevalence","Hospitalisation Rate","PubMed Count","Number of Genes"))

long <- long %>%
  mutate(
    # keep labels on the first column (Prevalence) for top ICDs
    label_left  = ifelse(metric == "Prevalence" & icd_code %in% label_keep, icd_code, NA),
    # labels on PubMed Count (last column)
    label_right = ifelse(metric == "PubMed Count" & icd_code %in% label_keep, icd_code, NA),
    # NEW: labels also on Number of Genes
    label_genes = ifelse(metric == "Number of Genes" & icd_code %in% label_keep, icd_code, NA)
  )
# ---- Plot ----
p <- ggplot(long, aes(x = metric, y = rank, group = icd_code)) +
  # Non-top ICDs: grey links
  geom_bump(data = subset(long, !is_top),
            colour = "grey70", linewidth = 0.5, alpha = 0.5) +
  
  # Top ICDs: colored links
  geom_bump(data = subset(long, is_top),
            aes(color = first_letter), linewidth = 1.1, alpha = 1) +
  
  # Points: all colored by first letter (alpha only for drawing, not legend)
  geom_point(aes(color = first_letter),
             size = 1.5, alpha = 0.9, show.legend = TRUE) +
  
  # Labels on ends
  geom_text_repel(
    data = subset(long, !is.na(label_left)),
    aes(label = label_left), nudge_x = -0.25, direction = "y", hjust = 1,
    size = 3.1, min.segment.length = 0.05, segment.alpha = 0.4, seed = 1
  ) +
  # Labels on Number of Genes
  geom_text_repel(
    data = subset(long, !is.na(label_genes)),
    aes(label = label_genes),
    nudge_x = 0.25, direction = "y", hjust = 0,
    size = 3.1, min.segment.length = 0.05, segment.alpha = 0.4, seed = 2
  ) +  scale_y_reverse(expand = expansion(mult = c(0.03, 0.08))) +
  scale_x_discrete(expand = expansion(add = 0.4)) +
  scale_color_manual(values = color_by_letter, drop = FALSE, guide = guide_legend(override.aes = list(size = 3, alpha = 1))) +
  labs(x = NULL, y = "Rank (1 = highest)", color = NULL) +   # remove legend title
  theme_minimal(base_size = 12) +
  theme(
    panel.grid.minor = element_blank(),
    axis.text.x = element_text(face = "bold"),
    legend.title = element_blank()  # extra safeguard to drop title
  ) + 
  scale_x_discrete(
    limits = c("Prevalence","Hospitalisation Rate","PubMed Count", "Number of Genes"),
    expand = expansion(add = 0.4)
  )


print(p)
OUT_PNG     <- "bump_icd3.png"
OUT_PDF     <- "bump_icd3.pdf"
# ---- Save ----
ggsave(OUT_PNG, p, width = 12, height = 7, dpi = 300)
ggsave(OUT_PDF, p, width = 12, height = 7)


