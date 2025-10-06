setwd("/Users/elmadervic/Documents/Github/GenEpi")





library(dplyr)
library(ggplot2)
library(tidyverse)
library(ggrepel)
# 1) Load your data (it's a CSV even if extension says .pickle)
df <- readr::read_csv(
  "output/jaccard_pairs_from_n_plets_elma_list_fdr_dict_male.csv",
  col_types = cols(key = col_integer(), value = col_double())
)

# 2) Ensure every unit around the circle is present (even if missing in data)
max_k <- max(df$key, na.rm = TRUE)
df_all <- tibble(key = 1:max_k) %>%
  left_join(df, by = "key")


# Row index for order around the circle
df$row <- seq_len(nrow(df))



legend_breaks <- c(0.05, 0.5, 0.75)
legend_labels <- c("0.05", "0.5", "0.75")


library(ggplot2)
library(ggrepel)

legend_breaks <- c(0.05, 0.5, 0.75)
legend_labels <- c("0.05", "0.5", "0.75")

p <- ggplot(df, aes(x = row, y = key, size = key, color = value)) +
  # main points
  geom_point(alpha = 0.8) +
  # green border for p < 0.05
  geom_point(
    data = subset(df, value < 0.05),
    shape = 21, stroke = 1.4, fill = NA, color = "#77DD77",
    aes(size = key)
  ) +
  # labels for key values
  geom_text_repel(
    aes(label = key),
    size = 3,
    color = "black",
    show.legend = FALSE
  ) +
  scale_size_continuous(range = c(3, 15), guide = "none") +
  scale_color_gradientn(
    colors = c("#1f77b4", "#ffeda0", "#e31a1c"),
    breaks = legend_breaks,
    labels = legend_labels,
    name = "p-value"
  ) +
  labs(x = NULL, y = "Trajectory size") +
  theme_minimal(base_size = 14) +
  theme(
    panel.grid.minor = element_blank(),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.title.x = element_blank(),
    panel.grid.major.x = element_blank(),
    legend.position = c(0.05, 0.95),
    legend.justification = c("left", "top"),
    legend.background = element_rect(fill = alpha("white", 0.7), color = NA)
  )


# save at width 4, height 1
ggsave("pvalue_male.png", p, width = 10, height = 8.3, dpi = 300)
ggsave("pvalue_male.pdf", p, width = 10, height = 8.3,)




# 1) Load your data (it's a CSV even if extension says .pickle)
df <- readr::read_csv(
  "output/jaccard_pairs_from_n_plets_elma_list_fdr_dict_female.csv",
  col_types = cols(key = col_integer(), value = col_double())
)

# 2) Ensure every unit around the circle is present (even if missing in data)
max_k <- max(df$key, na.rm = TRUE)
df_all <- tibble(key = 1:max_k) %>%
  left_join(df, by = "key")


# Row index for order around the circle
df$row <- seq_len(nrow(df))



legend_breaks <- c(0.05, 0.5, 0.75)
legend_labels <- c("0.05", "0.5", "0.75")


library(ggplot2)
library(ggrepel)

legend_breaks <- c(0.05, 0.5, 0.75)
legend_labels <- c("0.05", "0.5", "0.75")

p <- ggplot(df, aes(x = row, y = key, size = key, color = value)) +
  # main points
  geom_point(alpha = 0.8) +
  # green border for p < 0.05
  geom_point(
    data = subset(df, value < 0.05),
    shape = 21, stroke = 1.4, fill = NA, color = "#77DD77",
    aes(size = key)
  ) +
  # labels for key values
  geom_text_repel(
    aes(label = key),
    size = 3,
    color = "black",
    show.legend = FALSE
  ) +
  scale_size_continuous(range = c(3, 15), guide = "none") +
  scale_color_gradientn(
    colors = c("#1f77b4", "#ffeda0", "#e31a1c"),
    breaks = legend_breaks,
    labels = legend_labels,
    name = "p-value"
  ) +
  labs(x = NULL, y = "Trajectory size") +
  theme_minimal(base_size = 14) +
  theme(
    panel.grid.minor = element_blank(),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.title.x = element_blank(),
    panel.grid.major.x = element_blank(),
    legend.position = c(0.05, 0.95),
    legend.justification = c("left", "top"),
    legend.background = element_rect(fill = alpha("white", 0.7), color = NA)
  )


# save at width 4, height 1
ggsave("pvalue_female.png", p, width = 10, height = 8.3, dpi = 300)
ggsave("pvalue_female.pdf", p, width = 10, height = 8.3)




