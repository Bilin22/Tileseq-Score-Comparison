library(ggplot2)
library(dplyr)
library(tidyr)

# old data
map_2021 <- read.csv(file = "../Tileseq_Scores/GDI1/moving_window/data/urn_mavedb_00000066-a-1_scores.csv") %>% 
  select(hgvs_pro, score, score_sd) %>% 
  mutate(mut=substr(hgvs_pro, nchar(hgvs_pro)-2, nchar(hgvs_pro))) %>%
  mutate(wt=substr(hgvs_pro, 3, 5)) %>% 
  mutate(type= ifelse(grepl("=", mut), "synonymous", ifelse(grepl("Ter$", mut), 
                                                            "nonsense", "missense")))

# new data, make it match the format with 2021 ver on MAVEdb
map_2023 <- read.csv(file = "../Tileseq_Scores/GDI1/Score/data/select_t1_simple_aa_floored.csv", 
                     skip = 16) %>% 
  mutate(hgvs_pro = ifelse(grepl("=$", hgvs_pro), 
                           paste(substr(hgvs_pro, 1, 5), 
                           substr(hgvs_pro, 6, nchar(hgvs_pro) - 1), 
                           substr(hgvs_pro,3, 5), sep = ''), hgvs_pro)) %>% 
  mutate(mut=substr(hgvs_pro, nchar(hgvs_pro)-2, nchar(hgvs_pro))) %>%
  mutate(wt=substr(hgvs_pro, 3, 5)) %>% 
  mutate(type= ifelse(grepl("=", mut), "synonymous", ifelse(grepl("Ter$", mut), 
                                                            "nonsense", "missense")))
# VARITY_R Score & VARITY_ER Score

# a ref map for amino acid
aa <- c("A" = "Ala", "R" = "Arg", "N" = "Asn", "D" = "Asp", "C" = "Cys", 
        "E" = "Glu", "Q" = "Gln", "G" = "Gly", "H" = "His", "I" = "Ile",
        "L" = "Leu", "K" = "Lys", "M" = "Met", "F" = "Phe", "P" = "Pro",
        "S" = "Ser", "T" = "Thr", "W" = "Trp", "Y" = "Tyr", "V" = "Val")


VARITY_R.score <- read.csv(file = "../Tileseq_Scores/GDI1/moving_window/data/GDI1[P31150]_1-447_VARITY_R_20230727192515529432.csv") %>% 
  select(aa_pos, aa_ref, aa_alt, VARITY_R) %>% 
  drop_na() %>% 
  mutate(hgvs_pro = paste("p.", aa[aa_ref], aa_pos, aa[aa_alt], sep = '')) %>% 
  distinct() %>% 
  select(hgvs_pro, aa_pos, VARITY_R)

VARITY_ER.score <- read.csv(file = "../Tileseq_Scores/GDI1/moving_window/data/GDI1[P31150]_1-447_VARITY_ER_20230727192658498621.csv") %>% 
  select(aa_pos, aa_ref, aa_alt, VARITY_ER) %>% 
  drop_na() %>% 
  mutate(hgvs_pro = paste("p.", aa[aa_ref], aa_pos, aa[aa_alt], sep = '')) %>% 
  distinct() %>% 
  select(hgvs_pro, aa_pos, VARITY_ER)


# a match df for VARITY scores
VARITY_df <- merge(VARITY_R.score, VARITY_ER.score, by = "hgvs_pro") %>% 
  select(hgvs_pro, aa_pos.x, VARITY_R, VARITY_ER) %>% 
  rename(position = aa_pos.x)


# build a matched table
match_df <- merge(map_2023, map_2021, by = "hgvs_pro") %>% 
  rename(score.2023 = score.x, score.2021 = score.y) %>% 
  select(hgvs_pro, score.2023, score.2021, type)
  # mutate(plotname = ifelse(abs(score.2023 - score.2021) >= 0.4, hgvs_pro, "")) %>% 
  # mutate(position = ifelse(type == "synonymous", substr(hgvs_pro, 6, nchar(hgvs_pro) - 1), 
  #                          substr(hgvs_pro, 6, nchar(hgvs_pro) - 3)))


# match df between 2 map scores and VARITY
match_VARITY_df <- merge(match_df, VARITY_df, by = "hgvs_pro")


# select 3 columns
df <- match_VARITY_df %>% 
  select(position, score.2023, score.2021, VARITY_R, VARITY_ER)

# prepare the columns for moving window analysis
ws <- 30

new.columns <- t(sapply(df$position, function(v) {
  relevant <- which(abs(df$position - v) < ws / 2)
  score_2023 <- df$score.2023[relevant]
  score_2021 <- df$score.2021[relevant]
  VARITY_R.score <- df$VARITY_R[relevant]
  VARITY_ER.score <- df$VARITY_ER[relevant]
  rho <- cor(df[relevant,c("score.2023","score.2021")],method="spearman")[1,2]
  rho_VARITY_R_2023 <- cor(df[relevant,c("score.2023","VARITY_R")],method="spearman")[1,2]
  rho_VARITY_ER_2023 <- cor(df[relevant,c("score.2023","VARITY_ER")],method="spearman")[1,2]
  rho_VARITY_R_2021 <- cor(df[relevant,c("score.2021","VARITY_R")],method="spearman")[1,2]
  rho_VARITY_ER_2021 <- cor(df[relevant,c("score.2021","VARITY_ER")],method="spearman")[1,2]
  return(c(position = v, score2023.mean = mean(score_2023), score2021.mean = mean(score_2021), 
           VARITY_R.mean = mean(VARITY_R.score), VARITY_ER.mean = mean(VARITY_ER.score), 
           rho_VARITY_R_2023 = rho_VARITY_R_2023, rho_VARITY_ER_2023 = rho_VARITY_ER_2023, 
           rho_VARITY_R_2021 = rho_VARITY_R_2021, rho_VARITY_ER_2021 = rho_VARITY_ER_2021,
           rho = rho))
}))

# data frame for plotting
window.df <- as.data.frame(new.columns)

V_colors <- c("VARITY_R" = "magenta3", "VARITY_ER" = "slateblue3")
# moving window for VARITY and score 2023
ggplot(window.df, aes(x = position)) +
  geom_line(aes(y = rho_VARITY_R_2023, color = "VARITY_R")) +
  geom_line(aes(y = rho_VARITY_ER_2023, color = "VARITY_ER")) +
  scale_color_manual(values = V_colors) +
  theme_bw() +
  labs(title = "Correlation between 2023 Fitness Scores & VARITY Scores of GDI1",
       x = "position", y = "Spearman's rho", colour = "Scores") +
  theme(legend.position = c(0.9, 0.9), 
        legend.background = element_rect(fill = "white", colour = "grey22"))
# ggsave(filename = "../Tileseq_Scores/GDI1/moving_window/output/spearman_VARITY_2023.png", dpi = 700,
#        height = 6, width = 8)


# moving window for VARITY and score 2021
ggplot(window.df, aes(x = position)) +
  geom_line(aes(y = rho_VARITY_R_2021, color = "VARITY_R")) +
  geom_line(aes(y = rho_VARITY_ER_2021, color = "VARITY_ER")) +
  scale_color_manual(values = V_colors) +
  theme_bw() +
  labs(title = "Correlation between 2021 Fitness Scores & VARITY Scores of GDI1",
       x = "position", y = "Spearman's rho", colour = "Scores") +
  theme(legend.position = c(0.9, 0.9), 
        legend.background = element_rect(fill = "white", colour = "grey22"))
# ggsave(filename = "../Tileseq_Scores/GDI1/moving_window/output/spearman_VARITY_2021.png", dpi = 700,
#        height = 6, width = 8)


# moving window for mean fitness score 2023 and VARITY scores
colors_4 <- c("VARITY_R" = "magenta3", "VARITY_ER" = "slateblue3", 
              "Score_2023" = "darkgreen", "Score_2021" = "orange")

comparison <- ggplot(window.df, aes(x = position)) +
  geom_line(aes(y = score2023.mean , color = "Score_2023")) +
  geom_line(aes(y = score2021.mean , color = "Score_2021")) +
  geom_line(aes(y = VARITY_R.mean, color = "VARITY_R")) +
  geom_line(aes(y = VARITY_ER.mean, color = "VARITY_ER")) +
  scale_color_manual(values = colors_4) +
  theme_bw() +
  labs(title = "Average Fitness Scores & VARITY Scores of GDI1",
       x = "position", y = "Mean Fitness Score", colour = "Scores") 
comparison + guides(color = guide_legend(reverse = TRUE))
# ggsave(filename = "../Tileseq_Scores/GDI1/moving_window/output/scores_VARITY_year.png", dpi = 700,
#        height = 6, width = 8)



# 2. moving window for fitness score
# set the color manually
colors <- c("score2023" = "darkgreen", "score2021" = "orange")
score.window <- ggplot(window.df, aes(x = position)) +
  geom_line(aes(y = score2023.mean, color = "score2023")) +
  geom_line(aes(y = score2021.mean, color = "score2021")) +
  labs(title = "Average Fitness Scores of GDI1", x = "position", y = "mean fitness score", color = "Legend") +
  scale_color_manual(values = colors) +
  theme_bw() +
  theme(legend.position = c(0.9, 0.9), 
        legend.background = element_rect(fill = "white", colour = "grey22"))
score.window + guides(color = guide_legend(reverse = TRUE))
# ggsave(filename = "../Tileseq_Scores/GDI1/moving_window/output/scores.png", dpi = 700,
#        height = 6, width = 8)


# moving window for spearman's rho
ggplot(window.df, aes(x = position, y = rho)) +
  geom_line() +
  labs(title = "Spearman Correlation Coefficient between 2023 & 2021 Fitness Score of GDI1",
       x = "position", y = "Spearman's rho") +
  theme_bw()
# ggsave(filename = "../Tileseq_Scores/GDI1/moving_window/output/spearman.png", dpi = 700,
#        height = 6, width = 8)