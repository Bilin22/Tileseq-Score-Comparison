library(ggplot2)
library(dplyr)
library(tidyr)

# old score
map.2019 <- read.csv(file = "../Tileseq_Scores/SUMO1/moving_window/data/urn_mavedb_00000001-b-2_scores.csv") %>% 
  select(hgvs_pro, score, se) %>% 
  distinct() %>% 
  mutate(mut=substr(hgvs_pro, nchar(hgvs_pro)-2, nchar(hgvs_pro))) %>%
  mutate(wt=substr(hgvs_pro, 3, 5)) %>% 
  mutate(type= ifelse(grepl("=", mut), "synonymous", ifelse(grepl("Ter$", mut), 
                                                            "nonsense", "missense"))) %>% 
  mutate(position = ifelse(type == "synonymous", substr(hgvs_pro, 6, nchar(hgvs_pro) - 1), 
                           substr(hgvs_pro, 6, nchar(hgvs_pro) - 3)))

# new score
map.2023 <- read.csv(file = "../Tileseq_Scores/SUMO1/moving_window/data/select_t1_simple_aa.csv", 
                     skip = 16) %>% 
  mutate(mut=substr(hgvs_pro, nchar(hgvs_pro)-2, nchar(hgvs_pro))) %>%
  mutate(wt=substr(hgvs_pro, 3, 5)) %>% 
  mutate(type= ifelse(grepl("=", mut), "synonymous", ifelse(grepl("Ter$", mut), 
                                                            "nonsense", "missense"))) %>% 
  mutate(position = ifelse(type == "synonymous", substr(hgvs_pro, 6, nchar(hgvs_pro) - 1), 
                                                        substr(hgvs_pro, 6, nchar(hgvs_pro) - 3)))

# VARITY_R Score & VARITY_ER Score

# a ref map for amino acid
aa <- c("A" = "Ala", "R" = "Arg", "N" = "Asn", "D" = "Asp", "C" = "Cys", 
           "E" = "Glu", "Q" = "Gln", "G" = "Gly", "H" = "His", "I" = "Ile",
           "L" = "Leu", "K" = "Lys", "M" = "Met", "F" = "Phe", "P" = "Pro",
           "S" = "Ser", "T" = "Thr", "W" = "Trp", "Y" = "Tyr", "V" = "Val")


VARITY_R.score <- read.csv(file = "../Tileseq_Scores/SUMO1/moving_window/data/SUMO1[P63165]_1-101_VARITY_R_20230721193531406123.csv") %>% 
  select(aa_pos, aa_ref, aa_alt, VARITY_R) %>% 
  drop_na() %>% 
  mutate(hgvs_pro = paste("p.", aa[aa_ref], aa_pos, aa[aa_alt], sep = '')) %>% 
  distinct() %>% 
  select(hgvs_pro, aa_pos, VARITY_R)

VARITY_ER.score <- read.csv(file = "../Tileseq_Scores/SUMO1/moving_window/data/SUMO1[P63165]_1-101_VARITY_ER_20230721193815096805.csv") %>% 
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
match_df <- merge(map.2023, map.2019, by = "hgvs_pro") %>% 
  rename(score.2023 = score.x, se.2023 = se.x, score.2019 = score.y, se.2019 = se.y, 
         mut = mut.x, wt = wt.x, type = type.x) %>% 
  select(hgvs_pro, score.2023, se.2023, score.2019, se.2019, wt, mut, type) %>% 
  mutate(plotname = ifelse(abs(score.2023 - score.2019) >= 0.4, hgvs_pro, "")) %>% 
  mutate(position = ifelse(type == "synonymous", substr(hgvs_pro, 6, nchar(hgvs_pro) - 1), 
                           substr(hgvs_pro, 6, nchar(hgvs_pro) - 3)))


# match df between 2 map scores and VARITY





# select 3 columns
df <- match_df %>% 
  select(position, score.2023, score.2019) %>% 
  mutate(position = as.numeric(position))

ws <- 30

new.columns <- t(sapply(df$position, function(v) {
  relevant <- which(abs(df$position - v) < ws / 2)
  score_2023 <- df$score.2023[relevant]
  score_2019 <- df$score.2019[relevant]
  rho <- cor(df[relevant,c("score.2023","score.2019")],method="spearman")[1,2]
  return(c(position = v, score2023.mean = mean(score_2023), score2019.mean = mean(score_2019),
           rho = rho))
}))

# data frame for plotting
window.df <- as.data.frame(new.columns)

# moving window for fitness score
# set the color manually
colors <- c("score2023" = "blue3", "score2019" = "red3")
score.window <- ggplot(window.df, aes(x = position)) +
  geom_line(aes(y = score2023.mean, color = "score2023")) +
  geom_line(aes(y = score2019.mean, color = "score2019")) +
  labs(title = "Average Fitness Scores of SUMO1", x = "position", y = "mean fitness score", color = "Legend") +
  scale_color_manual(values = colors) +
  theme_bw() +
  theme(legend.position = c(0.9, 0.9), 
        legend.background = element_rect(fill = "white", colour = "grey22"))
score.window + guides(color = guide_legend(reverse = TRUE))
# ggsave(filename = "../Tileseq_Scores/SUMO1/moving_window/output/scores.png", dpi = 700,
#        height = 6, width = 8)


# moving window for spearman's rho
ggplot(window.df, aes(x = position, y = rho)) +
  geom_line() +
  labs(title = "Spearman Correlation Coefficient between 2023 & 2019 Fitness Score of SUMO1",
       x = "position", y = "Spearman's rho") +
  theme_bw()
# ggsave(filename = "../Tileseq_Scores/SUMO1/moving_window/output/spearman.png", dpi = 700,
#        height = 6, width = 8)


