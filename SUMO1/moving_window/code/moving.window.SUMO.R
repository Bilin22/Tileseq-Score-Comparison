library(ggplot2)
library(dplyr)
library(stringr)
library(ggpubr)
library(tidyr)
library(patchwork)
library(ggrepel)

# old data
map.2019 <- read.csv(file = "../Tileseq_Scores/SUMO1/moving_window//data/urn_mavedb_00000001-b-2_scores.csv") %>% 
  select(hgvs_pro, score, se) %>% 
  distinct() %>% 
  mutate(mut=substr(hgvs_pro, nchar(hgvs_pro)-2, nchar(hgvs_pro))) %>%
  mutate(wt=substr(hgvs_pro, 3, 5)) %>% 
  mutate(type= ifelse(grepl("=", mut), "synonymous", ifelse(grepl("Ter$", mut), 
                                                            "nonsense", "missense"))) %>% 
  mutate(position = ifelse(type == "synonymous", substr(hgvs_pro, 6, nchar(hgvs_pro) - 1), 
                           substr(hgvs_pro, 6, nchar(hgvs_pro) - 3)))

# new data, add a new column for amino acid position
map.2023 <- read.csv(file = "../Tileseq_Scores/SUMO1/moving_window/data/select_t1_simple_aa.csv", 
                     skip = 16) %>% 
  mutate(mut=substr(hgvs_pro, nchar(hgvs_pro)-2, nchar(hgvs_pro))) %>%
  mutate(wt=substr(hgvs_pro, 3, 5)) %>% 
  mutate(type= ifelse(grepl("=", mut), "synonymous", ifelse(grepl("Ter$", mut), 
                                                            "nonsense", "missense"))) %>% 
  mutate(position = ifelse(type == "synonymous", substr(hgvs_pro, 6, nchar(hgvs_pro) - 1), 
                                                        substr(hgvs_pro, 6, nchar(hgvs_pro) - 3)))

# try plotting the score vs. amino acid position
ggplot(map.2023, aes(x = position, y = score, color = type)) +
  geom_line(aes(group = type)) +
  scale_color_manual(values = c("grey59", "red2", "chartreuse3")) +
  facet_grid(~ type)+
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 60, hjust = 1)
  )


ggplot(map.2023, aes(x = position, y = score, color = type)) +
  # geom_point(aes(color = type, alpha = 0.8)) +
  geom_line(aes(group = type)) +
  scale_color_manual(values = c("grey59", "red2", "chartreuse3")) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 70))

# build a matched table
match_df <- merge(map.2023, map.2019, by = "hgvs_pro") %>% 
  rename(score.2023 = score.x, se.2023 = se.x, score.2019 = score.y, se.2019 = se.y, 
         mut = mut.x, wt = wt.x, type = type.x) %>% 
  select(hgvs_pro, score.2023, se.2023, score.2019, se.2019, wt, mut, type) %>% 
  mutate(plotname = ifelse(abs(score.2023 - score.2019) >= 0.4, hgvs_pro, "")) %>% 
  mutate(position = ifelse(type == "synonymous", substr(hgvs_pro, 6, nchar(hgvs_pro) - 1), 
                           substr(hgvs_pro, 6, nchar(hgvs_pro) - 3))) %>% 
  mutate(difference = abs(score.2023 - score.2019))

ggplot(match_df, aes(x = position, y = difference, color = type)) +
  # geom_point(aes(color = type, alpha = 0.8)) +
  geom_line(aes(group = type)) +
  scale_color_manual(values = c("grey59", "red2", "chartreuse3")) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 70))

