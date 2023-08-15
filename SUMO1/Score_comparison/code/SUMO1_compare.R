library(ggplot2)
library(dplyr)
library(stringr)
library(ggpubr)
library(tidyr)
library(patchwork)
library(ggrepel)

# old data
map_2019 <- read.csv(file = "../Tileseq_Scores/SUMO1/Score_comparison/data/urn_mavedb_00000001-b-2_scores.csv") %>% 
  select(hgvs_pro, score, se) %>% 
  distinct() %>% 
  mutate(mut=substr(hgvs_pro, nchar(hgvs_pro)-2, nchar(hgvs_pro))) %>%
  mutate(wt=substr(hgvs_pro, 3, 5)) %>% 
  mutate(type= ifelse(grepl("=", mut), "synonymous", ifelse(grepl("Ter$", mut), 
                                                            "nonsense", "missense")))

# new data
map_2023 <- read.csv(file = "../Tileseq_Scores/SUMO1/Score_comparison/data/select_t1_simple_aa.csv", 
                     skip = 16) %>% 
  mutate(mut=substr(hgvs_pro, nchar(hgvs_pro)-2, nchar(hgvs_pro))) %>%
  mutate(wt=substr(hgvs_pro, 3, 5)) %>% 
  mutate(type= ifelse(grepl("=", mut), "synonymous", ifelse(grepl("Ter$", mut), 
                                                            "nonsense", "missense")))

# build a matched table
match_df <- merge(map_2023, map_2019, by = "hgvs_pro") %>% 
  rename(score.2023 = score.x, se.2023 = se.x, score.2019 = score.y, se.2019 = se.y, 
         mut = mut.x, wt = wt.x, type = type.x) %>% 
  select(hgvs_pro, score.2023, se.2023, score.2019, se.2019, wt, mut, type) %>% 
  mutate(plotname = ifelse(abs(score.2023 - score.2019) >= 0.4, hgvs_pro, "")) 
# we label some pts if two scores disagree with each other by more than 0.4

# boxplot for se
se.2023 <- ggplot(match_df, aes(x = type, y = se.2023)) + 
  geom_boxplot(varwidth = TRUE, alpha = 0.8, aes(color = type)) +
  theme_bw() +
  scale_color_manual(values = c("grey59", "red2", "chartreuse3")) +
  labs(x = "mutation type", y = "standard error (2023)")


se.2019 <- ggplot(match_df, aes(x = type, y = se.2019)) + 
  geom_boxplot(varwidth = TRUE, alpha = 0.8, aes(color = type)) +
  theme_bw() +
  scale_color_manual(values = c("grey59", "red2", "chartreuse3")) +
  labs(x = "mutation type", y = "standard error (2019)")


se.boxplot <- se.2019 + se.2023 + plot_layout(guides = "collect")
# ggsave(file = "../Tileseq_Scores/SUMO1/Score_comparison/output/se_boxplot.png",
#        dpi = 700, width = 8, height = 5)


# scatter plot for scores, sep by mutation type
combined_mutation <- ggplot(match_df, aes(x = score.2023, y = score.2019)) + 
  geom_point(aes(color=type),  alpha=0.8) +
  geom_text_repel(aes(label = plotname),size = 1.8,max.overlaps = 8) +
  facet_wrap(~ type) +
  scale_color_manual(values = c("grey59", "red2", "chartreuse3")) +
  stat_cor(method = "spearman", cor.coef.name = c("rho", "R")) +
  theme_bw()+
  theme(legend.position = "bottom") +
  theme(legend.key.size = unit(0.5, "cm"))
# ggsave(file = "../Tileseq_Scores/SUMO1/Score_comparison/output/comparison_mutation.png",
#        dpi = 700, width = 8, height = 6)


# combined scatter plot (sep by mutation types)
overall_mutation <- ggplot(match_df, aes(x = score.2023, y = score.2019)) + 
  geom_point(aes(color=type), alpha=0.8) +
  stat_cor(method = "spearman", cor.coef.name = c("rho", "R")) +
  theme_bw()+
  theme(legend.position = "bottom") +
  theme(legend.key.size = unit(0.5, "cm")) +
  scale_color_manual(values = c("grey59", "red2", "chartreuse3"))
# ggsave(file = "../Tileseq_Scores/SUMO1/Score_comparison/output/comparison_overall.png",
#        dpi = 700, width = 7, height = 6)

# sep by aa
aa_mutation <- ggplot(match_df, aes(x = score.2023, y = score.2019)) + 
  geom_point(aes(color=type, alpha = 0.9)) +
  geom_text_repel(aes(label = plotname),size = 1.8,max.overlaps = 7) +
  scale_color_manual(values = c("grey59", "red2", "chartreuse3")) +
  # geom_smooth(method = "lm") +
  stat_cor(method = "spearman", cor.coef.name = c("rho", "R")) +
  # geom_rug(aes(color=wt))+
  facet_wrap(~ wt) +
  theme_bw()+
  theme(legend.position = "bottom") +
  theme(legend.key.size = unit(0.5, "cm"))
# ggsave(file = "../Tileseq_Scores/SUMO1/Score_comparison/output/comparison_aa.png",
#        dpi = 700, width = 10, height = 8)


