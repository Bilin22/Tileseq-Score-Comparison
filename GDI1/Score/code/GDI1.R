library(ggplot2)
library(dplyr)
library(stringr)
library(ggpubr)
library(tidyr)
library(patchwork)
library(ggrepel)

# old data
map_2021 <- read.csv(file = "../Tileseq_Scores/GDI1/Score/data/score2021.csv") %>% 
  select(hgvsp, score, se) %>% 
  distinct() %>% 
  mutate(mut=substr(hgvsp, nchar(hgvsp)-2, nchar(hgvsp))) %>%
  mutate(wt=substr(hgvsp, 3, 5)) %>% 
  mutate(type= ifelse(grepl("=", mut), "synonymous", ifelse(grepl("Ter$", mut), 
                                                            "nonsense", "missense"))) %>% 
  rename("hgvs_pro" = hgvsp)

# new data
map_2023 <- read.csv(file = "../Tileseq_Scores/GDI1/Score/data/select_t1_simple_aa_floored.csv", 
                     skip = 16) %>% 
  mutate(mut=substr(hgvs_pro, nchar(hgvs_pro)-2, nchar(hgvs_pro))) %>%
  mutate(wt=substr(hgvs_pro, 3, 5)) %>% 
  mutate(type= ifelse(grepl("=", mut), "synonymous", ifelse(grepl("Ter$", mut), 
                                                            "nonsense", "missense")))

# build a matched table
match_df <- merge(map_2023, map_2021, by = "hgvs_pro") %>% 
  rename(score.2023 = score.x, se.2023 = se.x, score.2021 = score.y, se.2021 = se.y, 
         mut = mut.x, wt = wt.x, type = type.x) %>% 
  select(hgvs_pro, score.2023, se.2023, score.2021, se.2021, wt, mut, type) %>% 
  mutate(plotname = ifelse(abs(score.2023 - score.2021 >= 0.4), hgvs_pro, ""))

# boxplot for se
se.2023 <- ggplot(match_df, aes(x = type, y = se.2023)) + 
  geom_boxplot(varwidth = TRUE, alpha = 0.8) +
  theme_light() +
  labs(x = "mutation type", y = "standard error (2023)") +
  theme(
    panel.grid.major.x = element_blank(),
    panel.border = element_blank(),
    axis.ticks.x = element_blank()
  )

se.2021 <- ggplot(match_df, aes(x = type, y = se.2021)) + 
  geom_boxplot(varwidth = TRUE, alpha = 0.8) +
  theme_light() +
  labs(x = "mutation type", y = "standard error (2021)") +
  theme(
    panel.grid.major.x = element_blank(),
    panel.border = element_blank(),
    axis.ticks.x = element_blank()
  )

se.boxplot <- se.2021 + se.2023
ggsave(file = "../Tileseq_Scores/GDI1/Score/output/se_boxplot.png", 
       dpi = 700, width = 9, height = 7)


# scatter plot for scores, sep by mutation type
combined_mutation <- ggplot(match_df, aes(x = score.2023, y = score.2021)) + 
  geom_point(aes(color=type, alpha=0.9)) +
  geom_text_repel(aes(label = plotname),size = 1.8,max.overlaps = 8) +
  # geom_smooth(method = "lm") +
  # stat_cor(method= "spearman", label.x = 2.5, label.y = -2, cor.coef.name = "rho",aes(after_stat(r.label))) +
  geom_rug(aes(color=type))+
  stat_cor(method = "spearman", cor.coef.name = c("rho", "R")) +
  facet_wrap(~ type) +
  theme_bw()+
  theme(legend.position = "bottom") +
  theme(legend.key.size = unit(0.5, "cm")) 
ggsave(file = "../Tileseq_Scores/GDI1/Score/output/score_scatter.png",
       dpi = 700, width = 8, height = 6)


# sep by aa
aa_mutation <- ggplot(match_df, aes(x = score.2023, y = score.2021)) + 
  geom_point(aes(color=wt)) +
  geom_text_repel(aes(label = plotname),size = 1.8,max.overlaps = 8) +
  # geom_smooth(method = "lm") +
  # stat_cor(method= "spearman", label.x = 2.5, label.y = -2, cor.coef.name = "rho",aes(after_stat(r.label))) +
  # geom_rug(aes(color=wt))+
  facet_wrap(~ wt) +
  theme_bw()+
  theme(legend.position = "bottom") +
  theme(legend.key.size = unit(0.5, "cm"))
ggsave(file = "../Tileseq_Scores/GDI1/Score/output/aa_scatter.png",
       dpi = 700, width = 10, height = 8)

# try convert old score file to mave db format
map_db_2021 <- read.csv("GDI1/Score/data/score2021.csv") %>% 
  rename("hgvs_pro" = hgvsp) %>% 
  select(hgvs_pro, score, sd, se)
write.csv(map_db_2021, "GDI1/Score/data/score_db_2021.csv", row.names = FALSE)
