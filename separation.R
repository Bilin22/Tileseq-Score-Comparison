library(tidyverse)
# library(ggplot2)


# score_df <- read.csv(file = 
#                        "~/Desktop/Github Projects/Tileseq_Scores/CALM1/Score_comparison/data/urn_mavedb_00000001-c-1_scores.csv")[, c(4, 5)] %>%
score_df <- read.csv(file =
                    "~/Desktop/Github Projects/Tileseq_Scores/SUMO1/Score_comparison/data/urn_mavedb_00000001-b-2_scores.csv")[, c(4, 5)] %>%
  distinct() %>% 
  mutate(mut=substr(hgvs_pro, nchar(hgvs_pro)-2, nchar(hgvs_pro))) %>%
  mutate(wt=substr(hgvs_pro, 3, 5)) %>% 
  mutate(type= ifelse(grepl("=", mut), "synonymous", ifelse(grepl("Ter$", mut), 
                                                            "nonsense", "missense"))) %>% 
  mutate(position = ifelse(type == "synonymous", substr(hgvs_pro, 6, nchar(hgvs_pro) - 1), 
                           substr(hgvs_pro, 6, nchar(hgvs_pro) - 3)))


 # add cohen's d
nonsense_mean <- mean(score_df[which(score_df$type == "nonsense"), 2])
synonymous_mean <- mean(score_df[which(score_df$type == "synonymous"), 2])
nonsense_var <- var(score_df[which(score_df$type == "nonsense"), 2])
synonymous_var <- var(score_df[which(score_df$type == "synonymous"), 2])
d_score <- round((synonymous_mean - nonsense_mean) / 
  sqrt((synonymous_var + nonsense_var) /2 ), 4)


# set y -axis abs value, increase the breaks
ggplot(score_df, aes(x = score, fill = type)) +
  geom_histogram(data = score_df %>% filter(type != "missense"), 
                 aes(y = ..density..), bins = 40, position = "dodge") +
  geom_histogram(data = score_df %>% filter(type == "missense"), 
                 aes(y = -..density..), bins = 40, position = "identity") +
  scale_fill_manual(values = 
                      c("nonsense" = "salmon", 
                        "synonymous" = "darkgreen", 
                        "missense" = "grey")) +
  geom_text(x = mean(range(score_df$score)), 
            y = 2, label = paste("Cohen's d: ", d_score))+
  theme_minimal() +
  theme(legend.position = c(0, 1), legend.justification = c(0, 1)) +
  guides(fill = guide_legend(title = NULL))+
  labs(y = "Density") +
  scale_y_continuous(labels = abs, 
                     breaks = scales::pretty_breaks(n = 10))  

ggsave("SUMO_old_separation.pdf", height = 5, width = 7)


