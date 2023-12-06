library(tidyverse)
library(ggplot2)
library(ggprism)
library(patchwork)

score_df <- read.csv(file = 
                    "~/Desktop/Github Projects/Tileseq_Scores/SUMO1/Score_comparison/data/urn_mavedb_00000001-b-2_scores.csv")[, c(4, 5)] %>%
  distinct() %>% 
  mutate(mut=substr(hgvs_pro, nchar(hgvs_pro)-2, nchar(hgvs_pro))) %>%
  mutate(wt=substr(hgvs_pro, 3, 5)) %>% 
  mutate(type= ifelse(grepl("=", mut), "synonymous", ifelse(grepl("Ter$", mut), 
                                                            "nonsense", "missense"))) %>% 
  mutate(position = ifelse(type == "synonymous", substr(hgvs_pro, 6, nchar(hgvs_pro) - 1), 
                           substr(hgvs_pro, 6, nchar(hgvs_pro) - 3)))
# nonsyn <- score_df %>% 
#   filter(type == "nonsense" | type == "synonymous")
# 
# mis <- score_df %>% 
#   filter(type == "missense")
# 
# 
# cols <- c("missense" = "grey50")
# col2 <- c("nonsense" = "salmon", "synonymous" = "darkgreen")
# p <- ggplot(data = nonsyn, mapping = aes(x = score)) +
#   geom_histogram(bins = 40, aes(fill = type), position = "dodge") +
#   theme_minimal() +
#   theme(legend.position = "bottom") 
# p2 <- p + scale_fill_manual(values = col2)
# 
# p3 <- ggplot(data = mis, mapping = aes(x = score)) + 
#   geom_histogram(bins = 40, aes(fill = type), position = "dodge") +
#   theme_minimal() +
#   scale_y_reverse()+
#   theme(legend.position = "bottom") 
# p4 <- p3 + scale_fill_manual(values = cols)
  

score_df <- score_df %>%
  mutate(reflection = ifelse(type == "missense", -1, 1)) 

# # add cohen's
nonsense_mean <- mean(score_df[which(score_df$type == "nonsense"), 2])
synonymous_mean <- mean(score_df[which(score_df$type == "synonymous"), 2])
nonsense_var <- var(score_df[which(score_df$type == "nonsense"), 2])
synonymous_var <- var(score_df[which(score_df$type == "synonymous"), 2])
d_score <- round((synonymous_mean - nonsense_mean) / 
  sqrt((synonymous_var + nonsense_var) /2 ), 4)



p <- ggplot(score_df, aes(x = score, fill = type)) +
  geom_histogram(data = score_df %>% filter(type != "missense"), 
                 aes(y = ..density..), bins = 40, position = "dodge") +
  geom_histogram(data = score_df %>% filter(type == "missense"), 
                 aes(y = -..density..), bins = 40, position = "identity") +
  scale_fill_manual(values = c("nonsense" = "salmon", "synonymous" = "darkgreen", "missense" = "grey")) +
  geom_text(x = 0, y = 2, label = paste("Cohen's d: ", d_score))+
  theme_minimal() +
  theme(legend.position = "bottom") +
  labs(y = "Density") +
  scale_y_continuous(labels = abs, breaks = scales::pretty_breaks(n = 10))  # Use absolute values for y-axis labels

# ggsave("SUMO1_old_separation.pdf", height = 5, width = 7)


