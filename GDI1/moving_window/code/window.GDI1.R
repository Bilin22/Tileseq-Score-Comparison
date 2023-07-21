library(ggplot2)
library(dplyr)

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
  mutate(plotname = ifelse(abs(score.2023 - score.2021) >= 0.4, hgvs_pro, "")) %>% 
  mutate(position = ifelse(type == "synonymous", substr(hgvs_pro, 6, nchar(hgvs_pro) - 1), 
                           substr(hgvs_pro, 6, nchar(hgvs_pro) - 3)))

# select 3 columns
df <- match_df %>% 
  select(position, score.2023, score.2021) %>% 
  mutate(position = as.numeric(position))

# prepare the columns for moving window analysis
ws <- 30

new.columns <- t(sapply(df$position, function(v) {
  relevant <- which(abs(df$position - v) < ws / 2)
  score_2023 <- df$score.2023[relevant]
  score_2021 <- df$score.2021[relevant]
  rho <- cor(df[relevant,c("score.2023","score.2021")],method="spearman")[1,2]
  return(c(position = v, score2023.mean = mean(score_2023), score2021.mean = mean(score_2021),
           rho = rho))
}))

# data frame for plotting
window.df <- as.data.frame(new.columns)


# moving window for fitness score
# set the color manually
colors <- c("score2023" = "blue3", "score2021" = "red3")
score.window <- ggplot(window.df, aes(x = position)) +
  geom_line(aes(y = score2023.mean, color = "score2023")) +
  geom_line(aes(y = score2021.mean, color = "score2021")) +
  labs(title = "Average Fitness Scores of GDI1", x = "position", y = "mean fitness score", color = "Legend") +
  scale_color_manual(values = colors) +
  theme_bw() +
  theme(legend.position = c(0.7, 0.9), 
        legend.background = element_rect(fill = "white", colour = "grey22"))

score.window + annotate("rect", xmin = 40, xmax = 90, ymin = 0, ymax = 1, alpha = .2, fill = "grey66") # annotation , added a shaded rect
# score.window + guides(color = guide_legend(reverse = TRUE))
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