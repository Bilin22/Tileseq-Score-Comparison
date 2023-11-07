library(tidyverse)
library(mavevis)
library(ggprism)
library(ggpubr)

# setwd
setwd("~/Desktop/Github Projects/Tileseq_Scores/SUMO1/SUMO1_2023-05-09-16-17-21_scores/")

# group into surface/core
struc <- mavevis::calc.strucfeats("3UIP", "B") %>% 
  drop_na()

newstruc <- struc %>% 
  select(res, pos, all.rel, rel.burial.UBC9, rel.burial.RAGP1, rel.burial.RBP2) %>% 
  mutate(InCore = ifelse(all.rel < 20, 1, 0)) %>% 
  mutate(InSurface = ifelse(all.rel > 50, 1, 0)) %>% 
  mutate(InUBC9 = ifelse(rel.burial.UBC9 >= 0.2, 1, 0)) %>% 
  mutate(InRBP2 = ifelse(rel.burial.RBP2 >= 0.2, 1, 0))

# histogram for each burials to pick the threshold
hist(newstruc$rel.burial.UBC9, main = "Histogram for UBC9 burial", 
     xlab = "burial")
# threshold -- 0.2

# hist(newstruc$rel.burial.RAGP1, main = "Histogram for RAGP1 burial", 
#      xlab = "burial")
# all zeros

# 
hist(newstruc$rel.burial.RBP2, main = "Histogram for UBC9 burial", 
     xlab = "burial")
# threshold -- 0.2

# amino acid in the core
core <- newstruc %>% 
  filter(InCore == 1) %>% 
  select(res, pos) %>% 
  mutate(type = "core")

# amino acid on the surface
surface <- newstruc %>% 
  filter(InSurface == 1) %>% 
  select(res, pos) %>% 
  mutate(type = "surface")

nearUBC9 <- newstruc %>% 
  filter(InUBC9 == 1) %>% 
  select(res, pos) %>% 
  mutate(type = "near UBC9")

nearRBP2 <- newstruc %>% 
  filter(InRBP2 == 1) %>% 
  select(res, pos) %>% 
  mutate(type = "near RBP2")

# labelled data
new <- bind_rows(core, surface, nearRBP2, nearUBC9)

# scores
score <- read.csv(file = "select_t1_simple_aa.csv",skip = 16) %>% 
  mutate(mut=substr(hgvs_pro, nchar(hgvs_pro)-2, nchar(hgvs_pro))) %>%
  mutate(wt=substr(hgvs_pro, 3, 5)) %>% 
  mutate(type= ifelse(grepl("=", mut), "synonymous", ifelse(grepl("Ter$", mut), 
                                                            "nonsense", "missense"))) %>% 
  mutate(pos = as.numeric(ifelse(type == "synonymous", substr(hgvs_pro, 6, nchar(hgvs_pro) - 1), 
                           substr(hgvs_pro, 6, nchar(hgvs_pro) - 3)))) %>% 
  select(pos, score, se) 

# and pool the fitness score for these columns -- median
# pooled <- score %>% 
#   group_by(pos) %>% 
#   mutate(score_median = median(score), se_median = median(se)) %>% 
#   select(pos, score_median, se_median) %>% 
#   distinct()

# join with the labelled data on "pos" column
withscore <- left_join(new, score)


# get the first and the third quantile, the median of the data
summarywithscore <-withscore %>% 
  group_by(type) %>% 
  summarise(median = median(score), 
            firstquantile = quantile(score, probs = c(0.25)), 
            thirdquantile = quantile(score, probs = c(0.75)))

# plot the dot plot
dot_graph <- ggplot(data = withscore, mapping = aes(x = type, y = score)) +
  geom_dotplot(binaxis = "y", stackdir = "center", dotsize=0.6, alpha = 0.6, binwidth = 1/30) +
  geom_hline(yintercept = 1, linetype = "dashed", col = "darkgreen", size = .5) +
  geom_hline(yintercept = 0, linetype = "dashed", col = "firebrick", size = .5) +
  labs(y = "functionaity score per AA") +
  theme_prism()


with_line <- dot_graph + 
  stat_summary(fun=median, geom="point", shape=18,
               size=2, color="gold2")
  
# ggsave(file = "comparison_overall.png",
#        dpi = 700, width = 6, height = 5)

# get the median of the withscore data

rbp2 <- as.numeric(unlist(withscore %>% 
  filter(type == "near RBP2") %>% 
  select(score)))

ubc9 <- as.numeric(unlist(withscore %>% 
                            filter(type == "near UBC9") %>% 
                            select(score)))

surface <- as.numeric(unlist(withscore %>% 
                            filter(type == "surface") %>% 
                            select(score)))

# wilcox.test() on the score data
RBP2vsSurface <- round(wilcox.test(x = rbp2, y = surface, 
                             alternative = "less")$p.value, 7)

UBC9vsSurface <- round(wilcox.test(x = ubc9, y = surface, 
                             alternative = "less")$p.value, 4)



stat.test <- tibble::tribble(~group1, ~group2, ~p.adj,
                            "near RBP2", "surface", RBP2vsSurface,
                            "near UBC9", "surface", UBC9vsSurface)
with_line +  stat_pvalue_manual(
  stat.test, 
  y.position = 1.6, step.increase = 0.1,
  label = "p.adj"
)

ggsave(file = "sumo_struc.png",
       dpi = 700, width = 6, height = 5)

# make it into a function