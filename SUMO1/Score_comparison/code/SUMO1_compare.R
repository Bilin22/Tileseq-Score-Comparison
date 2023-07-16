---
title: "Score_Comparison"
author: "Bilin"
date: "`r Sys.Date()`"
output: pdf_document
---

```{r, echo=FALSE, results=FALSE}
library(ggplot2)
library(dplyr)
library(stringr)
library(ggpubr)
library(tidyr)
library(patchwork)
```


```{r new}
df_new <- read.csv(file = "../data/select_t1_simple_aa.csv", skip = 16) %>% 
  select(hgvs_pro, score, se) %>% 
  mutate(mut=substr(hgvs_pro, nchar(hgvs_pro)-2, nchar(hgvs_pro))) %>%
  mutate(wt=substr(hgvs_pro, 3, 5)) %>% 
  mutate(type= ifelse(grepl("=", mut), "synonymous", ifelse(grepl("Ter$", mut), 
                                                            "nonsense", "missense")))

```

```{r old}
df_old <- read.csv(file = "../data/urn_mavedb_00000001-b-2_scores.csv") %>% 
  select(hgvs_pro, score, se) %>% 
  distinct() %>% 
  mutate(mut=substr(hgvs_pro, nchar(hgvs_pro)-2, nchar(hgvs_pro))) %>%
  mutate(wt=substr(hgvs_pro, 3, 5)) %>% 
  mutate(type= ifelse(grepl("=", mut), "synonymous", ifelse(grepl("Ter$", mut), 
                                                            "nonsense", "missense")))

```

```{r matching}
matched_df <- merge(df_new, df_old, by = "hgvs_pro") %>% 
  rename(new.score = score.x, new.se = se.x, old.score = score.y, old.se = se.y, 
         mut = mut.x, wt = wt.x, type = type.x) %>% 
  select(hgvs_pro, new.score, old.score, new.se, old.se, wt, mut, type)

```

```{r}
# plot out the se for both map

```





```{r}
new_old_combined <- ggplot(matched_df, aes(x = new.score, y = old.score)) + 
  geom_point(aes(color=type, alpha=0.9)) +
  # geom_smooth(method = "lm") +
  stat_cor(method= "spearman", label.x = -0.5, label.y = 1, cor.coef.name = "rho",aes(label = ..r.label..)) +
  geom_rug(aes(color=type))+
  xlab("score.2023") +
  ylab("score.2018") +
  theme_bw()+
  theme(legend.position = "bottom") +
  theme(legend.key.size = unit(0.5, "cm")) 


new_old_facet <- ggplot(matched_df, aes(x = new.score, y = old.score)) + 
  geom_point(aes(color=type, alpha = 0.9)) +
  facet_wrap(~type, scales = "fixed") +
  # geom_smooth(method = "lm", aes(color=wildtype))+
  stat_cor(method= "spearman", label.x = 0.4, label.y = -0.7, cor.coef.name = "rho", aes(label = ..r.label..)) +
  geom_rug(aes(color=type))+
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45), legend.position = "none")

new_old_facet + new_old_combined

# ggsave(filename = "../output/SUMO1_compare_mut.jpg", dpi = 600, width = 13, height = 6)
```

```{r}
new_old_combined_aa <- ggplot(matched_df, aes(x = new.score, y = old.score)) + 
  geom_point(aes(color=wt, alpha=0.8)) +
  # geom_smooth(method = "lm") +
  stat_cor(method= "spearman", label.x = -0.5, label.y = 1, cor.coef.name = "rho",aes(label = ..r.label..)) +
  geom_rug(aes(color=wt))+
  theme_bw()+
  theme(legend.position = "bottom") +
  theme(legend.key.size = unit(0.5, "cm")) 


new_old_facet_aa <- ggplot(matched_df, aes(x = new.score, y = old.score)) + 
  geom_point(aes(color=wt, alpha = 0.8)) +
  facet_wrap(~wt, scales = "fixed") +
  # geom_smooth(method = "lm", aes(color=wildtype))+
  stat_cor(method= "spearman", label.x = 0, label.y = -1, cor.coef.name = "rho", aes(label = ..r.label..)) +
  geom_rug(aes(color=wt))+
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45), legend.position = "none")

new_old_facet_aa + new_old_combined_aa

# ggsave(filename = "../output/SUMO1_compare_aa.jpg", dpi = 600, width = 14, height = 6)
```

