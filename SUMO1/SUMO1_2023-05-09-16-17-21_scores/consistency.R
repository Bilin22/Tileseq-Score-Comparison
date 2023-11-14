# load required packages
library(tidyverse)
library(mavevis)
library(ggprism)
library(ggpubr)


# parameters needed: (pdb.acc) PDB accession ,
# (main.chain) chain identifier
# use histogram to manually pick the threshold for burial
# default value of threshold = 0.1 # is it make sense to have the default value??
# (score) a string indicates the score file


# consistency_check(pdb.acc = "1KPS",
# main.chain = "A", threshold = 0.1, 
# score_file = "~/Desktop/Github Projects/Tileseq_Scores/UBE2I/UBE2I_2023-05-26-14-05-02_scores/select_t1_simple_aa.csv")


consistency_check <- function(pdb.acc, main.chain, threshold = 0.1, score_file){
  # get the df
  struc <- struc <- mavevis::calc.strucfeats(pdb.acc, main.chain) %>% 
    drop_na()
  # select columns and create indicators
  newstruc <- struc %>% 
    select(contains(c("pos", "all.rel", "rel.burial")))
  # core
  core <- subset(newstruc, all.rel < 20) %>% 
    mutate(type = "core")
  # surface
  surface <- subset(newstruc, all.rel > 50) %>% 
    mutate(type = "surface")
  
  labelled <- bind_rows(core, surface) #core + surface
  
  # looping over the column indices contains "rel.burial" 
  for (i in c(grep("rel.burial", colnames(newstruc))[1]: 
              length(as_vector(colnames(newstruc))))) {
    char <- nchar(colnames(newstruc[i]))
    name <- substr(colnames(newstruc[i]), 12, char)
    new <- newstruc %>% 
      filter(newstruc[i] >= threshold) %>% 
      mutate(type = paste("near", name))
    labelled <- bind_rows(labelled, new)
  }
  # now labelled contains all labelled data
  
  score_df <- read.csv(file = score_file ,skip = 16) %>% 
    mutate(mut=substr(hgvs_pro, nchar(hgvs_pro)-2, nchar(hgvs_pro))) %>%
    mutate(wt=substr(hgvs_pro, 3, 5)) %>% 
    mutate(type= ifelse(grepl("=", mut), "synonymous", ifelse(grepl("Ter$", mut), 
                                                              "nonsense", "missense"))) %>% 
    mutate(pos = as.numeric(ifelse(type == "synonymous", substr(hgvs_pro, 6, nchar(hgvs_pro) - 1), 
                                   substr(hgvs_pro, 6, nchar(hgvs_pro) - 3)))) %>% 
    select(pos, score, se) 
  

  # join with the labelled data on "pos" column
  withscore <- left_join(labelled, score_df)
  
  # plot the dotplot
  dot_plot <- ggplot(data = withscore, mapping = aes(x = type, y = score)) +
    geom_dotplot(binaxis = "y", stackdir = "center", dotsize=0.6, alpha = 0.6, binwidth = 1/30) +
    geom_hline(yintercept = 1, linetype = "dashed", col = "darkgreen", size = .5) +
    geom_hline(yintercept = 0, linetype = "dashed", col = "firebrick", size = .5) +
    labs(y = "functionaity score per AA") +
    theme_prism()
  
  # add median
  dot_plot_with_median <- dot_plot + 
    stat_summary(fun=median, geom="point", shape=18,
                 size=2, color="gold2")
  
  # Wilcoxon test between near* groups and surface
  p_val_df <- pairwise.wilcox.test(withscore$score, as.factor(withscore$type),
                       alternative = "greater",
                       p.adjust.method = "none") %>% 
    broom::tidy() %>% 
    mutate(p.value = round(p.value, 6))
  
  plot <- dot_plot_with_median +  stat_pvalue_manual(
    p_val_df, 
    y.position = 1.6, step.increase = 0.1,
    label = "p.value"
  )
  
  return(plot)

}

