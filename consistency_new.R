# load required packages
library(tidyverse)
library(mavevis)
library(ggprism)
library(ggpubr)
# library(beeswarm)
library(ggbeeswarm)


# parameters needed: (pdb.acc) PDB accession ,
# (main.chain) chain identifier
# use histogram to manually pick the (threshold) for burial
# default value of threshold = 0.1 # is it make sense to have the default value??
# (score) a string indicates the score file location

# example

# consistency_check(pdb.acc = "1KPS",
# main.chain = "A",
# score.file = "~/Desktop/Github Projects/Tileseq_Scores/UBE2I/UBE2I_2023-05-26-14-05-02_scores/select_t1_simple_aa.csv")

# core_threshold
# surface_threshold

consistency_check <- function(pdb.acc, main.chain, 
                              surface.threshold = 50, 
                              core.threshold = 20, 
                              burial.threshold = 0.1, 
                              score.file){
  # get the df
  struc <- mavevis::calc.strucfeats(pdb.acc, main.chain) %>% 
    drop_na()
  # select columns and create indicators
  newstruc <- struc %>% 
    select(contains(c("pos", "all.rel", "rel.burial"))) %>% 
    mutate(core = ifelse(all.rel < core.threshold, 1, 0)) %>% 
    mutate(surface = ifelse(all.rel > surface.threshold, 1, 0))
    
  for (i in c(grep("rel.burial", names(newstruc)))) {
    x <- ifelse(newstruc[i] >= burial.threshold, 1, 0)
    colnames(x) <- paste("near.",substr(names(newstruc)[i], 12, nchar(names(newstruc)[i])), sep = "")
    newstruc <- cbind(newstruc, x)
  }
  
  # read the score df
  score_df <- read.csv(file = score.file ,skip = 16) %>% 
    mutate(mut=substr(hgvs_pro, nchar(hgvs_pro)-2, nchar(hgvs_pro))) %>%
    mutate(wt=substr(hgvs_pro, 3, 5)) %>% 
    mutate(type= ifelse(grepl("=", mut), "synonymous", ifelse(grepl("Ter$", mut), 
                                                              "nonsense", "missense"))) %>% 
    mutate(pos = as.numeric(ifelse(type == "synonymous", substr(hgvs_pro, 6, nchar(hgvs_pro) - 1), 
                                   substr(hgvs_pro, 6, nchar(hgvs_pro) - 3)))) %>% 
    select(pos, score) 
  
  # check 
  # join with the labelled data on "pos" column
  withscore <- merge(newstruc, score_df, by = "pos", all = TRUE)

  core <- subset(withscore, core == 1) %>% 
    select(pos, score) %>% 
    mutate(type = "core")
  
  surface <- subset(withscore, surface == 1) %>% 
    select(pos, score) %>% 
    mutate(type = "surface")
  
  final <- rbind(core, surface)
  
  for (i in c(grep("near", names(withscore)))){
    print(names(withscore)[i])
    n <- subset(withscore, withscore[,i] == 1) %>% 
      select(pos, score) %>% 
      mutate(type = substr(names(withscore)[i], 6, nchar(names(withscore)[i])))
    final <- rbind(final, n)
  }
  
  compare_means(score~type, data = final, 
                ref.group = "surface", method = "wilcox.test")
  # plot the dotplot
  beeswarm_plot <- ggplot(data = final, mapping = aes(x = type, y = score)) +
    geom_beeswarm(cex = 1.1,na.rm = TRUE, aes(color = type), show.legend = FALSE) +
    # geom_dotplot(binaxis = "y", stackdir = "center", dotsize=0.6, alpha = 0.6, binwidth = 1/30) +
    geom_hline(yintercept = 1, linetype = "dashed", col = "darkgreen", size = .5) +
    geom_hline(yintercept = 0, linetype = "dashed", col = "firebrick", size = .5) +
    stat_summary(fun.y = median, 
                 fun.ymin = median, 
                 fun.ymax = median, 
                 geom = "crossbar", 
                 width = 0.4) +
    stat_compare_means(label = "p.signif", 
                        method = "wilcox.test", ref.group = "surface")+
    labs(y = "functionaity score per AA") +
    theme(legend.position = "none")+
    theme_prism()
  

  

  return(beeswarm_plot)
  
}


