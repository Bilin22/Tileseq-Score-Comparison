# load required packages
library(tidyverse)
library(mavevis)
library(ggbeeswarm)
library(argparser)



# parameters needed: (pdb.acc) PDB accession ,
# (main.chain) chain identifier
# use histogram to manually pick the (threshold) for burial
# default value of threshold = 0.1 # is it make sense to have the default value??
# (score) a string indicates the score file location


# consistency_check("3UIP", "B",
#                   map =
#                     "~/Desktop/Github Projects/Tileseq_Scores/SUMO1/SUMO1_2023-05-09-16-17-21_scores/select_t1_simple_aa.csv")
# # example

# consistency_check(pdb.acc = "1KPS",
# main.chain = "A",
# map = "~/Desktop/Github Projects/Tileseq_Scores/UBE2I/UBE2I_2023-05-26-14-05-02_scores/select_t1_simple_aa.csv")

# core_threshold
# surface_threshold

# Rscript biochemistry.R
# "~/Desktop/Github Projects/Tileseq_Scores/SUMO1/SUMO1_2023-05-09-16-17-21_scores/select_t1_simple_aa.csv"
# --pdb "3UIP" --chain "B"


#process command line arguments
p <- arg_parser(
  "Generate a beeswarm plot for biochemistry consistency check",
  name = "biochemistry.R"
)
p <- add_argument(p, "--name", help = "name of the gene")
p <- add_argument(p, "--pdb", help = "a string represents PDB accession")
p <- add_argument(p, "--chain", help = "a string represents chain identifier")
p <- add_argument(p, "--surfaceThreshold", 
                  help = "greater than the threshold amount (in percentage) 
                  of surface area accessible to solvent, default = 50", 
                  default = 50)
p <- add_argument(p, "--coreThreshold", help = "less than the 
                  threshold amount (in percentage) of surface area
                  accessible to solvent, default = 20", 
                  default = 50)
p <- add_argument(p, "--burialThreshold", help = "threshold for burial, 
                  default = 0.1", default = 0.1)
p <- add_argument(p, "map", help = "the location of map file in MAVEdb format.")
p <- add_argument(p, "--outputFile", help = "output PDF file name", 
                  default = "beeswarm_plot.pdf")
args <- parse_args(p)



consistency_check <- function(name, pdb.acc, main.chain, 
                              surface.threshold = 50, 
                              core.threshold = 20, 
                              burial.threshold = 0.1, 
                              map, 
                              output.file = 'beeswarm_plot.pdf'){
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
  print(paste("Reading score file:", map))
  score_df <- read.csv(file = map ,skip = 16) %>%
    mutate(mut=substr(hgvs_pro, nchar(hgvs_pro)-2, nchar(hgvs_pro))) %>%
    mutate(wt=substr(hgvs_pro, 3, 5)) %>% 
    mutate(type= ifelse(grepl("=", mut), "synonymous", ifelse(grepl("Ter$", mut), 
                                                              "nonsense", "missense"))) %>% 
    mutate(pos = as.numeric(ifelse(type == "synonymous", substr(hgvs_pro, 6, nchar(hgvs_pro) - 1), 
                                   substr(hgvs_pro, 6, nchar(hgvs_pro) - 3)))) %>% 
    select(pos, score) 
  # score_df <- read.csv(file = map) %>% 
      
  
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
    n <- subset(withscore, withscore[,i] == 1) %>% 
      select(pos, score) %>% 
      mutate(type = substr(names(withscore)[i], 6, nchar(names(withscore)[i])))
    final <- rbind(final, n)
  }
  
  # compare_means(score~type, data = final,
  #               ref.group = "surface", method = "wilcox.test")
  # plot the beeswarm plot
  beeswarm_plot <- ggplot(data = final, mapping = aes(x = type, y = score)) +
    geom_beeswarm(cex = 0.9,na.rm = TRUE, aes(color = type),
                  show.legend = FALSE) +
    geom_hline(yintercept = 1, linetype = "dashed", col = "darkgreen", size = .5) +
    geom_hline(yintercept = 0, linetype = "dashed", col = "firebrick", size = .5) +
    stat_summary(fun.y = median, 
                 fun.ymin = median, 
                 fun.ymax = median, 
                 geom = "crossbar", 
                 width = 0.4) +
    stat_compare_means(label = "p.signif", 
                        method = "wilcox.test", ref.group = "surface")+
    theme(legend.position = "none") +
    ggtitle(label = paste("gene:",name,
                          "chain:",main.chain, "pdb.access:",pdb.acc)) +
    theme_classic()
    # theme_light()
    # theme_prism()
  

  ggsave(output.file, beeswarm_plot, device = "pdf", 
         scale = 0.8)
}

consistency_check(args$name,args$pdb.acc, args$main.chain,
                  args$surfaceThreshold, args$coreThreshold,
                  args$burialThreshold,
                  args$map, args$outputFile)

# consistency_check(args$name,args$pdb.acc, args$main.chain, 
#                   args$surface.threshold, args$core.threshold, 
#                   args$burial.threshold, 
#                   args$map, args$output.file)

