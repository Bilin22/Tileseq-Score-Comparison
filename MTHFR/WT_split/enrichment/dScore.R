library(tidyverse)
#set working directory
setwd("/Users/admin/Desktop/Github Projects/Tileseq_Scores/MTHFR/WT_split/enrichment")

# read the dataset: f12WT
f12WT <- read.csv(file = "f12WT_t1_enrichment.csv", skip = 17) %>% 
  select(hgvsp, type, logPhi, logPhi.se) %>% 
  filter(is.finite(logPhi))


# nonsense
f12WT_nonsense <- f12WT %>% 
  filter(type == "nonsense")
f12WT_nonsense_mean <- mean(f12WT_nonsense$logPhi)
f12WT_nonsense_variance <- var(f12WT_nonsense$logPhi)

# synonymous
f12WT_synonymous <- f12WT %>% 
  filter(type == "synonymous")
f12WT_synonymous_mean <- mean(f12WT_synonymous$logPhi)
f12WT_synonymous_variance <- var(f12WT_synonymous$logPhi)

# Cohen's d for Welch test
f12WT_dscore <- (f12WT_synonymous_mean - f12WT_nonsense_mean) / 
  sqrt((f12WT_synonymous_variance + f12WT_nonsense_variance) /2 )


# a function for calculating the d scoress
# enrichment data should be in ''
# need tidyverse/ dplyr
calc_dScore <- function(enrichment_data){
  cleaned_df <- read.csv(enrichment_data, skip = 17) %>% 
    select(hgvsp, type, logPhi, logPhi.se) %>% 
    filter(is.finite(logPhi))
  synonymous_df <- cleaned_df %>% 
    filter(type == "synonymous")
  nonsense_df <- cleaned_df %>%
    filter(type == "nonsense")
  synonymous_mean <- mean(synonymous_df$logPhi)
  nonsense_mean <- mean(nonsense_df$logPhi)
  synonymous_var <- var(synonymous_df$logPhi)
  nonsense_var <- var(nonsense_df$logPhi)
  return( (synonymous_mean - nonsense_mean)/
            sqrt((synonymous_var + nonsense_var)/2)
    
  ) 
  
}