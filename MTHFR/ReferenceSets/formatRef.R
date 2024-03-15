library(tidyverse)

setwd("~/Desktop/Github Projects/Tileseq_Scores/MTHFR/referenceSets")

# read the file, 89 carrying missense variants in both alleles
missense <- read.csv(file = "literatureRef.csv")[1:13] %>% 
  filter(Type_a1 == "missense", Type_a2 == "missense")


early <- missense %>% 
  filter(Category == "early")

early_a1 <- early %>% 
  select(Sample, Category, Nucleotide_change_a1, Protein_consequence_a1)

early_a2 <- early %>% 
  select(Sample, Category, Nucleotide_change_a2, Protein_consequence_a2)

late <- missense %>% 
  filter(Category == "late")

late_a1 <- late %>% 
  select(Sample, Category, Nucleotide_change_a1, Protein_consequence_a1)

late_a2 <- late %>% 
  select(Sample, Category, Nucleotide_change_a2, Protein_consequence_a2)

# turn early and late into factor
factor(match_a1$Category)
factor(match_a2$Category)

match_a1 <- rbind.data.frame(early_a1, late_a1) %>% 
  arrange(Nucleotide_change_a1)

match_a2 <- rbind.data.frame(early_a2, late_a2) %>% 
  arrange(Nucleotide_change_a2)

# for (c in distinct(match_a1 %>% select(Nucleotide_change_a1))) {print(c)} 
# all unique changes

# vector contains all unique variants
codon_a1 <- unique(match_a1$Nucleotide_change_a1)

codon_a2 <- unique(match_a2$Nucleotide_change_a2)

category_df_a1 <- data.frame(matrix(ncol = 4, nrow = 0))
colnames(category_df_a1) <- c("hgvsc","hgvsp" ,"early", "late")

for (c in codon_a1){
  subset <- match_a1 %>% 
    filter(Nucleotide_change_a1 == c)
  early_count <- as.numeric(table(subset$Category)["early"])
  late_count <- as.numeric(table(subset$Category)["late"])
  hgvsp <- as.character(subset$Protein_consequence_a1[1])
  rbind.data.frame(category_df_a1, data.frame("hgvsc" = c, "hgvsp" = hgvsp,"early" = early_count, "late" = late_count)) -> category_df_a1
}


category_df_a2 <- data.frame(matrix(ncol = 4, nrow = 0))
colnames(category_df_a2) <- c("hgvsc","hgvsp" ,"early", "late")

for (c in codon_a2){
  subset <- match_a2 %>% 
    filter(Nucleotide_change_a2 == c)
  early_count <- as.numeric(table(subset$Category)["early"])
  late_count <- as.numeric(table(subset$Category)["late"])
  hgvsp <- as.character(subset$Protein_consequence_a2[1])
  rbind.data.frame(category_df_a2, data.frame("hgvsc" = c, "hgvsp" = hgvsp,"early" = early_count, "late" = late_count)) -> category_df_a2
}

both_alleles <- full_join(category_df_a1, category_df_a2, by = "hgvsc")



# sub NA to 0
both_alleles[is.na(both_alleles)] <- 0

both_alleles$category <- NA


both_alleles$early_sum <- both_alleles$early.x + both_alleles$early.y
both_alleles$late_sum <- both_alleles$late.x + both_alleles$late.y

# Identify rows where (early.x + early.y) > (late.x + late.y)
positive_rows <- both_alleles$early_sum >= both_alleles$late_sum

both_alleles$category[positive_rows] <- "Positive"

early_onset_positive_ref <- both_alleles %>% 
  filter(category == "Positive")

write.csv(early_onset_positive_ref, "early_onset_positive_ref.csv", row.names = FALSE)

# the gnomad negative ref sets from the script
gnomADRef <- read.csv("MTHFR_refVars.csv") %>% 
  filter(referenceSet == "Negative")

# change mis-annotated values
write.csv(gnomADRef, "NegativeRef.csv", row.names = FALSE)

# newer version 
negativeRef <- read.csv("NegativeRef.csv")
positiveRef <- read.csv("PositiveRef.csv") %>% 
  select(hgvsc, hgvsp, category) %>% 
  rename("referenceSet"= "category")
  
Ref <- full_join(negativeRef, positiveRef, by = "hgvsc")

write.csv(Ref, "combined_ref.csv", row.names = FALSE)
