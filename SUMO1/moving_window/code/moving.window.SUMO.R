library(ggplot2)
library(dplyr)
library(stringr)
library(ggpubr)
library(tidyr)
library(patchwork)
library(ggrepel)

# old match_df
map.2019 <- read.csv(file = "../Tileseq_Scores/SUMO1/moving_window/data/urn_mavedb_00000001-b-2_scores.csv") %>% 
  select(hgvs_pro, score, se) %>% 
  distinct() %>% 
  mutate(mut=substr(hgvs_pro, nchar(hgvs_pro)-2, nchar(hgvs_pro))) %>%
  mutate(wt=substr(hgvs_pro, 3, 5)) %>% 
  mutate(type= ifelse(grepl("=", mut), "synonymous", ifelse(grepl("Ter$", mut), 
                                                            "nonsense", "missense"))) %>% 
  mutate(position = ifelse(type == "synonymous", substr(hgvs_pro, 6, nchar(hgvs_pro) - 1), 
                           substr(hgvs_pro, 6, nchar(hgvs_pro) - 3)))

# new match_df, add a new column for amino acid position
map.2023 <- read.csv(file = "../Tileseq_Scores/SUMO1/moving_window/data/select_t1_simple_aa.csv", 
                     skip = 16) %>% 
  mutate(mut=substr(hgvs_pro, nchar(hgvs_pro)-2, nchar(hgvs_pro))) %>%
  mutate(wt=substr(hgvs_pro, 3, 5)) %>% 
  mutate(type= ifelse(grepl("=", mut), "synonymous", ifelse(grepl("Ter$", mut), 
                                                            "nonsense", "missense"))) %>% 
  mutate(position = ifelse(type == "synonymous", substr(hgvs_pro, 6, nchar(hgvs_pro) - 1), 
                                                        substr(hgvs_pro, 6, nchar(hgvs_pro) - 3)))

# build a matched table
match_df <- merge(map.2023, map.2019, by = "hgvs_pro") %>% 
  rename(score.2023 = score.x, se.2023 = se.x, score.2019 = score.y, se.2019 = se.y, 
         mut = mut.x, wt = wt.x, type = type.x) %>% 
  select(hgvs_pro, score.2023, se.2023, score.2019, se.2019, wt, mut, type) %>% 
  mutate(plotname = ifelse(abs(score.2023 - score.2019) >= 0.4, hgvs_pro, "")) %>% 
  mutate(position = ifelse(type == "synonymous", substr(hgvs_pro, 6, nchar(hgvs_pro) - 1), 
                           substr(hgvs_pro, 6, nchar(hgvs_pro) - 3))) %>% 
  mutate(difference = abs(score.2023 - score.2019))

# select 3 columns
df <- match_df %>% 
  select(hgvs_pro, score.2023, score.2019, position) %>% 
  mutate(position = as.numeric(position))


## moving window -- Jochen's script

#assume we have a match_df.frame `match_df` containing hgvs.pro, fitness.score, and varity.score columns.
#extract the amino acid positions from the hgvs column:

#also assume we have a window size parameter:
window.size <- 2
#iterate over the range of AA positions
result <- lapply(range(df$position,na.rm=TRUE), function(pos.i) {
  #find relevant rows (within window size distance of pos.i)
  relevantRows <- which(abs((df$position) - pos.i) < window.size)
  #calculate the correlation between scores of the relevant entries
  rho <- cor(df[relevantRows,c("score.2023","score.2019")],method="spearman")[1,2]
  return(c(pos=pos.i,rho=rho))
})
#result is currently a list of vectors, which we can turn into a matrix for convenience:
result <- do.call(rbind,result)

#now we can draw that line plot
plot(result, type="l", xlab="Position", ylab="Spearman's rho")


# moving window for fitness score
ws <- 25

new.columns <- t(sapply(df$position, function(v) {
  relevant <- which(abs(df$position - v) < ws / 2)
  score_2023 <- df$score.2023[relevant]
  score_2019 <- df$score.2019[relevant]
  rho <- cor(df[relevant,c("score.2023","score.2019")],method="spearman")[1,2]
  return(c(position = v, score2023.mean = mean(score_2023), score2019.mean = mean(score_2019),
           rho = rho))
}))

window.df <- as.data.frame(new.columns)

colors <- c("score2023" = "blue3", "score2019" = "red3")

score.window <- ggplot(window.df, aes(x = position)) +
  geom_line(aes(y = score2023.mean, color = "score2023")) +
  geom_line(aes(y = score2019.mean, color = "score2019")) +
  labs(title = "Average Fitness Scores of SUMO1", x = "position", y = "mean fitness score", color = "Legend") +
  scale_color_manual(values = colors) +
  theme_bw() +
  theme(legend.position = c(0.9, 0.9), 
        legend.background = element_rect(fill = "white", colour = "grey22"))
score.window + guides(color = guide_legend(reverse = TRUE))
# ggsave(filename = "../Tileseq_Scores/SUMO1/moving_window/output/scores.png", dpi = 700,
#        height = 6, width = 8)

# moving window for spearman rho














# par(mar=c(6, 3, 3, 3) + 2)
# plot.new()

#set axis for fitness scores
plot.window(xlim = c(min(new.columns[, "position"]), max(new.columns[, "position"])),
            #change ylim here for left axis fitness
            ylim = c(0, 1),
            xlab = "AA position", ylab = "score.2023", main = "Moving Window Analysis")
#plot fitness scores
axis(side = 2, ylab = "score.2023", col = "red",col.axis = 'red', las = 2, line = 0)
lines(new.columns[, "position"], new.columns[, "score2023.mean"], col = "red", type = "l", lwd = 2)

#create new plot for stability, overlaid on the first fitness plot 
par(new = TRUE)
#set axis for stability scores
plot.window(xlim = c(min(new.columns[, "VARIANT"]), max(new.columns[, "VARIANT"])),
            #change ylime for right axis stability
            ylim = c(-1, 0),
            xlab = "AA position", ylab = "Stability", main = "Moving Window Analysis")
#Plot stability scores
axis(side = 4, ylab = "Stability", col = "blue", col.axis = 'blue', las = 2, line = 0)
lines(new.columns[, "VARIANT"], new.columns[, "stability.mean"], type = "l", lwd = 2, col = "blue")

axis(side = 1, xlab = "VARIANT", col = "black", las = 2, line = 0, at = seq(min(new.columns[, 'VARIANT']), max(new.columns[, 'VARIANT']), by = 10))
box()
#legend("bottomright", legend = c("Fitness", "Stability"), fill = c("red", "blue"), cex = 1, x.intersp = 0.2)
title(main = "Moving Window Analysis of Fitness and Stability Scores", cex.main = 1.2, cex.lab = 1.2)
# Add labels to x and y axes
xlabel <- expression("AA position")
ylabel_left <- expression(paste("Mean ", italic("Fitness")))
ylabel_right <- expression(paste("Mean ", italic("Stability")))

mtext(xlabel, side = 1, line = 3, cex = 1.2, font = 2)
mtext(ylabel_left, side = 2, line = 3, cex = 1.2, font = 2, col = 'red')
#mtext does not allow for rotation
mtext(ylabel_right, side = 4, line = 3, cex = 1.2, font = 2, col = 'blue')


#pearson correlation between moving window of fitness scores and stability scores
stats <- cor.test(new.columns[,'fitness.mean'], new.columns[,'stability.mean'], method = 'pearson')

p_value <- stats$p.value
correlation <- stats$estimate







# try plotting the score vs. amino acid position
# ggplot(map.2023, aes(x = position, y = score, color = type)) +
#   geom_line(aes(group = type)) +
#   scale_color_manual(values = c("grey59", "red2", "chartreuse3")) +
#   facet_grid(~ type)+
#   theme_bw() +
#   theme(
#     axis.text.x = element_text(angle = 60, hjust = 1)
#   )


# ggplot(map.2023, aes(x = position, y = score, color = type)) +
#   # geom_point(aes(color = type, alpha = 0.8)) +
#   geom_line(aes(group = type)) +
#   scale_color_manual(values = c("grey59", "red2", "chartreuse3")) +
#   theme_classic() +
#   theme(axis.text.x = element_text(angle = 70))


# ggplot(match_df, aes(x = position, y = difference, color = type)) +
#   # geom_point(aes(color = type, alpha = 0.8)) +
#   geom_line(aes(group = type)) +
#   scale_color_manual(values = c("grey59", "red2", "chartreuse3")) +
#   theme_classic() +
#   theme(axis.text.x = element_text(angle = 70))
