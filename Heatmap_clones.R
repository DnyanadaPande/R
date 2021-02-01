#
# Dnyanada Pande, June 2020
# Compare the clones found in the control files and make a heatmap.
# 

# Load the required packages.
library(ggplot2)
library(pheatmap)

dir <- "/Volumes/dnyanada/final_files/"
setwd(dir)

std1 <- read.table("standard124.tsv", header = T, stringsAsFactors = F, sep = '\t')
std1 <- std1[1:100,]
head(std1)
clone1 <- std1$id
ct1 <- rep(1,length(clone1))
std124 <- cbind.data.frame(clone1,ct1)
colnames(std124) <- c("clone","chip124")

std2 <- read.table("standard125.tsv", header = T, stringsAsFactors = F, sep = '\t')
std2 <- std2[1:100,]
head(std2)
clone2 <- std2$id
ct2 <- rep(1,length(clone2))
std125 <- cbind.data.frame(clone2,ct2)
colnames(std125) <- c("clone","chip125")

# Merge the 2 files by the clone id and make a data frame that indicates 
# the presence (1) or absence (0) of the clone in the 2 controls.
df <- merge(std124, std125, by = 'clone', all.x = T, all.y = T,)
df[is.na(df)] <- 0

df <- as.data.frame(df)
df2 <- df[,-1]
rownames(df2) <- df[,1]
df2 <- as.matrix(df2)
write.csv(df, "standard_clones.csv")

# Heatmap
pn <- paste0(dir,"Standard_top100.pdf")
pdf(pn)
pheatmap(df2, color = c("lightsteelblue1", "steelblue3"), show_rownames = T, 
         border_color = c("black"), cluster_rows = F, cluster_cols = F,
         legend = T, legend_breaks = c(0,1), fontsize_row = 5,
         fontsize_col = 10, main = "Z13264 standard - Top 100")
dev.off()

