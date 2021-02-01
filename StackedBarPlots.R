#
# Dnyanada P, 07/09/2019
# Script to take an input file with an alignment sequence, frequencies, CIGAR strings and
# use pattern matching to extract specific patterns.
#
# Loading the required libraries.
library("data.table")
library("ggplot2")
library("RColorBrewer")
library("scales")
library("stringr")
library("rebus")
library("grid")

basedir <- "/Volumes/dnyanada/GE_files/"
setwd(basedir)

# Make a data frame with the file names and labels.
tests <- data.frame(c(
  "09_S9_ge.csv",
  "08_S8_ge.csv",
  "07_S7_ge.csv",
  "06_S6_ge.csv"
),c(
  "BM_D9",
  "BM_D12",
  "BM_D27",
  "BM_D36"), stringsAsFactors = F)
colnames(tests) <- c("file","label")
tests

# Initializing an empty matrix
df2 <- data.frame(matrix(ncol = 3, nrow = 0))
x <- c("name", "freq", "mut")
colnames(df2) <- x

# Calculating the frequencies
for (i in 1:length(tests$file)){   
  print(tests$file[i])
  
  ge_file <- read.csv(tests$file[i],header = TRUE)
  
  dt <- cbind.data.frame(ge_file$Frequency,ge_file$CIGAR)
  colnames(dt) <- c("Frequency","Cigar")
  print(head(dt))
  
  ge_freq  <- numeric()
  filen    <- character()
  mut_lab  <- character()
  tot_fq   <- 0.000
  d1_freq  <- 0.000
  d2_freq  <- 0.000
  d3_freq  <- 0.000
  d4_freq  <- 0.000
  d5_freq  <- 0.000
  d6_freq  <- 0.000
  d7_freq  <- 0.000
  d8_freq  <- 0.000
  d9_freq  <- 0.000
  d10_freq <- 0.000
  d11_freq <- 0.000
  d12_freq <- 0.000
  d13_freq <- 0.000
  other    <- 0.000
  
  for (j in 1:length(dt$Frequency)){ #
    # Pattern match to isolate CIGARS within the base pair range for editing and 
    # just deletions ((d/dd)Match - (d/dd)Deletion - (d/dd)Match).
    pat <- grepl("^\\dM\\d\\dD\\d\\dM$|^\\d\\dM\\d\\dD\\d\\dM$|^\\dM\\dD\\d\\dM$|^\\d\\dM\\dD\\d\\dM$|^\\d\\dM\\d\\dD\\dM$|^\\dM\\d\\dD\\dM$|^\\d\\dM\\dD\\dM$|^\\dM\\dD\\dM$",dt$Cigar[j])
    
    if (all(pat) == TRUE){
      print(dt$Cigar[j])
      tot_fq <- tot_fq + dt$Frequency[j]
      
      # Extracting the deletion sizes.
      pat2 <- str_extract(dt$Cigar[j],"\\d\\dD|\\dD")
      
      # Calculating the frequency for specific deletion pattern.
      if (pat2 == "1D" & pat2 != "11D"){
        print(pat2)
        d1_freq <- d1_freq + dt$Frequency[j]
      }
      else if (pat2 == "2D" & pat2 != "12D"){
        d2_freq <- d2_freq + dt$Frequency[j]
      }
      else if (pat2 == "3D" & pat2 != "13D"){
        d3_freq <- d3_freq + dt$Frequency[j]
      }
      else if (pat2 == "4D"){
        d4_freq <- d4_freq + dt$Frequency[j]
      }
      else if (pat2 == "5D"){
        d5_freq <- d5_freq + dt$Frequency[j]
      }
      else if (pat2 == "6D"){
        d6_freq <- d6_freq + dt$Frequency[j]
      }
      else if (pat2 == "7D"){
        d7_freq <- d7_freq + dt$Frequency[j]
      }
      else if (pat2 == "8D"){
        d8_freq <- d8_freq + dt$Frequency[j]
      }
      else if (pat2 == "9D"){
        d9_freq <- d9_freq + dt$Frequency[j]
      }
      else if (pat2 == "10D"){
        d10_freq <- d10_freq + dt$Frequency[j]
      }
      else if (pat2 == "11D"){
        d11_freq <- d11_freq + dt$Frequency[j]
      }
      else if (pat2 == "12D"){
        d12_freq <- d12_freq + dt$Frequency[j]
      }
      else if (pat2 == "13D"){
        d13_freq <- d13_freq + dt$Frequency[j]
      }
      else{
        other <- other + dt$Frequency[j]
      }
    }
    else {
      pat3 <- grepl("D|I",dt$Cigar[j]) 
      if (all(pat3) == TRUE){
        other <- other + dt$Frequency[j]
        tot_fq <- tot_fq + dt$Frequency[j]
      }
    }
  }
  
  print(d1_freq)
  print(tot_fq)
  
  ge_freq <- c(d13_freq,d12_freq,d11_freq,d10_freq,d9_freq,d8_freq,d7_freq,d6_freq,
               d5_freq,d4_freq,d3_freq,d2_freq,d1_freq,other,tot_fq)
  filen <- c(rep(tests$label[i],15))
  mut_lab <- c("13D","12D","11D","10D","09D","08D","07D","06D","05D","04D","03D","02D",
               "01D","Other","Total frequency")
  df.data <- data.frame(filen, ge_freq, mut_lab, stringsAsFactors=TRUE)
  df.data$mut_lab <- factor(df.data$mut_lab, levels = c("13D","12D","11D","10D","09D",
                   "08D","07D","06D","05D","04D","03D","02D","01D","Other","Total frequency"))
  df2 <- rbind.data.frame(df2,df.data)
}

write.csv(df2,"Deletion_frequencies_BM_July7_2020.csv")

# Plotting the frequencies and deletion patterns.
pn <- (paste0(basedir,"Plot1.pdf"))
pdf(pn)
legend_title <- "Deletion size"
df2$group <- ""
df2[1:196,]$group  <- "HSC"              
df2[197:280,]$group <- "Lineage"
ggplot(df2,aes(x=filen,y=ge_freq,fill=mut_lab)) + geom_bar(position = 'stack',stat='identity') +
  scale_y_continuous(labels = percent, limits = c(0,1), expand = c(0,0)) +
  xlab("Source") + 
  ylab("Mutation frequency") +
  labs(title = paste0("A BM"), fill = legend_title) +
  theme(axis.text.x = element_text(face="bold", color="#993333", 
                                   size=14, angle=90,hjust=0)) +
  facet_grid(~group,scales="free",space="free") + 
  scale_fill_manual(values=c("steelblue4","springgreen3","mediumpurple1","lightgoldenrod1",
                             "magenta","aquamarine","red4","palevioletred3", "yellow2",
                             "olivedrab3","lightskyblue","orange","lightpink","snow4"),
                    breaks=c("13D","12D","11D","10D","09D","08D","07D","06D","05D","04D",
                             "03D","02D","01D","Other")) 

ggplot(df2,aes(x=filen,y=ge_freq,fill=mut_lab)) + geom_bar(position = 'fill',stat='identity') +
  scale_y_continuous(labels = percent, limits = c(0,1), expand = c(0,0)) +
  xlab("Source") + 
  ylab("Mutation frequency") +
  labs(title = paste0("A Mutations"), fill = legend_title) +
  theme(axis.text.x = element_text(face="bold", color="#993333", 
                                   size=14, angle=90,hjust=0)) +
  facet_grid(~group,scales="free",space="free") + 
  scale_fill_manual(values=c("steelblue4","springgreen3","mediumpurple1","lightgoldenrod1",
                             "magenta","aquamarine","red4","palevioletred3", "yellow2",
                             "olivedrab3","lightskyblue","orange","lightpink","snow4"),
                    breaks=c("13D","12D","11D","10D","09D","08D","07D","06D","05D","04D",
                             "03D","02D","01D","Other")) 
dev.off()

