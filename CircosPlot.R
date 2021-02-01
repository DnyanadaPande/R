#
# Dnyanada Pande, October 2020
# 
# This script reads the chromosomal integration site positions from csv files and
# is used to make a circos template file to generate a histogram to chromosomal
# integration positions across the rhesus macaque/human/canine genome.

# Setting the file paths for the directories.
circos_dir <- '/Volumes/Dnyanada_HD/RIS/Circos/'
setwd(circos_dir)
test_dir <- '/Volumes/Dnyanada_HD/RIS/Circos/'

files <- list.files(path=circos_dir, pattern="File", full.names=TRUE, recursive=FALSE)
files <- basename(files)
files
tests = data.frame(files,stringsAsFactors = FALSE)
tests
colnames(tests) = c("file")

# Reading the input files and adding the chromosome end position.
lt = lapply(c(1:length(tests$file)),function(t) {
  
  fn = tests[t,1]
  print(fn)
  df = read.csv(fn, header = T, stringsAsFactors = F, comment.char = "#")
  print(head(df)) 
  colnames(df) <- c("chr","start","strand")
  head(df)
  
  df <- df[-c(3)]
  # Removing the entries that do not align with a specific chromosome.
  df <- df[!(df$chr == "chrUn" | df$start == "chrUn"),]
  
  #Calculating the end positions based on the start position of the integration site.
  end <- numeric()
  for(j in 1:length(df$chr)){
    ep <- as.numeric(df$start[j]) + 3
    end <- c(end,ep)
  }
  
  df$end <- end
  head(df)
  write.csv(df,paste0("Genes_",fn),row.names = F)
  
})

# Combine all the files in one dataframe for the plot.
files <- list.files(path=circos_dir, pattern="Genes_", full.names=TRUE, recursive=FALSE) 
files <- basename(files)
files
tests = data.frame(files,stringsAsFactors = FALSE)
tests
colnames(tests) = c("file")

df2 <- data.frame(matrix(nrow = 0, ncol = 3))
colnames(df2) <- c('chr','start','end')

for (j in 1:length(tests$file)){
  
  d <- read.csv(tests$file[j], header = T, stringsAsFactors = F)
  print(head(d))
  df2 <- rbind.data.frame(df2,d)
}

#Selecting the unique chromosomal integration positions.
df2.uniq <- unique.data.frame(df2)

write.csv(df2.uniq,paste0(circos_dir,'A_pb_wbc_integrations.csv'),row.names = F)


# Read the circos input file and the total insertion sites for Rhesus macaque.
# Make a dataframe, update the chromosomal bins & increment the counts.
# Circos rhesus macaque reference file format:
# chr	start	  end	    ct
# rm1	0	      1000000	
# rm1	1000001	2000000	

rh <- read.csv(paste0(test_dir,"Circos_template.csv"), header = T, stringsAsFactors = F)
rh$ct <- rep(0,length(rh$ct))
head(rh)

df3 <- read.csv(paste0(circos_dir,"A_pb_wbc_integrations.csv"), header = T, stringsAsFactors = F)
head(df3)

# Increment the counts by determining in which bin the integration site falls.
for (i in 1:length(df3$start)){
  for (j in 1:length(rh$chr)){
    if (df3$chr[i] == rh$chr[j]){
      print(rh$chr[j])
      
      if (df3$start[i] >= rh$start[j] & df3$start[i] <= rh$end[j]){
        rh$ct[j] <- rh$ct[j] + 1
        print(rh$ct[j])
      }
    }
  }
}

write.csv(rh, "A_total_PB_WBC_circos.csv", row.names = F)




