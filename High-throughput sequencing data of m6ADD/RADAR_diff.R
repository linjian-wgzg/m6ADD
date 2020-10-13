if (!requireNamespace("BiocManager", quietly = TRUE))
install.packages("BiocManager")
BiocManager::install("devtools")
library(devtools)
install_github("scottzijiezhang/RADAR")
library(RADAR)
install.packages("Rtools")
library(devtools)
install_github("scottzijiezhang/RADAR")
library("RADAR")
radar <- countReads(
  
  samplenames = c("ctrl1","ctrl2","ctrl3","ctrl4","ctrl5","ctrl6","ctrl7","case1","case2","case3","case4","case5","case6"),
  
  gtf = "/Note file path/gencode.v32.annotation.gtf",
  
  bamFolder = "./bam file path",
  
  modification = "m6A",
  
  outputDir = "/Output path",
  
  threads = 10
  
)

radar <- normalizeLibrary( radar )

radar <- adjustExprLevel( radar )

variable(radar) <- data.frame( Group = c("Ctl","Ctl","Ctl","Ctl","Ctl","Ctl","Ctl","Case"£¬"Case"£¬"Case"£¬,"Case","Case","Case") )

radar <- filterBins( radar ,minCountsCutOff = 15)

radar <- diffIP_parallel(radar,thread = 10)

radar <- reportResult( radar, cutoff = 0.1, Beta_cutoff = 0.5 )

res <- results( radar )

write.table(res,file='/Output path/RADAR_diff_site.xls',sep='\t',quote=F)

