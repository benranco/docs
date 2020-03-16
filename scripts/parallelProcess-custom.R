args <- commandArgs(trailingOnly = TRUE)
path <- args[1]
file <- "report_p35.Rds" #args[2]
haploidOrDiploid <- args[3]

#path <- "/home/gosuzombie/Desktop/region_38/run"
#file <- "report_p432.Rds"

#message(paste0("looking up NA data on ", file))
options(stringsAsFactors = FALSE, warn = 1)
message("running parallel process CUSTOM")

report <- readRDS(paste0(path, "/reporttemp/", file))

if(!("COMBINED" %in% colnames(report)))
{
  s <- 4
}else 
{
  s <- 5
}

increment <- 1
#len <- table(is.na(report))["TRUE"][[1]]

for(a in 1:nrow(report))
{
  if(s <= ncol(report))
  {
    for(b in s:ncol(report))
    {
      if(is.na(report[a,b]))
      {
        fn <- paste(substr(colnames(report)[b], 0 , nchar(colnames(report)[b]) - 4), 
                    "_sorted.bam", sep = "")
        cmd <- paste0(path, "/tools/samtools-1.3.1/samtools tview ", path ,"/dataTemp/single/", fn ,
                      " ", path, "/reference/formatted_output.fasta -d T", 
                      ' -p \"', paste(report[a, "CHROM"], report[a, "POS"], sep = ":"), '"')
        out <- as.matrix(read.delim(pipe(cmd), sep = "\n"))
        
        if(substring(out[1,1], 1 ,1) == report[a, "REF"])
        {
          if(nrow(out) >= 2)
          {
            if(substring(out[2,1], 1 ,1) == "." || substring(out[2,1], 1 ,1) == ",")
            {
              # 1 == haploid, 2 == diploid. If it's haploid, we follow the format in the .tab file of "A/",
              # whereas if it's diploid we follow the format in the .tab file of "A/A".
              if (haploidOrDiploid == 1) {
                report[a,b] <- paste0(report[a, "REF"], "/")
              } else {
                report[a,b] <- paste0(report[a, "REF"], "/", report[a, "REF"])
              }
            }
          }
        }
        #message(increment / len)
        increment <- increment + 1
      }
    }
  }
}

saveRDS(report, file = paste0(path, "/reporttemp/", substr(file, 0, nchar(file) - 4), "_filled.Rds"))
