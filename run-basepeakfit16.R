library(tidyverse)
library(ggrepel)

source("fun/fun-peakfit16.R")

as_num <- function(x){ return(as.numeric(as.character(x))) }

rel_int_thresh = 0.05
maxwidth = 0.5


args = commandArgs(trailingOnly=TRUE)  

#args = c("T:/a/a166/submission3/a166-example3/210421OFc1_PM38_025.tsv")
args <- gsub("\\\\", "/", args, perl=TRUE)

argerror=0

if(exists("args") ) {
  
  bpfile = args[1]
  
  if(file.exists(bpfile)){
    cat(paste(" - base peak file", bpfile, " ... exists\n"))
  } else {
    cat(paste(" - base peak file", bpfile, " DOES NOT EXIST!\n"))
    argerror=1
  }
}
  
if(argerror > 0) {
  cat("Usage: run-peakfit*.R raw=rawline-file.raw\n")
  cat("-- raw      ... rawline file\n")

} else {
  
  # out file name
  peakqc1.file = gsub(".tsv", ".bpqc.txt", bpfile)
  
  # load base peak data
  data <- read_tsv(bpfile,
                   col_types = "ldcdccddddddddddd") %>%
    select(-`NULL`) %>%
    filter(msLevel == "ms1") 


  int_thresh = max(data$basePeakInt) * rel_int_thresh

  data <- data %>%
    mutate(use = if_else(basePeakInt >= !!int_thresh,1,0),
           use = as.factor(use))
  
  # scan for peak data 
  {
    index <- data$index
    bpcmz <- data$basePeakMZ
    peak  = 1
    peaks <- c(peak)    
    
    for(i in 2:length(index)){
        if(abs(data$basePeakMZ[i] - data$basePeakMZ[i-1] ) > 0.1){ peak = peak + 1 }
        peaks <- c(peaks, peak)
    }
    data$peak <- peaks
  }

  { # preparing peak data table
    # peak ids are defined by m/z and RT
    peakids <- data %>%
      filter(use == 1) %>%
      group_by(peak) %>%
      mutate(rank_int     = rank(-basePeakInt) ) %>%
      filter(rank_int == 1) %>%
      mutate(peakid = str_c("pk", round(basePeakMZ,1),"m/z @ ", round(rt/60,2),"min", sep = "")) %>%
      select(peak, peak_mz = basePeakMZ, peak_rt = rt, peakid)

    peak_data <- data %>%
      mutate(rt = round(rt/60,3)) %>%
      select(peak, rt, basePeakInt,  index) %>%
      left_join(peakids, by = "peak") %>%
      filter(!is.na(peakid)) %>%
      group_by(peak, peakid) %>%
      summarize(n           = n_distinct(index),
                Times       = paste0(rt, collapse=","),
                Intensities = paste0(basePeakInt, collapse="," ),
                .groups = "drop")
  }
  
  { # run peak fit
    # check if npick, nquan, and maxwidth are defined
    npick=0
    nquan=0
    
    start.time <- proc.time()
    R <- peak_data %>% 
      mutate(npick = 0, nquan = 0) %>%
      group_by(peak, peakid, n) %>% 
      nest() %>% 
      mutate(., data = map(data, peak.analysis.safe)) %>% 
      unnest(cols = c(data))  
    run.time = proc.time() - start.time
  }



  { # exporting peak basepeak QC file:
  
    R2 <- R %>% 
      mutate(t0 = round(t0, 2), tend = round(tend, 2), rt = round(rt, 2), height = round(height, 2)) %>%
      mutate(w50 = round(w50, 3), sym50 = round(sym50, 3)) %>% 
      mutate(a = w50 - w50/(sym50+1), b = w50 - a  ) %>%
      mutate(dfront = round((rt - t0)/a,2), 
             dtail = round((tend - rt)/b,2),
             start_excess = round(start - t0, 3),
             stop_excess  = round(tend - stop, 3),
             flag = ifelse(stop_excess < 0, "TT", ifelse(start_excess < 0, "FT",""))
      ) %>%
      mutate(height_logN = round(height_logN, 3), rt_logN = round(rt_logN, 3), width_logN = round(width_logN, 4), sym_logN = round(sym_logN, 4)) %>%   
      select(peak, peakid, n, t0, tend, rt, height, w50, sym50, start, stop, A, dfront, dtail, height_logN, rt_logN, width_logN, sym_logN, start_logN, stop_logN, A_logN, start_excess, stop_excess )
    
    R2 %>%
      write_tsv(peakqc1.file)
    
    R2 %>%
      filter(!is.na(rt_logN)) %>%
      arrange(rt_logN) %>%
      select(peakid, rt_logN, height_logN, width_logN, sym_logN, A_logN)
    
  }
}
