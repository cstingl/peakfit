{ # PARAMETERS ----
  library(tidyverse)
  
  # curr.skyfile        = "T:/testarea/PeakFit.rproj/data/peakfit-test-data.sky"
  # runner              = "T:/testarea/skyline/runner/SkylineDailyRunner.exe"
  # curr.peakqc1.file   = "T:/testarea/PeakFit.rproj/data/peakfit-test-data.peakqc.txt"
  # curr.work.dir       = str_c(dirname(skyfile), "/")
}

# export peak boundaries report (skyrunner) ----
make.pbrep <- function(curr.skyfile, curr.work.dir, runner){
  
  pbrep.file   = paste(curr.work.dir, gsub(".sky", ".pbrep1.csv", basename(curr.skyfile)), sep="")
  
  command = paste(runner, 
                  paste("--in",curr.skyfile, sep="="),
                  paste("--report-file", pbrep.file, sep="="),
                  '--report-name="Peak Boundaries"')
#  return(   gsub("/", "\\\\", command)    )
  system(gsub("/", "\\\\", command) )
  
  if(file.exists(pbrep.file)){
    return(pbrep.file)
    #return(command)
  } else {
    #return(NA)
    return(command)
  }

}


update_pbrep <- function(curr.peakqc1.file, curr_xsigma, curr.skyfile, curr.work.dir, curr.label.p1, runner){ # update peak boundaries ----
  
  pbrep.file.1 = make.pbrep(curr.skyfile, curr.work.dir,  runner)

  if(!is.na( pbrep.file.1 )){
    pbrep <- read_csv(pbrep.file.1) %>%
      left_join(read_tsv(curr.peakqc1.file) %>% 
                  filter(IsotopeLabelType == curr.label.p1) %>%
                  # mutate(t_start = round( rt_logN - curr_xsigma/2 * width_logN / sym_logN, 3),
                  #        t_stop  = round( rt_logN + curr_xsigma/2 * width_logN * sym_logN, 3) ) %>%
                  mutate(t_start = round( start_logN, 3),
                         t_stop  = round( stop_logN, 3) ) %>%
                  mutate(t_start  = ifelse(t_start < t0, t0, t_start),
                         t_stop  = ifelse(t_start > tend, tend, t_stop)) %>%
                  mutate(t_start  = ifelse(t_start >= t_stop, NA, t_start),
                         t_stop   = ifelse(t_start >= t_stop, NA, t_stop)) %>%
                  select(`File Name` = FileName, 
                         `Peptide Modified Sequence` = PeptideModifiedSequence, 
                         rt_logN, width_logN/sym_logN, t_start, t_stop) ,
                by = c("File Name", "Peptide Modified Sequence")) %>%
      mutate(`Min Start Time` = ifelse(is.na(t_start) | is.na(t_stop), `Min Start Time`, t_start),
             `Max End Time`   = ifelse(is.na(t_start) | is.na(t_stop) , `Max End Time`, t_stop)) %>%
      select(-rt_logN, -width_logN, -t_start, -t_stop) %>%
      write_csv(gsub(".pbrep1.csv", ".pbrep2.csv", pbrep.file.1))
    
    check_overlap <- read_csv(pbrep.file.1) %>% 
      inner_join(read_tsv(curr.peakqc1.file) %>% select(`File Name` = FileName, 
                                                     `Peptide Modified Sequence` = PeptideModifiedSequence) ,
                 by = c("File Name", "Peptide Modified Sequence")) %>% pull(`File Name`) %>% length
    return(str_c(pbrep.file.1, " --> ", curr.peakqc1.file, "(n = ",check_overlap,")"))
  } else {
    return(-1)
  }
  
}




