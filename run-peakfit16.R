#                 _    __ _ _
#  _ __  ___ __ _| |__/ _(_) |_
# | '_ \/ -_) _` | / /  _| |  _|
# | .__/\___\__,_|_\_\_| |_|\__|
# |_|

# peakfit, version: 210428, author: Christoph Stingl


library(tidyverse)
source("fun/fun-make-chrom16.R")
source("fun/fun-peakfit16.R")
source("fun/fun-gen-pbrep16.R")

{ # PARAMETERS ----
  # Location of skylinerunner
  runner  = "C:/Tools/pwiz/SkylineDailyRunner.exe"  # location of skyline(daily)runner
  
  # Default parameters
  default.label.p1 = "light"
  default.npick = 3
  default.nquan = 1
  default.maxwidth = 0.5
  default.mslevel = "MS2"
  default.pbtype = "relheight"
  default.pbfactor = 0.01
}
# end parameters.

args = commandArgs(trailingOnly=TRUE)  

argerror=0

if(exists("args") ) {
  args <- gsub("\\\\", "/", args, perl=TRUE)
  cat("# parameters:\n")
  for (arg in args){
    if(grepl("=", arg)){
      pv <- strsplit(arg, split="=", fixed=TRUE)[[1]]
      para = tolower(pv[1])
      if(para == "sky"){
        skyfile = pv[2]

        if(file.exists(skyfile)){
          cat(paste(" - skyfile", skyfile, " exists\n"))
        } else {
          cat(paste(" - skyfile", skyfile, " DOES NOT EXIST!\n"))
          argerror=1
        }
      
      } else if (para == "workdir"){
        
        work.dir = pv[2]
        
        if(dir.exists(work.dir)){
          cat(paste(" - workdir:", work.dir))
        } else {
          cat(paste(" - work dir", work.dir, " DOES NOT EXIST!\n"))
          argerror=1
        }
        
      } else if (para == "label"){
        label.p1 = tolower(pv[2])
        
        if(label.p1 == "light" | label.p1 == "heavy"){
          cat(paste(" - label used for 1st pass: ", label.p1,"\n"))
        } else {
          cat(paste(" - undefined label", label.p1, "! Use light or heavy.\n"))
          argerror=1
        }
      
      } else if (para == "npick") {
        npick = pv[2]
      } else if (para == "nquan") {
        nquan = pv[2]
      } else if (para == "maxw"){
        maxwidth = pv[2] 
      } else if (para == "mslevel") {
        mslevel = pv[2] 
      } else if (para == "pbtype") {
        pbtype = pv[2] 
      } else if (para == "pbfactor") {
        pbfactor = as.numeric(as.character(pv[2]))
      } else {
        cat(paste("!! UNKNOWN PARAMETER",para, "(", pv[2],")!\n"))
      }
        
    } else {
      cat(paste("!Undefined parameter: ", para,"\n"))
      argerror=1
    }
    
    cat(paste("--",arg,"\n"))
  }
} else {
  cat("No parameters specifed!\n")
  argerror=1
}
  
if(argerror == 0){
  if(!exists("work.dir")){
    work.dir = dirname(skyfile)
  }
  
  if(!grepl("/$", work.dir)){
    work.dir = paste(work.dir,"/", sep="")
  }
  
  if(!exists("label.p1")){
    label.p1 = default.label.p1
  }
  if(!exists("npick")){
    npick = default.npick
  }
  if(!exists("nquan")){
    nquan = default.nquan
  }
  if(!exists("maxwidth")){
    maxwidth = default.maxwidth
  }  
  if(!exists("mslevel")){
    mslevel = default.mslevel
  }
  if(!exists("pbtype")){
    pbtype = default.pbtype
  }
  if(!exists("pbfactor")){
    pbfactor = default.pbfactor
  }
  
  
  cat(paste("-- sky     :", skyfile,"\n"))
  cat(paste("-- dir     :", work.dir,"\n"))
  cat(paste("-- label   :", label.p1,"\n"))
  cat(paste("-- npick   :", npick,"\n"))
  cat(paste("-- nquan   :", nquan,"\n"))
  cat(paste("-- maxw    :", maxwidth,"\n"))
  cat(paste("-- mslevel :", mslevel,"\n"))
  cat(paste("-- pbtype  :", pbtype,"\n"))
  cat(paste("-- pbfactor:", pbfactor,"\n"))
} else{
  cat("Usage: run-peakfit*.R sky=skyline-file.sky dir=d:/workdir/ label=heavy\n")
  cat("-- sky      ... skyline file\n")
  cat("-- dir      ... working directory (storage of peakqc and chrom files)\n")
  cat("-- label    ... filter of label in first pass; if not specified, all peaks will be processed (not 2nd pass)\n")
  cat("-- nquan    ... traces required (i > 0) for peak picking (default = 3)\n")
  cat("-- nquan    ... traces required (i > 0) for peak intergration (default = 1)\n")
  cat("-- maxw     ... maximum peak width at baseline, in minutes (default = 0.5)\n")
  cat("-- mslevel  ... MS level of data: MS1 or MS2 (default = MS2)\n")
}
  
if(argerror == 0){
  exp.chrom.file   = paste(work.dir, gsub(".sky", ".chrom", basename(skyfile)), sep="")
  
  if(file.exists(exp.chrom.file)){
    cat(paste("--> ", exp.chrom.file, " exists (skipping generation)."))
    chrom.file = exp.chrom.file
  } else {
    cat(paste("--> Chromatograms '", exp.chrom.file, "' does not exists (expected).", sep = ""))
    chrom.file <- make.chrom(skyfile, work.dir, mslevel) 
    cat(paste("--> Chromatograms written to '", chrom.file,"'", sep = ""))
  }
} else {
  cat("An error accured!")
}
  
  

      
if( is.na(chrom.file) ){
  print("Generation of chromatography file (chrom file) failed!\n")
} else {
  peakqc1.file = gsub("\\.chrom", ".peakqc1.txt", chrom.file)

  chrom <- read_tsv(chrom.file, col_types = "ccddcdcdcc")
  
  start.time <- proc.time()
  R <- chrom %>% 
    mutate(npick = npick, nquan = nquan) %>%
    select(FileName,PeptideModifiedSequence,PrecursorCharge,IsotopeLabelType,Times,Intensities) %>%
    group_by(FileName,PeptideModifiedSequence,PrecursorCharge,IsotopeLabelType) %>% 
    nest() %>% 
    mutate(., data = map(data, peak.analysis.safe)) %>% 
    unnest(cols = c(data))  
  run.time = proc.time() - start.time

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
    select(FileName, PeptideModifiedSequence, PrecursorCharge, IsotopeLabelType, t0, tend, rt, height, w50, sym50, start, stop, A, dfront, dtail, height_logN, rt_logN, width_logN, sym_logN, start_logN, stop_logN, A_logN, start_excess, stop_excess, flag )
  
  R2 %>%
    write_tsv(peakqc1.file)
  
  print(paste(length(R2$FileName), "peaks in  ", round(run.time[[3]],1), " seconds -> ", round(run.time[[3]]*1000/length(R$FileName), 1)," ms per peak.") )
  print(paste("Result stored in", peakqc1.file) )
}

if(file.exists(peakqc1.file)){
  print(paste("-- ", peakqc1.file," generated.", sep = "") )
  
  update_pbrep(peakqc1.file, xsigma, skyfile, work.dir, label.p1, runner)
  
  cat("--> Generate peak boudaries report.\n")
  
} else {
  print(paste("First pass failed! No '", peakqc1.file,"'!", sep = "") )
}


cat("--> finished.\n")

