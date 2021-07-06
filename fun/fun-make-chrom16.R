# extract chromatogram from skyline file using skyline runner


make.chrom <- function(skyfile, curr.work.dir, chrom.type = "PRM"){


  chrom.file   = paste(curr.work.dir, gsub(".sky", ".chrom", basename(skyfile)), sep="")
  
  if(chrom.type == "MS2" | is.na(chrom.type)){
    command = paste(runner, 
                    paste("--in",skyfile, sep="="),
                    paste("--chromatogram-file",chrom.file, sep="="),
                    "--chromatogram-products")
  } else if (chrom.type == "MS1"){
    command = paste(runner, 
                    paste("--in",skyfile, sep="="),
                    paste("--chromatogram-file",chrom.file, sep="="),
                    "--chromatogram-precursors")
  }
  cat(command)
  system(gsub("/", "\\\\", command) )
  
  if(file.exists(chrom.file)){
    return(chrom.file)
  } else {
    return(NA)
  }
}
