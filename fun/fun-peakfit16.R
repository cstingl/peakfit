# redo pt classification

library(tidyverse)
library(minpack.lm)
library(broom)
library(purrr)


peak.ntop50.min = 2

count_above <- function(x, t = 0){ sum( x > t) }

integrate_peak <- function(peak.data, int.start, int.stop, npick, nquan){
  data.A <- peak.data %>% filter(rt >= int.start & rt <= int.stop) %>% select(rt, int)
  return(  sum(diff(data.A$rt) * (head(data.A$int,-1) +  tail(data.A$int,-1))/2)  )  # reference: https://stackoverflow.com/questions/4954507/calculate-the-area-under-a-curve
  
}

rt_at_height_logN <- function(curr_logN.mod, qry_height, start_rt, stop_rt, n_bin, match_accuracy){  # v16
  n_iter = 0    
  repeat{
    test.x <- seq(start_rt, stop_rt, length.out = n_bin)
    test.y <- curr_logN.mod(logN.x = test.x)
    test.y[is.na(test.y)] <- 0

    bin = which(abs(test.y-qry_height)==min(abs(test.y-qry_height)))

    ref_err = abs((test.y[bin] - qry_height)/qry_height)
    
    if(bin > 1)     {  start_rt = test.x[bin-1] }
    if(bin < n_bin) {  stop_rt  = test.x[bin+1] }
    
    
    if(ref_err < match_accuracy | n_iter > 1000){
      h  = test.y[bin]
      rt = test.x[bin]
      break
    } 
    
    n_iter = n_iter + 1
  }
  return(tibble(rt=rt, height=h, height_deviation = ref_err, iteration = n_iter))
}



# peak.analysis.logN <- possibly( function(m.data, rt, height, width, sym){
# 
#     model <- try(   nlsLM( y ~ m.height * exp(- (log(2)/((log(m.sym))^2) )  *   (log(((x - m.rt)/m.width )   *   (m.sym - 1)/m.sym   +   1 ))^2 # ),                                                data = m.data, 
#                                                 start = list(m.rt = peak.rt, m.width = peak.width, m.height = peak.height , m.sym = peak.sym))
#                                 , silent = TRUE)
# 
# }, NA)

#                   _                       _           _                  __
#  _ __   ___  __ _| | __  __ _ _ __   __ _| |_   _ ___(_)___   ___  __ _ / _| ___
# | '_ \ / _ \/ _` | |/ / / _` | '_ \ / _` | | | | / __| / __| / __|/ _` | |_ / _ \
# | |_) |  __/ (_| |   < | (_| | | | | (_| | | |_| \__ \ \__ \_\__ \ (_| |  _|  __/
# | .__/ \___|\__,_|_|\_(_)__,_|_| |_|\__,_|_|\__, |___/_|___(_)___/\__,_|_|  \___|
# |_|                                         |___/

peak.analysis.safe <- possibly( function(.x, peak.model){
  peak.model = "logN"
  
  peak <- .x %>% 
    select(Times, Intensities) %>% separate_rows(Times, Intensities, sep = ",") %>% 
    mutate(Times = as.numeric(Times), Intensities = as.numeric(Intensities)) %>%
    group_by(rt = Times) %>% 
    summarize(int_sum = sum(Intensities), int_prod = prod(Intensities), n_traces = count_above(Intensities, 0), .groups = "drop_last" ) %>%
    mutate(int.pick = if_else(n_traces >= npick, int_sum, 0),
           int = if_else(n_traces >= nquan, int_sum, 0))
  
#  colnames(peak) <- c("rt","int")
#  peak <- peak %>% select(rt, int) %>% arrange(rt)
  
  ptm <- proc.time()
  {  # determine peak characteristics from raw data
      peak.height = peak %>% summarize(max = max(int.pick), .groups = "drop_last") %>% pull(max)
      peak <- peak %>% mutate(apex  = if_else(int.pick == peak.height,1,0))
      
      peak.rt = peak %>% filter(apex == 1) %>% top_n(1, rt) %>% pull(rt)
      
      peak <- peak %>% mutate(top50 = if_else(int > peak.height * 0.50 & rt > peak.rt - maxwidth/2 & rt < peak.rt + maxwidth/2, 1, 0))
      peak <- peak %>% mutate(use   = if_else(int > peak.height * 0.01 & rt > peak.rt - maxwidth/2 & rt < peak.rt + maxwidth/2, 1, 0))
      
      
      peak.ntop50 = peak %>% filter(top50 == 1) %>% summarize(n = length(int), .groups = "drop_last") %>% pull(n)
      
      # last data point below 50% height
      front.low   = peak %>% filter(rt <= peak.rt & top50 == 0) %>% select(rt, int) %>% top_n(1, rt)
      front.high  = peak %>% filter(rt <= peak.rt & top50 == 1) %>% select(rt, int) %>% top_n(1, -rt)
      tail.high   = peak %>% filter(rt >= peak.rt & top50 == 1) %>% select(rt, int) %>% top_n(1, rt)
      tail.low    = peak %>% filter(rt >= peak.rt & top50 == 0) %>% select(rt, int) %>% top_n(1, -rt)
      data.tstart = peak %>% summarize(min = min(rt), .groups = "drop_last") %>% pull(min)
      data.tend   = peak %>% summarize(max = max(rt), .groups = "drop_last") %>% pull(max)
      
      
      # RT at 50 %
      
      front.k    = (front.high %>% pull(int)  - front.low %>% pull(int) )/ (front.high %>% pull(rt)  - front.low %>% pull(rt) )
      front.d    = front.low %>% pull(int) - (front.low %>% pull(rt)) * front.k
      front.rt50 = (peak.height * 0.5 - front.d)/front.k 
    
      tail.k     = (tail.high %>% pull(int)  - tail.low %>% pull(int) )/ (tail.high %>% pull(rt)  - tail.low %>% pull(rt) )
      tail.d     = tail.low %>% pull(int) - (tail.low %>% pull(rt)) * tail.k
      tail.rt50  = (peak.height * 0.5 - tail.d)/tail.k 
      
      if(length(front.rt50) == 1 & length(tail.rt50) == 1){
        peak.width = tail.rt50 - front.rt50  
        peak.sym   = (peak.rt - front.rt50) / (tail.rt50 - peak.rt)
        peak.log = "front/tail:ok"
      } else if (length(front.rt50) == 1){  # estimating peak with and symetry
        peak.width = 2 * (peak.rt -  front.rt50 )
        peak.sym   = 1.05
        peak.log = "incomplete peak data: tail missing, width and symetry estimated."
      } else if (length(tail.rt50) == 1){
        peak.width = 2 * (tail.rt50 - peak.rt )
        peak.sym   = 1.05
        peak.log = "incomplete peak data: front missing, width and symetry estimated."
      } else {
        peak.width = NA
        peak.sym   = NA
        peak.log = "incomplete peak data: width and symetry cannot be determined."
      }
      
      if(!is.na(peak.sym)){
        peak.margins.raw  =  c(peak.rt - 2*peak.width/peak.sym, peak.rt + 2*peak.width*peak.sym)
        A_raw  = integrate_peak(peak, peak.margins.raw[1], peak.margins.raw[2])
      } else {
        peak.margins.raw  =  c(NA, NA)
        A_raw = NA
      }
    
  }
  
  if(peak.height > 0 & peak.ntop50 >= peak.ntop50.min & !is.na(peak.width) ){
  
      m.data <- as.data.frame(  peak %>% filter(use == 1) %>% select(x = rt, y = int) %>% filter(x > peak.rt - 3*peak.width & x < peak.rt + 3*peak.width) )
    
      if(peak.sym <= 1){
        peak.sym.start = 1.01
      } else {
        peak.sym.start = peak.sym
      }
      
      
      #    __                   __                           _
      #   / /  ___   __ _    /\ \ \___  _ __ _ __ ___   __ _| |
      #  / /  / _ \ / _` |  /  \/ / _ \| '__| '_ ` _ \ / _` | |
      # / /__| (_) | (_| | / /\  / (_) | |  | | | | | | (_| | |
      # \____/\___/ \__, | \_\ \/ \___/|_|  |_| |_| |_|\__,_|_|
      #             |___/

      model.logN <- try(   nlsLM( y ~ m.height * exp(- (log(2)/((log(m.sym))^2) )  *   (log(((x - m.rt)/m.width )   *   (m.sym^2 - 1)/m.sym   +   1 ))^2 ),                                                
                                  data = m.data, 
                             #start = list(m.rt = peak.rt, m.width = peak.width, m.height = peak.height , m.sym = peak.sym.start))
                             start = list(m.rt = peak.rt, m.width = peak.width, m.height = peak.height , m.sym = 1.05))
                      , silent = TRUE)

      

      
      if(class(model.logN) == "nls") { #!is.null(model.logN) | is.na(model.logN) ){
        
        peak.height.mod = tidy(model.logN) %>% filter(term == "m.height") %>% pull(estimate)
        peak.rt.mod     = tidy(model.logN) %>% filter(term == "m.rt") %>% pull(estimate)
        peak.width.mod	= tidy(model.logN) %>% filter(term == "m.width") %>% pull(estimate)
        peak.sym.mod    = tidy(model.logN) %>% filter(term == "m.sym") %>% pull(estimate)

        logN.HatW = 2
        logN.mod  <- function(logN.x){  peak.height.mod * exp(- (log(logN.HatW)/((log(peak.sym.mod))^2) )  *   (log(((logN.x - peak.rt.mod)/(peak.width.mod) )   *   (peak.sym.mod^2 - 1)/peak.sym.mod   +   1 ))^2 )  }  
        
        # calculate peak boundaries at 1% of peak height

        if(pbtype == "width"){
          peak.margins.logN =  c(peak.rt.mod - pbfactor/2*peak.width.mod/peak.sym.mod, peak.rt.mod + pbfactor/2*peak.width.mod*peak.sym.mod)
          # print(str_c("margins at PBfactor = ", pbfactor, " -> ", peak.margins.logN[1], "..", peak.margins.logN[2]))
        } else {
          curr_qry_height = peak.height.mod * pbfactor
          peak.margins.logN <- c( rt_at_height_logN(curr_logN.mod = logN.mod, 
                                                    qry_height = curr_qry_height, 
                                           start_rt  = peak.rt.mod-(2*peak.width.mod/peak.sym.mod), 
                                           stop_rt   = peak.rt.mod, 
                                           n_bin     = 20,  
                                           match_accuracy = 0.01) %>% pull(rt),
                                  rt_at_height_logN(curr_logN.mod = logN.mod, 
                                                    qry_height = curr_qry_height, 
                                           start_rt  = peak.rt.mod, 
                                           stop_rt   = peak.rt.mod+(2*peak.width.mod*peak.sym.mod), 
                                           n_bin     = 20,  
                                           match_accuracy = 0.01) %>% pull(rt))
          # print(str_c("margins at PBfactor = ", pbfactor, " -> ", peak.margins.logN[1], "..", peak.margins.logN[2]))

        }
        A_logN = integrate(logN.mod, peak.margins.logN[1], peak.margins.logN[2])$value

        peak.result <- tibble(t0     = data.tstart,
                              tend   = data.tend,
                              rt     = peak.rt,
                              height = peak.height,
                              w50    = peak.width,
                              sym50  = peak.sym,
                              start  = peak.margins.raw[1],
                              stop   = peak.margins.raw[2],
                              A      = A_raw,
                              #log    = peak.log,
                              height_logN = peak.height.mod,
                              rt_logN     = peak.rt.mod,
                              width_logN	= peak.width.mod,
                              sym_logN    = peak.sym.mod,
                              start_logN  = peak.margins.logN[1],
                              stop_logN   = peak.margins.logN[2],
                              A_logN      = A_logN
                              )
      } else { # if no model could be generated
        peak.result <- tibble(t0     = data.tstart,
                              tend   = data.tend,
                              rt     = peak.rt,
                              height = peak.height,
                              w50    = peak.width,
                              sym50  = peak.sym,
                              start  = peak.margins.raw[1],
                              stop   = peak.margins.raw[2],
                              A      = A_raw,
                              #log    = peak.log,
                              height_logN = NA,
                              rt_logN     = NA,
                              width_logN	= NA,
                              sym_logN    = NA,
                              start_logN  = NA,
                              stop_logN   = NA,
                              A_logN      = NA
        )
        
      }

    
  } else {
      peak.result <- tibble(t0     = NA,
                            tend   = NA,
                            rt     = peak.rt,
                            height = peak.height,
                            w50    = NA,
                            sym50  = NA,
                            start  = NA,
                            stop   = NA,
                            A      = 0,
                            height_logN = NA,
                            rt_logN     = NA,
                            width_logN	= NA,
                            sym_logN    = NA,
                            start_logN  = NA,
                            stop_logN   = NA,
                            A_logN      = NA
                            )
  }
  return(peak.result)
  
} ,

  peak.result <- tibble(t0     = NA,
                        tend   = NA,
                        rt     = NA,
                        height = NA,
                        w50    = NA,
                        sym50  = NA,
                        start  = NA,
                        stop   = NA,
                        A      = NA,
                        #log    = peak.log,
                        height_logN = NA,
                        rt_logN     = NA,
                        width_logN	= NA,
                        sym_logN    = NA,
                        start_logN  = NA,
                        stop_logN   = NA,
                        A_logN      = NA
  )
)

