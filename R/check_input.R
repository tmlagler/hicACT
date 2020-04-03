check_input <- function(indata, kb, h, pthres, bin1, bin2, ccount, pval, ignore_warnings){

  # check if bin IDs are positive integers
  if(any(indata[,eval(bin1)]%%1!=0)>0 | any(indata[,eval(bin1)]<=0)){
    stop("Fragment identifiers must be positive integers")
  }
  if(any(indata[,eval(bin2)]%%1!=0)>0 | any(indata[,eval(bin2)]<=0)){
    stop("Fragment identifiers must be positive integers")
  }

  # check if contact counts are integers >=0
  if(any(indata[,eval(ccount)]%%1!=0)>0 | any(indata[,eval(ccount)]<0)){
    stop("Observed contact counts must be integers >=0")
  }

  # check if raw p-values are in range
  if(any(indata[,eval(pval)]<0) | any(indata[,eval(pval)]>1)){
    stop("p-values must be in [0, 1]")
  }

  # check if p-value threshold is reasonable
  if(pthres>=1){
    warning("HiC-ACT applied to all p-values")
  }else if(pthres<=1.0e-10){
    warning("Consider appling a less stringent p-value threshold")
  }else if(pthres>.3 & pthres<1){
    warning("Consider appling a more stringent p-value threshold")
  }

  # check if any missing values in data
  if(any(is.na(indata))){
    stop("Remove rows with missing values")
  }

  # check if Kb resolution is integer and in reasonable range
  if(kb%%1!=0 | kb<=0 | kb > 200 ){
    stop("Kb resolution should be positive integer in (0, 200]")
  }

  # if h is specified, check if positive integer
  # also check if reasonable value
  if(is.null(h)==F){
    if(h%%1!=0 | h<=0){
      stop("h must be a postive integer")
    }
    if(kb %in% 20:25 & (h >= 20 | h <= 5)){
      if(ignore_warnings==F){
        stop("Smoothing window (h) may not be of appropriate size for specified data resolution. To proceed, specify 'ignore_warnings=T' as function parameter")
      }else{warning("Smoothing window (h) may not be of appropriate size for specified data resolution")}
    }
    if(kb %in% 16:19 & (h >= 30 | h <= 10)){
      if(ignore_warnings==F){
        stop("Smoothing window (h) may not be of appropriate size for specified data resolution. To proceed, specify 'ignore_warnings=T' as function parameter")
      }else{warning("Smoothing window (h) may not be of appropriate size for specified data resolution")}
    }
    if(kb %in% 10:15 & (h >= 40 | h <= 15)){
      if(ignore_warnings==F){
        stop("Smoothing window (h) may not be of appropriate size for specified data resolution. To proceed, specify 'ignore_warnings=T' as function parameter")
      }else{warning("Smoothing window (h) may not be of appropriate size for specified data resolution")}
    }
    if(kb %in% 5:9 & (h >= 50 | h <= 20)){
      if(ignore_warnings==F){
        stop("Smoothing window (h) may not be of appropriate size for specified data resolution. To proceed, specify 'ignore_warnings=T' as function parameter")
      }else{warning("Smoothing window (h) may not be of appropriate size for specified data resolution")}
    }
    if(kb %in% 1:4 & (h >= 60 | h <= 25)){
      if(ignore_warnings==F){
        stop("Smoothing window (h) may not be of appropriate size for specified data resolution. To proceed, specify 'ignore_warnings=T' as function parameter")
      }else{warning("Smoothing window (h) may not be of appropriate size for specified data resolution")}
    }
  }

}
