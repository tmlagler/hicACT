##' HiC-ACT: Aggregated Cauchy Test approach for improved chromatin interaction calling from Hi-C data
##' @param infile output file from a Hi-C chromatin interaction calling method, such as Fit-Hi-C (column names required)
##' @param kb data resolution in Kb, e.g. 10
##' @param h smoothing parameter, see 'Details' for suggestions
##' @param pthres p-value threshold for selecting which p-values to adjust using HiC-ACT
##' @param outdir output directory path
##' @param outname output file name
##' @param bin1col column number containing the first bin IDs (defaults to 2 for Fit-Hi-C)
##' @param bin2col column number containing the second bin IDs (defaults to 4 for Fit-Hi-C)
##' @param ccountcol column number containing the observed contact counts (defaults to 5 for Fit-Hi-C)
##' @param pcol column number containing the p-values (defaults to 6 for Fit-Hi-C)
##' @param ignore_warnings default is FALSE
##' @section Details: smoothing parameter (h) is based on data resolution. Suggestions: (kb=10, h=20), (kb=20, h=12), (kb=25, h=10)
##' @return A zipped text file containing the original input file with an appendended column of HiC-ACT p-values for all p-values processed (those less than the specified threshold)
##' @import data.table
##' @import fastmatch
##' @importFrom R.utils gzip
##' @export hicACT

hicACT <- function(infile, kb, h, pthres, outdir, outname="ACT_adjusted",
                   bin1col=2, bin2col=4, ccountcol=5, pcol=6,
                   ignore_warnings=F){

  # read in Hi-C peak calling output (eg from Fit-Hi-C)
  indata <- data.table::fread(infile, header=T)

  # force numeric inputs to be numeric if read in as strings
  kb <- as.numeric(kb)
  h <- as.numeric(h)
  pthres <- as.numeric(pthres)
  bin1col <- as.numeric(bin1col)
  bin2col <- as.numeric(bin2col)
  ccountcol <- as.numeric(ccountcol)
  pcol <- as.numeric(pcol)

  # store necessary column names for data.table operations
  cnames <- names(indata)
  bin1 <- as.name(cnames[bin1col])
  bin2 <- as.name(cnames[bin2col])
  ccount <- as.name(cnames[ccountcol])
  pval <- as.name(cnames[pcol])

  # check data input and specified parameters
  #check_input(indata, kb, h, pthres, bin1, bin2, ccount, pval, ignore_warnings)

  # remove all zero valued contactCounts
  data <- indata[eval(ccount) > 0]

  # adds simple identifiers for each fragment
  Kb <- kb*1e3
  setDT(data)[, i := (eval(bin1)+Kb/2)/Kb][, j := (eval(bin2)+Kb/2)/Kb]

  # select bin pairs for which to compute HiC-ACT test statistic and p-value
  outdata <- data[eval(pval) < pthres]
  ij_set <- outdata[,c("i","j")]


  # compute test statistic Tact and return corresponding p-value
  Tact <- function(i,j){

    # reduce size of data for quicker computation
    reduce <- abs(data$i - i) <= h
    data_red <- data[reduce,]

    # identify possible (m,n) pairs
    m_pos <- seq(i-h,i+h)
    n_pos <- seq(j-h,j+h)
    pair_pos <- data.table(expand.grid(m_pos, n_pos))
    colnames(pair_pos) <- c("m", "n")
    inrange <- abs(i-pair_pos$m)+abs(j-pair_pos$n) <= h
    pairs <- pair_pos[inrange,]

    # only (m,n) pairs with pvalues in data
    reduce.dim <- (data_red$i %fin% pairs$m) & (data_red$j %fin% pairs$n)
    small_data <- data_red[reduce.dim,]
    is.included <- paste(small_data$i, small_data$j, sep=":") %fin% paste(pairs$m, pairs$n, sep=":")
    p_pairs <- small_data[is.included,]

    m <- p_pairs$i
    n <- p_pairs$j
    pvals <- p_pairs[,eval(pval)]

    maxp <- max(pvals[pvals!=1])
    minp <- min(pvals[pvals!=0])

    # Gaussian kernel weight function
    d2 <- (m-i)^2 + (n-j)^2
    weights <- exp(-d2/2)/sum(exp(-d2/2))

    # Compute T_ACT statistic
    # check if pvals exactly equal to 0 or to 1
    # estimate tan() for very small pvalues
    is.zero <- pvals == 0
    is.one <- pvals == 1

    pvals[is.zero] <- minp
    pvals[is.one] <- maxp

    is.small <- pvals < 1e-16
    Tact <- sum(weights[is.small]/pvals[is.small]/pi)
    Tact <- Tact + sum(weights[!is.small]*tan((0.5-pvals[!is.small])*pi))

    # Compute p-value for statistic
    # check if the test statistic is very large
    pvalue <- ifelse(Tact > 1e+15,
                     (1/Tact)/pi, # use approximation for large p
                     1-pcauchy(Tact)) # use Cauchy distribution otherwise
    return(pvalue)
  }

  # compute all desired Tact pvalues
  n_pairs <- dim(ij_set)[1]
  Tact_pvals <- array(NA)
  for(k in 1:n_pairs){
    Tact_pvals[k] <- Tact(as.numeric(ij_set$i[k]), as.numeric(ij_set$j[k]))
  }

  # add HiC-ACT p-values to input file
  # only includes bin pairs of interest (post filtering)
  outdata <- data[eval(pval) < pthres]
  setDT(outdata)[, ACT_pvalue := Tact_pvals]
  # drop (i, j) columns
  outdata[, c("i","j"):=NULL]

  # write results to zipped text file
  data.table::fwrite(outdata, paste0(outdir, outname, ".txt"), sep='\t', col.names = T, row.names = F)
  R.utils::gzip(paste0(outdir, outname, ".txt"),
       destname=paste0(outdir, outname, ".txt.gz"),
       overwrite=TRUE)

}

