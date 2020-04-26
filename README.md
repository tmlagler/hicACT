# hicACT

## Introduction
<code>hicACT</code> is an R package to implement HiC-ACT. HiC-ACT, an aggregated Cauchy test (ACT) based approach, improves the detection of significant chromatin interactions by post-processing calling results from methods relying on the assumption that all chromatin interactions are statistically independent. HiC-ACT uses a linear combination of transformed *p*-values with non-negative weights to compute a smoothed *p*-value.

## Installation
We recommend installing from Github for the latest version of the code:
```r
install.packages("devtools")
devtools::install_github("tmlagler/hicACT")
library(hicACT)
```
## Model
Let $p_{ij}$ represent the *p*-value for chromatin interaction between bin *i* and bin *j* from a specific Hi-C peak calling method. Consider the null hypothesis that bin pair (*i*, *j*) is compatible with random chromatin looping. Define the HiC-ACT test statistic $T_{ACT_{ij}}$ as

<div style="text-align: center"><a href="https://www.codecogs.com/pdf.latex?eqnedit.php?latex=T_{ACT_{ij}}=\sum_{0&space;\leq&space;|m-i|&plus;|n-j|&space;\leq&space;h}&space;w_{mn}&space;tan\big\{(0.5-p_{mn})\pi\big\}" target="_blank"><img src="https://latex.codecogs.com/png.latex?T_{ACT_{ij}}=\sum_{0&space;\leq&space;|m-i|&plus;|n-j|&space;\leq&space;h}&space;w_{mn}&space;tan\big\{(0.5-p_{mn})\pi\big\}" title="T_{ACT_{ij}}=\sum_{0 \leq |m-i|+|n-j| \leq h} w_{mn} tan\big\{(0.5-p_{mn})\pi\big\}" /></a></div>

Here, $h$ is the local smoothing bandwidth. We take $w_{mn}$ to be the Gaussian kernel weight function, defined as:

<div style="text-align: center">$w_{mn}=\frac{\exp\big\{{\frac{-d^2_{mn}}{2}}\big\}}{\sum_{0 \leq |m'-i|+|n'-j| \leq h}\exp\big\{{\frac{d^2_{m'n'}}{2}}\big\}}$
$\text{ st } \sum w_{mn}=1,$  $\text{ and where } d^2_{mn} = (m-i)^2+(n-j)^2$</div>

Note that the p-value for the pair itself will also contribute to the statistic and thus the smoothed p-value.

$T_{ACT_{ij}}$ approximately follows a standard Cauchy distribution (*see HiC-ACT paper for details*). Therefore, the p-value for $T_{ACT_{ij}}$ can be approximated by:

<div style="text-align: center">$p^*_{ij} \approx 0.5-\big(tan^{-1}\{T_{ACT_{ij}}\}\big) \pi^{-1}$</div>

We can interpret $p_{ij}^*$ as the local neighborhood smoothed p-value. Intuitively, for a biologically meaningful chromatin interaction, all bin pairs in its neighborhood are more likely to have significant p-values. Thus, the combined p-value $p_{ij}^\*$ tends to be more significant and is driven by small p-values in its neighborhood. 

## Usage
```
hicACT(infile, kb, h, pthres,
    outdir, outname="ACT_adjusted",
    bin1col=2, bin2col=4, ccountcol=5, pcol=6,
    ignore_warnings=F)
```
Required paramters:

- **infile**: output file from a Hi-C chromatin interaction calling method, such as Fit-Hi-C (column names required)
- **kb**: data resolution in Kb, e.g. 10
- **h**: smoothing parameter, see 'Details' for suggestions
- **pthres**: p-value threshold for selecting which p-values to adjust using HiC-ACT

Optional parameters:

- **outdir**: output directory path (defaults to working directory)
- **outname**: desired output file name (defaults to "ACT_adjusted")
- **bin1col**: column number containing the first bin IDs (defaults to 2 for Fit-H-C)
- **bin2col**: column number containing the second bin IDs (defaults to 4 for Fit-Hi-C)
- **ccountcol**: column number containing the observed contact counts (defaults to 5 for Fit-Hi-C)
- **pcol**: column number containing the p-values (defaults to 6 for Fit-Hi-C)
- **ignore_warnings**: default is FALSE

## Suggestions
We suggest setting the smoothing parameter ($h$) based on the resolution of the data.

| Data Resolution (Kb) | Smoothing Parameter (h) |
|:--------------------:|:-----------------------:|
| 5 | 40 |
| 10 | 20 |
| 20 | 12 |
| 25 | 11 |
| 40 | 5 |

We also suggest setting the *p*-value threshold parameter based on the size of the data. For instance, in data with 0.5 billion raw reads <code>pthres</code> could be set to 0.1, but in larger data sizes (e.g. 2 billion raw reads), choosing <code>pthres</code>$\approx 1.0e-5$ is more reasonable and helps keep the computation time low.

## Example

```
hicACT(infile="fithic_chr22.ds0.1.spline_pass2.res10000.significances.txt.gz",
       kb=10, h=20, pthres=1.0e-16)
```



