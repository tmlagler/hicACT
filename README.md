# hicACT

## Introduction
<code>hicACT</code> is an R package to implement HiC-ACT. HiC-ACT, an aggregated Cauchy test (ACT) based approach, improves the detection of significant chromatin interactions by post-processing calling results from methods relying on the assumption that all chromatin interactions are statistically independent. HiC-ACT uses a linear combination of transformed *p*-values with non-negative weights to compute a smoothed *p*-value.

hicACT is maintained by Taylor Lagler [tmlagler@live.unc.edu] and Yun Li [yunli@med.unc.edu].

## Installation
We recommend installing from Github for the latest version of the code:
```r
install.packages("devtools")
devtools::install_github("tmlagler/hicACT")
library(hicACT)
```
## Model
Let <a href="https://www.codecogs.com/eqnedit.php?latex=\inline&space;\fn_phv&space;p_{ij}" target="_blank"><img src="https://latex.codecogs.com/svg.latex?\inline&space;\fn_phv&space;p_{ij}" title="p_{ij}" /></a> represent the *p*-value for chromatin interaction between bin *i* and bin *j* from a specific Hi-C peak calling method. Consider the null hypothesis that bin pair (*i*, *j*) is compatible with random chromatin looping. Define the HiC-ACT test statistic <a href="https://www.codecogs.com/eqnedit.php?latex=\inline&space;\fn_phv&space;T_{ACT_{ij}}" target="_blank"><img src="https://latex.codecogs.com/svg.latex?\inline&space;\fn_phv&space;T_{ACT_{ij}}" title="T_{ACT_{ij}}" /></a> as

<a href="https://www.codecogs.com/eqnedit.php?latex=\inline&space;\dpi{100}&space;\fn_phv&space;\large&space;T_{ACT_{ij}}=\sum_{0&space;\leq&space;|m-i|&plus;|n-j|&space;\leq&space;h}&space;w_{mn}&space;tan\big\{(0.5-p_{mn})\pi\big\}" target="_blank"><img src="https://latex.codecogs.com/svg.latex?\inline&space;\dpi{100}&space;\fn_phv&space;\large&space;T_{ACT_{ij}}=\sum_{0&space;\leq&space;|m-i|&plus;|n-j|&space;\leq&space;h}&space;w_{mn}&space;tan\big\{(0.5-p_{mn})\pi\big\}" title="\large T_{ACT_{ij}}=\sum_{0 \leq |m-i|+|n-j| \leq h} w_{mn} tan\big\{(0.5-p_{mn})\pi\big\}" /></a>

Here, *h* is the local smoothing bandwidth. We take <a href="https://www.codecogs.com/eqnedit.php?latex=\inline&space;\fn_phv&space;w_{mn}" target="_blank"><img src="https://latex.codecogs.com/svg.latex?\inline&space;\fn_phv&space;w_{mn}" title="w_{mn}" /></a> to be the Gaussian kernel weight function, defined as:

<a href="https://www.codecogs.com/eqnedit.php?latex=\inline&space;\fn_phv&space;\large&space;w_{mn}=\frac{\exp\big\{{\frac{-d^2_{mn}}{2}}\big\}}{\sum_{0&space;\leq&space;|m'-i|&plus;|n'-j|&space;\leq&space;h}\exp\big\{{\frac{d^2_{m'n'}}{2}}\big\}}" target="_blank"><img src="https://latex.codecogs.com/svg.latex?\inline&space;\fn_phv&space;\large&space;w_{mn}=\frac{\exp\big\{{\frac{-d^2_{mn}}{2}}\big\}}{\sum_{0&space;\leq&space;|m'-i|&plus;|n'-j|&space;\leq&space;h}\exp\big\{{\frac{d^2_{m'n'}}{2}}\big\}}" title="\large w_{mn}=\frac{\exp\big\{{\frac{-d^2_{mn}}{2}}\big\}}{\sum_{0 \leq |m'-i|+|n'-j| \leq h}\exp\big\{{\frac{d^2_{m'n'}}{2}}\big\}}" /></a>
<a href="https://www.codecogs.com/eqnedit.php?latex=\inline&space;\fn_phv&space;\large&space;\text{&space;st&space;}&space;\sum&space;w_{mn}=1,&space;\text{&space;and&space;where&space;}&space;d^2_{mn}&space;=&space;(m-i)^2&plus;(n-j)^2" target="_blank"><img src="https://latex.codecogs.com/svg.latex?\inline&space;\fn_phv&space;\large&space;\text{&space;st&space;}&space;\sum&space;w_{mn}=1,&space;\text{&space;and&space;where&space;}&space;d^2_{mn}&space;=&space;(m-i)^2&plus;(n-j)^2" title="\large \text{ st } \sum w_{mn}=1, \text{ and where } d^2_{mn} = (m-i)^2+(n-j)^2" /></a>

Note that the p-value for the bin pair itself will also contribute to the statistic and thus the smoothed p-value.

<a href="https://www.codecogs.com/eqnedit.php?latex=\inline&space;\fn_phv&space;T_{ACT_{ij}}" target="_blank"><img src="https://latex.codecogs.com/svg.latex?\inline&space;\fn_phv&space;T_{ACT_{ij}}" title="T_{ACT_{ij}}" /></a> approximately follows a standard Cauchy distribution (*see HiC-ACT paper for details*). Therefore, the p-value for <a href="https://www.codecogs.com/eqnedit.php?latex=\inline&space;\fn_phv&space;T_{ACT_{ij}}" target="_blank"><img src="https://latex.codecogs.com/svg.latex?\inline&space;\fn_phv&space;T_{ACT_{ij}}" title="T_{ACT_{ij}}" /></a> can be approximated by:

<a href="https://www.codecogs.com/eqnedit.php?latex=\inline&space;\fn_phv&space;\large&space;p^*_{ij}&space;\approx&space;0.5-\big(tan^{-1}\{T_{ACT_{ij}}\}\big)&space;\pi^{-1}" target="_blank"><img src="https://latex.codecogs.com/svg.latex?\inline&space;\fn_phv&space;\large&space;p^*_{ij}&space;\approx&space;0.5-\big(tan^{-1}\{T_{ACT_{ij}}\}\big)&space;\pi^{-1}" title="\large p^*_{ij} \approx 0.5-\big(tan^{-1}\{T_{ACT_{ij}}\}\big) \pi^{-1}" /></a>

We can interpret <a href="https://www.codecogs.com/eqnedit.php?latex=\inline&space;\fn_phv&space;p_{ij}^*" target="_blank"><img src="https://latex.codecogs.com/svg.latex?\inline&space;\fn_phv&space;p_{ij}^*" title="p_{ij}^*" /></a> as the local neighborhood smoothed p-value. Intuitively, for a biologically meaningful chromatin interaction, all bin pairs in its neighborhood are more likely to have significant p-values. Thus, the combined p-value <a href="https://www.codecogs.com/eqnedit.php?latex=\inline&space;\fn_phv&space;p_{ij}^*" target="_blank"><img src="https://latex.codecogs.com/svg.latex?\inline&space;\fn_phv&space;p_{ij}^*" title="p_{ij}^*" /></a> tends to be more significant and is driven by small p-values in its neighborhood. Through the properties of the aggregated Cauchy combination test, HiC-ACT specifically accounts for the inherent correlation between the contact frequency of neighboring pairs and maintains the benefit of not needing information about the correlation structure.

Further details on the Cauchy combination test underlying the HiC-ACT method can be found in the papers by Liu et. al ([10.1080/01621459.2018.1554485](https://doi.org/10.1080/01621459.2018.1554485), [10.1016/j.ajhg.2019.01.002](https://doi.org/10.1016/j.ajhg.2019.01.002)). Details on the HiC-ACT method can be found in our AJHG paper [10.1016/j.ajhg.2021.01.009](https://doi.org/10.1016/j.ajhg.2021.01.009).

## Usage
```
hicACT(infile, kb, h, pthres,
    outdir, outname="ACT_adjusted",
    bin1col=2, bin2col=4, ccountcol=5, pcol=6,
    ignore_warnings=F)
```
Required paramters:

- **infile**: output file from a Hi-C chromatin interaction calling method, such as Fit-Hi-C/FitHiC2 (column names required)
- **kb**: data resolution in Kb, e.g. 10
- **h**: smoothing parameter, see 'Details' for suggestions
- **pthres**: p-value threshold for selecting which p-values to adjust using HiC-ACT

Optional parameters:

- **outdir**: output directory path (defaults to working directory)
- **outname**: desired output file name (defaults to "ACT_adjusted")
- **bin1col**: column number containing the first bin IDs (defaults to 2 for Fit-Hi-C/FitHiC2)
- **bin2col**: column number containing the second bin IDs (defaults to 4 for Fit-Hi-C/FitHiC2)
- **ccountcol**: column number containing the observed contact counts (defaults to 5 for Fit-Hi-C/FitHiC2)
- **pcol**: column number containing the p-values (defaults to 6 for Fit-Hi-C/FitHiC2)
- **ignore_warnings**: default is FALSE

Note that hicACT() is intended to be run separately for each chromosome.

## Suggestions
We suggest setting the smoothing parameter (*h*) based on the resolution of the data.

| Data Resolution (Kb) | Smoothing Parameter (*h*) |
|:--------------------:|:-----------------------:|
| 5 | 40 |
| 10 | 20 |
| 20 | 12 |
| 25 | 11 |
| 40 | 5 |

We also suggest setting the *p*-value threshold parameter based on the size of the data. For instance, we suggest choosing<code>pthres</code><a href="https://www.codecogs.com/eqnedit.php?latex=\inline&space;\fn_phv&space;\small&space;\approx&space;1.0e-3" target="_blank"><img src="https://latex.codecogs.com/svg.latex?\inline&space;\fn_phv&space;\small&space;\approx&space;1.0e-3" title="\small \approx 1.0e-3" /></a> for
data with fewer than 1 billion raw reads, and choosing <code>pthres</code><a href="https://www.codecogs.com/eqnedit.php?latex=\inline&space;\fn_phv&space;\small&space;\approx&space;1.0e-6" target="_blank"><img src="https://latex.codecogs.com/svg.latex?\inline&space;\fn_phv&space;\small&space;\approx&space;1.0e-6" title="\small \approx 1.0e-6" /></a> for data with greater than 1 billion raw reads. This pre-filtering helps keep the computation time low. 

## Example
The test data used in the following example is the output of FitHiC2 applied to the GM12878 10Kb Hi-C data (Rao, et al. 2014) chromosome 22 down-sampled to ~0.5 billion raw reads, then subsetted to only include 10,000 pairwise interactions. The purpose of this test data is to test that the package has been installed correctly. Local run time for the following example is ~10 seconds. <code>hicACT</code> should not be run locally on large data sets.
```r
library(hicACT)
test_file <- paste0(path.package("hicACT"), "/test_data.txt")
head(fread(test_file))
hicACT(infile=test_file, kb=10, h=20, thres=0.1)
```

## Citation
Lagler TM, Abnousi A, Hu M, Yang Y, Li Y. HiC-ACT: improved detection of chromatin interactions from Hi-C data via aggregated Cauchy test. Am J Hum Genet. 2021 Feb 4;108(2):257-268. doi: [10.1016/j.ajhg.2021.01.009](https://doi.org/10.1016/j.ajhg.2021.01.009). PMID: 33545029.


