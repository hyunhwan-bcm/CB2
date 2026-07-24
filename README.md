[![](https://cranlogs.r-pkg.org/badges/CB2)](https://cran.r-project.org/package=CB2)
[![](https://img.shields.io/cran/l/CB2.svg)](https://cran.r-project.org/package=CB2/LICENSE)


## CB<sup>2</sup> <img src="man/figures/hexsticker.png" align="right" height="140"/>

CB<sup>2</sup>(CRISPRBetaBinomial) is a new algorithm for analyzing CRISPR data based on beta-binomial distribution. 
We provide CB<sup>2</sup> as a R package, and the interal algorithms of CB<sup>2</sup> are also implemented in [CRISPRCloud](https://crispr.nrihub.org/).

## Update

### Continuous and multivariable phenotypes

Original CB<sup>2</sup> asks whether guide abundance differs between two
groups. The additive CB<sup>2</sup>-Reg API asks whether guide abundance follows
a quantitative or adjusted sample design. It preserves the original
beta-binomial treatment of sequencing and between-library variation while
replacing the group comparison with an R model matrix. Existing two-group
functions are unchanged.

The relationship is methodological, not a claim that the two functions are
numerically identical. A saturated two-cell generalized least-squares
calculation reproduces original CB<sup>2</sup> after its weighted group
proportions and variances have already been computed. This is a legacy
compatibility identity, not strict nesting of the default `bbreg()` estimator.
Under common dispersion, local alternatives, and increasing numbers of
independent libraries, the default logit-scale statistic is first-order
equivalent. Increasing sequencing depth at fixed biological replication is not
sufficient. Finite-sample results can differ because CB<sup>2</sup> estimates
group-specific weights and uses Welch–Satterthwaite degrees of freedom, while
`bbreg()` uses one guide-wise dispersion across the design and residual sample
degrees of freedom. Its model-based covariance treats that dispersion estimate
as a fixed plug-in value and does not formally propagate its uncertainty.
Control-tail calibration requires all controls to share one residual degree of
freedom, which is satisfied when every guide uses the same complete design.

Use CB<sup>2</sup>-Reg for dose, time, ordered phenotypes, batch or donor
adjustment, interactions, and named contrasts:

```r
fit <- bbreg(
  count = guide_count,
  total = library_size,
  formula = ~ scale(dose) + batch,
  data = design
)

screen <- bb_screen(
  counts = count_matrix,
  totals = library_size,
  data = design,
  formula = ~ scale(dose) + batch,
  term = "scale(dose)",
  gene = guide_annotation$gene,
  ncores = 4
)
```

`bbreg()` fits a beta-binomial logit mean model and tests coefficients with a
Student t reference based on sample-level residual degrees of freedom.
`bb_contrast()` tests named linear contrasts, and `bb_screen()` applies the
model guide by guide with Benjamini-Hochberg correction. Weighted IRLS
cross-products and solves use the package's RcppArmadillo layer. See the
`barcs-regression` vignette for a complete example.

When a screen includes prespecified negative-control guides,
`bb_calibrate_controls()` can estimate a conservative empirical-null scale
from their t-statistic tail. It preserves raw inferential columns and does not
change coefficient estimates.

MAGeCK is kept external: use the official `mageck mle` executable for a
negative-binomial sensitivity analysis with the same sample design.

In short: use original CB<sup>2</sup> when the estimand is a two-condition
difference; use CB<sup>2</sup>-Reg when the estimand is a coefficient or
contrast; use a specialist model when the libraries are correlated partitions
or repeated measurements that require an explicit joint likelihood.

### Oct 1, 2025

Update the C++ dependency 

### Jun 7, 2022

A bug fix regarding issue #14. Thanks @DaneseAnna for reporting the issue.

### Dec 4, 2020

If you are experiencing an issue during the installation, it would be possible due to `multtest` package hasn't been installed. If so, please use the following snippet to install the package. 

```r
install.package("BiocManager") # can be omitted if you have installed the package
install.packages("multtest")
```

### May 26, 2020

* Regarding issue #9, CB<sup>2</sup> now provides logFC of gene-level analysis with two different modes. The default option is the same as the previous version, and setting `logFC` parameter value of `measure_gene_stats` to `gene` will provide the `logFC` calculate by gene-level CPMs.

### April 14, 2020

* Regarding issue #6, now users can use `join_count_and_design` function.

### December 16, 2019

* Regarding issue #4, CB<sup>2</sup> now supports gzipped FASTQ file.
* Regarding issue #5, `calc_mappability()` provide `total_reads` and `mapped_reads` columns.

### July 2, 2019 

There are several updates.

* We have change the function name for the sgRNA-level test to `measure_sgrna_stats`. The original name `run_estimation` has been *deprecated*.
* CB<sup>2</sup> now supports a `data.frame` with character columns. In other words, you can use 

## How to install

Currently CB<sup>2</sup> is now on `CRAN`, and you can install it using `install.package` function.

```r
install.package("CB2")
```

Installation Github version of CB<sup>2</sup> can be done using the following lines of code in your R terminal.

```r
install.packages("devtools")
devtools::install_github("hyunhwan-jeong/CB2")
```

Alternatively, here is a one-liner command line for the installation.

```
Rscript -e "install.packages('devtools'); devtools::install_github('hyunhwan-jeong/CB2')"
```

## A simple example how to use CB<sup>2</sup> in R

```r
FASTA <- system.file("extdata", "toydata",
                     "small_sample.fasta",
                     package = "CB2")
df_design <- data.frame()
for(g in c("Low", "High", "Base")) {
  for(i in 1:2) {
    FASTQ <- system.file("extdata", "toydata",
                         sprintf("%s%d.fastq", g, i), 
                         package = "CB2")
    df_design <- rbind(df_design, 
      data.frame(
        group = g, 
        sample_name = sprintf("%s%d", g, i),
        fastq_path = FASTQ, 
        stringsAsFactors = F)
      )
  }
}

MAP_FILE <- system.file("extdata", "toydata", "sg2gene.csv", package="CB2")
sgrna_count <- run_sgrna_quant(FASTA, df_design, MAP_FILE)
  
sgrna_stat <- measure_sgrna_stats(sgrna_count$count, df_design, 
                                  "Base", "Low", 
                                  ge_id = "gene",
                                  sg_id = "id")
gene_stat <- measure_gene_stats(sgrna_stat)
```

Or you could run the example with the following commented code.

```r
sgrna_count <- run_sgrna_quant(FASTA, df_design)
sgrna_stat <- measure_sgrna_stats(sgrna_count$count, df_design, "Base", "Low")
gene_stat <- measure_gene_stats(sgrna_stat)
```

More detailed tutorial is available [here](https://CRAN.R-project.org/package=CB2/vignettes/cb2-tutorial.html)!
