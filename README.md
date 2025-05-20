
<!-- README.md is generated from README.Rmd. Please edit that file -->

# voomCLR

`voomCLR` allows effective differential cell composition analysis in
cell type count data. It leverages compositional transformations, and
adopts bias correction on the estimated effect sizes to correct for
compositional bias induced by such transformation. The uncertainty
involved in estimating the bias can be propagated in the statistical
inference via bootstrapping. Additionally, it accommodates proper
modeling of the mean-variance structure of the counts.

`voomCLR` relies on the `limma` package, and in fact re-uses code chunks
from the `limma` R package, which is available on Bioconductor at
<https://bioconductor.org/packages/release/bioc/html/limma.html>.

## Installation

Install the development version from GitHub using

``` r
devtools::install_github("koenvandenberge/voomCLR")
```

# Getting started

The vignette of voomCLR walks you through the basics of using `voomCLR`.
