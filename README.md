# distQTL
Distribution-based eQTL Mapping by Single-Cell RNA Sequencing


## Author
Alexander Coulter


## Maintainer
Alexander Coulter <coultera@tamu.edu>


## Description
distQTL is an `R` package for applying **dist**ribution-based **Q**uantitative
**T**rait **L**ocus mapping using single-cell RNA sequence data.  distQTL takes
cell-level expression data and SNP genotype data, and applies Fréchet regression
for distribution responses with corresponding partial F-test inference methods. 
This allows QTL mapping and characterization of distribution-wide effects.
distQTL sidesteps the donor level cell-cell correlation issue in single-cell RNA
sequence data, and obtains fast and numerically stable inference results. 
Unlike mixed effects approaches which fit a random intercept per donor, distQTL
computation burden scales linearly with number of donors, and does not suffer
convergence issues.


## Installation

As of April 1, 2025, distQTL passes R CMD checks on MacOS.  We will shortly
expand functionality to Windows OS.

To use `distQTL`, you need to install [`R`](https://cran.r-project.org/). To
enhance your user experience, you may use some IDE for it
(e.g. [`RStudio`](https://www.rstudio.com/)).

The development version of
[`distQTL`](https://github.com/alexandercoulter/distQTL) is available on GitHub.
You can download it with the help of the `devtools` package in `R` as follow:

```r
install.packages("devtools")
install.packages("data.table")
devtools::install_github("alexandercoulter/distQTL")
```

To build vignettes upon install, use
```r
devtools::install_github("alexandercoulter/distQTL", build_vignettes = TRUE)
```


## Demo Code & Vignettes
* [intro-distQTL vignette](https://alexandercoulter.github.io/distQTL/intro-distQTL.html) (loads .html file), with explanation and OneK1K subset example.

## Community Guidelines, Questions, and Issues

1.  If you have any questions with the package, feel free to email the maintainer at coultera@tamu.edu. We will try our best to reply as soon as we can.
2.  Report issues or problems with the software using GitHub’s [issue tracker](https://github.com/alexandercoulter/distQTL/issues).


## Citation


## Common questions
