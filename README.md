# distQTL
Distribution-based eQTL Mapping by Single-Cell RNA Sequencing


## Author
Alexander Coulter


## Maintainer
Alexander Coulter <coultera@tamu.edu>


## Description
distQTL is a package for applying **dist**ribution-based **Q**uantitative **T**rait **L**ocus mapping using single-cell RNA sequence data.  distQTL takes cell-level expression data and SNP genotype data, and applies Fréchet regression for distribution responses with corresponding partial F-test inference methods.  This allows QTL mapping and characterization of distribution-wide effects.  distQTL sidesteps the donor level cell-cell correlation issue in single-cell RNA sequence data, and obtains fast and numerically stable inference results.  Unlike mixed effects approaches which fit a random intercept per donor, distQTL computation burden scales linearly with number of donors, and does not suffer convergence issues.


## Installation
```r
install.packages("devtools")
install.packages("data.table")
devtools::install_github("alexandercoulter/distQTL")
```


## Demo code & Vignettes
* distQTL
  * [Vignette]


## Questions & Issues
If you have any questions with the package, feel free to email the maintainer at coultera@tamu.edu. We will try our best to reply as soon as we can.


## Citation


## Common questions
