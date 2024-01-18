# Online Supplementary Materials for PhD Thesis Related to SIEVE
### R Codes and Related Supplementary Materials

Hongxiang Li

SIEVE: A statistical method for building a unified framework for the simultaneous testing of differential expression, variability and skewness of genes using RNA-Seq data. This framework adopts a compositional data analysis approach to modelling RNA-Seq count data, applies the centered log-ratio (CLR) transformation to convert them into continuous variables, and uses a skew-normal distribution to model them. 

Comparison of SIEVE's performance on differential expression (DE) and differential variability (DV) with existing methods involves the use of both simulated and real RNA-Seq data. Additionally, to demonstrate the practical usefulness of the newly developed DE, DV, and DS (differential skewness) methods by applying them to a real RNA-Seq data.


### Installation:
 `devtools::install_github("Divo-Lee/SIEVE")`
 
 
### Dependencies:
 `SIEVE` `R` package depends on the following packages: `sn`, `vioplot`, `stats`, `utils`, `grDevices`, `graphics`.
