# Overview
The mtp package provides functions to apply
various forms of statistical control on the
occurrence of Type I errors (incorrect rejections)
when conducting multiple statistical hypothesis tests.

Different error criteria (FWER, FRX) to assess the likelihood
of committing Type I errors, as well as different
methods (step-down, step-up) to adjust the
significance level of each test.

# Installation

To install **mtp** from R:

```R
# install/load R package devtools
install.packages("devtools")
library(devtools)

# install/load R package mtp from github
install_github("gleday/mtp")
library(mtp)
```

# Why using multiple testing procedures?
When performing a single statistical test,
the probability of committing an error when rejecting
the null hypothesis, is bounded by $\alpha$,
which is conventionally set at 5\%.

When performing multiple statistical tests,
the risk of committing errors increases with
the number of tests and is practically guaranteed
when the number of tests is large.

Multiple testing procedures allow to exert
statistical control on the occurrence of
Type I errors by capping the likelihood of
committing errors (in the sense of a criterion).
This is achieved by adjusting the significance
level of each test (different methods are available).

Overall, multiple testing procedures offer:
- to improve the balance between incorrectly and
correctly rejected null hypotheses (i.e. between
Type I errors and statistical power).
- statistical guarantees regarding the extent of
incorrectly rejected null hypotheses

# Other R packages

If the present package does not meet your needs,
perhaps one of the following packages will:

* [cherry](https://cran.r-project.org/package=cherry):
Multiple Testing Methods for Exploratory Research
* [DiscreteFDR](https://cran.r-project.org/package=DiscreteFDR):
FDR Based Multiple Testing Procedures with Adaptation for Discrete Tests
* [DiscreteFWER](https://cran.r-project.org/package=DiscreteFWER):
FWER-Based Multiple Testing Procedures with Adaptation for Discrete Tests
* [FAMT](https://cran.r-project.org/package=FAMT):
Factor Analysis for Multiple Testing
* [fdrci](https://cran.r-project.org/package=fdrci):
Permutation-Based FDR Point and Confidence Interval Estimation]
* [fdrtool](https://cran.r-project.org/package=fdrtool):
Estimation of (Local) False Discovery Rates and Higher Criticism
* [flip](https://cran.r-project.org/package=flip):
Multivariate Permutation Tests
* [FixSeqMTP](https://cran.r-project.org/package=FixSeqMTP):
Fixed Sequence Multiple Testing Procedures
* [pwrFDR](https://cran.r-project.org/package=pwrFDR):
FDR Power
* [FDX](https://cran.r-project.org/package=FDX):
False Discovery Exceedance Controlling Multiple Testing Procedures
* [gMCP](https://cran.r-project.org/package=gMCP):
Graph Based Multiple Comparison Procedures
* [graphicalMCP](https://cran.r-project.org/package=graphicalMCP):
Graphical Multiple Comparison Procedures
* [hommel](https://cran.r-project.org/package=hommel):
Methods for Closed Testing with Simes Inequality, in Particular Hommel's Method
* [MHTdiscrete](https://cran.r-project.org/package=MHTdiscrete):
Multiple Hypotheses Testing for Discrete Data
* [MHTmult](https://cran.r-project.org/package=MHTmult):
Multiple Hypotheses Testing for Multiple Families/Groups Structure
* [multcomp](https://cran.r-project.org/package=multcomp):
Simultaneous Inference in General Parametric Models
* [multtest](https://bioconductor.org/packages/multtest/):
Resampling-based multiple hypothesis testing
* [multxpert](https://cran.r-project.org/package=multxpert):
Common Multiple Testing Procedures and Gatekeeping Procedures
* [mutoss](https://cran.r-project.org/package=mutoss):
Unified Multiple Testing Procedures
* [nparcomp](https://cran.r-project.org/package=nparcomp):
Multiple Comparisons and Simultaneous Confidence Intervals
* [qvalue](https://bioconductor.org/packages/qvalue/):
Q-value estimation for false discovery rate control
* [sgof](https://cran.r-project.org/package=sgof):
Multiple Hypothesis Testing
* [someMTP](https://cran.r-project.org/package=someMTP):
Some Multiple Testing Procedures
* [StepwiseTest](https://cran.r-project.org/package=StepwiseTest):
Multiple Testing Method to Control Generalized Family-Wise Error Rate and False Discovery Proportion


