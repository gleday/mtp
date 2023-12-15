# Overview
The mtp package provides functions to apply
various forms of statistical control on the
occurrence of Type I errors (incorrect rejections)
when conducting multiple statistical hypothesis tests.

mtp offers to use different error criteria
(familywise error rate, false rejection exceedance, ...)
to assess the likelihood of committing Type I errors,
as well as different methods (step-down, step-up)
to adjust the significance level of each test.

## Installation

To install **mtp** from R:

```R
# Install/load R package devtools
install.packages("devtools")
library(devtools)

# Install/load R package beam from github
install_github("gleday/mtp")
library(mtp)
```

## Why using multiple testing procedures?
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
level of each test (possibly using different methods).

The benefits of multiple testing procedures include:
- improved balance between incorrectly and
correctly rejected null hypotheses (i.e. between
Type I errors and statistical power).
- offer statistical guarantees regarding the extent of
incorrectly rejected null hypotheses
