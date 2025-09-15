#' Estimate the number or proportion of null hypotheses
#'
#' @param lambda \[ `numeric()` \]\cr
#' A numeric vector containing one or more values
#' between 0 and 1 for the tuning parameter.
#'
#' @inheritParams fwer_bon
#'
#' @description
#' The functions provides Storey's estimator of
#' the number and proportion of true null hypotheses.
#'
#' @details
#' Storey's estimators of the
#' number and proportion of true null
#' hypotheses are respectively:
#' \deqn{\displaystyle{
#' \hat{m}_0 = \min\left(
#'  \frac{\sum_{j=1}^{m}{1_{p_j > \lambda} + 1}}
#'  {(1-\lambda)}, m
#' \right)
#' }}
#'
#' and
#'
#' \deqn{\displaystyle{
#' \hat{\pi}_0 = \frac{\hat{m}_0}{m}
#' }}
#'
#' @return
#' \[ `numeric()` \] \cr \cr
#' A numeric vector the same length as `lambda`.
#'
#' @family adaptive
#'
#' @author Gwenaël G.R. Leday
#'
#' @references
#' Storey, J. D. (2002). A direct approach to false discovery rates.
#' Journal of the Royal Statistical Society:
#' Series B (Statistical Methodology), 64(3), 479-498.
#'
#' @export
m0 <- function(p_values, lambda = 0.5) {

  # check arguments
  .check_p_values()
  .check_lambda()

  # get Storey's estimator
  m <- length(p_values)
  m0_hat <- (sum(p_values > lambda) + 1) / (1 - lambda)
  m0_hat <- pmin(m0_hat, m)

  # output
  if (length(lambda) == 1) {
    return(m0_hat)
  }
  names(lambda) <- lambda
  map_dbl(lambda, m0, p_values = p_values)
}
