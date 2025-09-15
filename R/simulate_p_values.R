#' Simulate p-values
#'
#' @param m \[ `integer(1)` \]\cr
#' Number of hypotheses.
#' @param m0 \[ `integer(1)` \]\cr
#' Number of null hypotheses.
#' @param s1 \[ `numeric(1)` \]\cr
#' A non-negative number specifying the first shape parameter of
#' the Beta distribution representing the alternative hypothesis.
#' @param s2 \[ `numeric(1)` \]\cr
#' A non-negative number specifying the second shape parameter of
#' the Beta distribution representing the alternative hypothesis.
#'
#' @description
#' Generate p-values from a two-component mixture model with
#' uniform and Beta distributions.
#'
#' @details
#' The generating model is:
#'
#' \deqn{\displaystyle{
#' \pi_0 Beta(1, 1) + (1 - \pi_0) Beta(s_1, s_2)
#' }}
#' where \eqn{w_0 = m_0/m}.
#'
#' @return
#' \[ `tibble` \] \cr \cr
#' A tibble/data.frame with components:
#' * `hypothesis`: \[`integer`\] type of hypothesis
#' (0 = null and 1 = alternative).
#' * `p_value`:  \[`numeric`\] generated p-value.
#'
#' @references
#' Allison, D. B., Gadbury, G. L., Heo, M.,
#' Fernández, J. R., Lee, C. K., Prolla, T. A., &
#' Weindruch, R. (2002). A mixture model approach for
#' the analysis of microarray gene expression data.
#' Computational Statistics & Data Analysis, 39(1), 1-20.
#'
#' @export
simulate_p_values <- function(m, m0, s1, s2) {
  # true null and alternative hypotheses
  h <- rep.int(0L:1L, times = c(m0, m - m0))

  # simulate p-values
  p <- numeric(m)
  p[h == 0] <- rbeta(m0, shape1 = 1, shape2 = 1)
  p[h == 1] <- rbeta(m - m0, shape1 = s1, shape2 = s2)

  tibble::tibble(hypothesis = h, p_value = p)
}
