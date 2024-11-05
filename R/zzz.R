#----------------------#
# Internal R functions
#----------------------#

.check_p <- function() {
  p <- get("p", envir = parent.frame())
  assert_that(all(map_lgl(p, is.number)))
  assert_that(all(is.finite(p)))
  assert_that(all(between(p, 0, 1)))
}

.check_k <- function() {
  k <- get("k", envir = parent.frame())
  assert_that(is.count(k + 1))
  assert_that(is.finite(k))
  m <- get("m", envir = parent.frame())
  assert_that(k >= 0)
  assert_that(k < m)
}

.check_d <- function() {
  d <- get("d", envir = parent.frame())
  assert_that(is.number(d))
  assert_that(is.finite(d))
  assert_that(d >= 0)
  assert_that(d < 1)
}

.check_lambda <- function() {
  lambda <- get("lambda", envir = parent.frame())
  assert_that(all(map_lgl(lambda, is.number)))
  assert_that(all(is.finite(lambda)))
  assert_that(all(between(lambda, 0, 1)))
}

.check_output <- function() {
  output <- get("output", envir = parent.frame())
  assert_that(is.string(output))
  assert_that(is.element(output, set = c("p", "a")))
}

.check_alpha <- function() {
  alpha <- get("alpha", envir = parent.frame())
  if (!is.null(alpha)) {
    assert_that(is.number(alpha))
    assert_that(is.finite(alpha))
    assert_that(between(alpha, 0, 1))
  }
}

.check_gamma <- function() {
  gamma <- get("gamma", envir = parent.frame())
  if (!is.null(gamma)) {
    assert_that(is.number(gamma))
    assert_that(is.finite(gamma))
    m <- get("m", envir = parent.frame())
    assert_that(between(gamma, 1, m))
  }
}
