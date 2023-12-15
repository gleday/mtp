#----------------------#
# Internal R functions
#----------------------#

.check_p_value <- function() {
  p_value <- get("p_value", envir = parent.frame())
  assert_that(all(map_lgl(p_value, is.number)))
  assert_that(all(is.finite(p_value)))
  assert_that(all(between(p_value, 0, 1)))
}

.check_k <- function() {
  k <- get("k", envir = parent.frame())
  assert_that(is.count(k))
  assert_that(is.finite(k))
}

.check_d <- function() {
  d <- get("d", envir = parent.frame())
  assert_that(is.number(d))
  assert_that(is.finite(d))
  assert_that(between(d, 0, 1))
}

.check_lambda <- function() {
  lambda <- get("lambda", envir = parent.frame())
  assert_that(all(map_lgl(lambda, is.number)))
  assert_that(all(is.finite(lambda)))
  assert_that(all(between(lambda, 0, 1)))
}

.check_return <- function() {
  .return <- get(".return", envir = parent.frame())
  assert_that(is.string(.return))
  assert_that(is.element(.return, set = c("p", "a")))
}
