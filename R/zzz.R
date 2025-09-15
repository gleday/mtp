#----------------------#
# Internal R functions
#----------------------#

.check_p_values <- function() {
  p <- get("p_values", envir = parent.frame())
  assert_that(
    all(map_lgl(p, is.number)),
    all(is.finite(p)),
    all(p >= 0),
    all(p <= 1)
  )
}

.check_k <- function() {
  k <- get("k", envir = parent.frame())
  m <- get("m", envir = parent.frame())
  assert_that(
    is.count(k + 1),
    is.finite(k),
    k >= 0,
    k < m
  )
}

.check_d <- function() {
  d <- get("d", envir = parent.frame())
  assert_that(
    is.number(d),
    is.finite(d),
    d >= 0,
    d < 1
  )
}

.check_lambda <- function() {
  lambda <- get("lambda", envir = parent.frame())
  assert_that(
    all(map_lgl(lambda, is.number)),
    all(is.finite(lambda)),
    all(lambda >= 0),
    all(lambda <= 1)
  )
}

.check_output <- function() {
  values <- c("p_values", "c_values", "decisions", "factors")
  output <- get("output", envir = parent.frame())
  output <- ifelse(
    are_equal(output, values),
    output[1L],
    output
  )
  assert_that(
    is.string(output),
    is.element(output, set = values)
  )
}

.check_c_value <- function() {
  cv <- get("c_value", envir = parent.frame())
  output <- get("output", envir = parent.frame())
  ifelse(
    is.null(cv),
    assert_that(
      is.element(output, set = c("p_values", "factors")),
      msg = "A value for `c_value` is required."
    ),
    assert_that(
      is.number(cv),
      is.finite(cv)
    )
  )
}

.check_weights <- function() {
  weights <- get("weights", envir = parent.frame())
  assert_that(
    all(map_lgl(weights, is.number)),
    all(is.finite(weights)),
    all(weights >= 0),
    all(weights <= 1),
    sum(weights) == 1
  )
}
