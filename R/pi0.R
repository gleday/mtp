#' @rdname m0
#'
#' @export
pi0 <- function(p, lambda = 0.5) {
  m0(p = p, lambda = lambda) / length(p)
}
