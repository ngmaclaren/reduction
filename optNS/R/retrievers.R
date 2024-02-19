## Convenience functions for getting information back

#' Retrieve information from lists of sentinel node sets
#'
#' Convenience functions to retrieve information from lists of sentinel node sets, typically the output of make_dataset().
#'
#' @param dl The list of node sets
#'
#' @return A matrix (vs, ks, ws if not NULL) or numeric vector (error).
#' @name convenience
NULL

#' @rdname convenience
#' @export
get_vs <- function(dl) sapply(dl, `[[`, "vs")

#' @rdname convenience
#' @export
get_error <- function(dl) sapply(dl, `[[`, "error")

#' @rdname convenience
#' @export
get_ks <- function(dl) sapply(dl, `[[`, "ks")

#' @rdname convenience
#' @export
get_ws <- function(dl) sapply(dl, `[[`, "ws")
