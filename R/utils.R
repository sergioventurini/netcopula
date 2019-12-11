#' Auxiliary function to recursively check NAs in a list.
#'
#' \code{check_list_na()} compares two lists and fills in the missing
#'   elements in the first with those included in the second. The
#'   comparison is recursive in the sense that the process is repeated for
#'   all lists included in those given.
#'
#' @param orig A list whose content must be checked.
#' @param des A list to use as a reference with which compare the first one.
#'
#' @return A list with all elements added.
#'
#' @author Sergio Venturini \email{sergio.venturini@unito.it}
#'
#' @export
check_list_na <- function(orig, des) {
  check_it <- function(o, d) {
    d.nm <- names(d)
    d.na <- is.na(match(d.nm, names(o)))
    if (any(d.na))
      o <- c(o, d[d.nm[which(d.na)]])

    return(o)
  }

  if (!is.list(orig))
    stop("the 'orig' argument must be a list")
  if (!is.list(des))
    stop("the 'des' argument must be a list")

  orig_new <- check_it(orig_new <- orig, des)

  for (el in 1:length(orig_new)) {
    if (is.list(orig_new[[el]]))
      orig_new[[el]] <- check_list_na(orig_new[[el]], des[[el]])
  }

  return(orig_new)
}

#' Check for suggested package (requireNamespace) and throw error if necessary
#'
#' @noRd
#' @param pkg Package name as a string.
#' @param min_version Optionally, a minimum version number as a string.
#' @return TRUE, invisibly, if no error is thrown.
#'
suggested_package <- function(pkg, min_version = NULL) {
  stopifnot(length(pkg) == 1, is.character(pkg))
  if (!requireNamespace(pkg, quietly = TRUE)) {
    stop(
      "Please install the ",
      pkg, " package to use this function.",
      call. = FALSE
    )
  }

  if (!is.null(min_version)) {
    stopifnot(is.character(min_version))
    if (utils::packageVersion(pkg) < package_version(min_version)) {
      stop(
        "Version >=", min_version, " of the ",
        pkg, " package is required to use this function.",
        call. = FALSE
      )
    }
  }

  invisible(TRUE)
}

#' Explicit and/or regex parameter selection
#'
#' @noRd
#' @param explicit Character vector of selected parameter names.
#' @param patterns Character vector of regular expressions.
#' @param complete Character vector of all possible parameter names.
#' @return Characeter vector of combined explicit and matched (via regex)
#'   parameter names, unless an error is thrown.
#'
select_pars <- function(explicit = character(), patterns = character(), complete = character()) {
  stopifnot(is.character(explicit),
            is.character(patterns),
            is.character(complete))

  if (!length(explicit) && !length(patterns))
    return(complete)

  if (length(explicit)) {
    if (!all(explicit %in% complete)) {
      not_found <- which(!explicit %in% complete)
      stop(
        "Some 'pars' don't match parameter names: ",
        paste(explicit[not_found], collapse = ", ")
      )
    }
  }

  if (!length(patterns)) {
    return(unique(explicit))
  } else {
    regex_pars <-
      unlist(lapply(seq_along(patterns), function(j) {
        grep(patterns[j], complete, value = TRUE)
      }))
    if (!length(regex_pars))
      stop("no matches for 'regex_pars'.", call. = FALSE)
  }

  unique(c(explicit, regex_pars))
}

choose_colors <- function(n) {
  all_clrs <- unlist(bayesplot::color_scheme_get())
  clrs <- switch(
    as.character(n),
    "1" = get_color("m"),
    "2" = get_color(c("l", "d")),
    "3" = get_color(c("l", "m", "d")),
    "4" = all_clrs[-c(2, 4)],
    "5" = all_clrs[-3],
    "6" = all_clrs,
    rep_len(all_clrs, n)
  )
  unname(rev(clrs))
}

# Access a subset of the scheme colors
#
# @param level A character vector of level names (see scheme_level_names()). The
#   abbreviations "l", "lh", "m", "mh", "d", and "dh" can also be used instead
#   of the full names.
# @return A character vector of color values.
#
# [Source: bayesplot]
get_color <- function(levels) {
  sel <- which(!levels %in% scheme_level_names())
  if (length(sel)) {
    levels[sel] <- sapply(levels[sel], full_level_name)
  }
  stopifnot(all(levels %in% scheme_level_names()))
  color_vals <- bayesplot::color_scheme_get()[levels]
  unlist(color_vals, use.names = FALSE)
}

full_level_name <- function(x) {
  switch(x,
         l = "light", lh = "light_highlight",
         m = "mid", mh = "mid_highlight",
         d = "dark", dh = "dark_highlight")
}

# Color scheme level names
#
# [Source: bayesplot]
scheme_level_names <- function() {
  c("light",
    "light_highlight",
    "mid",
    "mid_highlight",
    "dark",
    "dark_highlight")
}
