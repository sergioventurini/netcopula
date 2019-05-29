#' @export
build_output <- function(a) {
  nparams <- ifelse(length(dim(a)) < 3, 1, dim(a)[length(dim(a))])
  output <- matrix(NA, ncol = 8, nrow = nparams)
  if (length(dim(a)) == 2) {
    a <- array(a, c(dim(a), 1))
  }
  for (i in 1:nparams) {
    ai <- a[, , i, drop = FALSE]
    quantiles <- quantile(as.vector(ai), probs = c(.025, .25, .5, .75, .975))
    output[i, ] <- c(mean(ai), median(ai), sd(as.vector(ai)), quantiles)
  }
  dimnames(output) <- list(dimnames(a)[[3]], c("mean", "median", "sd", "2.5%", "25%", "50%", "75%", "97.5%"))
  return(output)
}

#' @export
sigma2corr <- function(sigma) {
  nn <- length(sigma)
  M <- (-1 + sqrt(1 + 8*nn))/2
  sigma_mat <- matrix(NA, nrow = M, ncol = M)
  sigma_mat[lower.tri(sigma_mat, diag = TRUE)] <- sigma
  sigma_mat[upper.tri(sigma_mat)] <- sigma_mat[lower.tri(sigma_mat)]
  corr <- cov2cor(sigma_mat)
  if (M > 1) {
    out <- corr[lower.tri(corr)]
  } else {
    out <- 1
  }

  return(out)
}

#' @export
add_corr <- function(res_mcmc) {
  Sigma_M <- res_mcmc[, grep("Sigma_M", dimnames(res_mcmc)[[2]]), drop = FALSE]
  Corr_M <- t(apply(Sigma_M, 1, sigma2corr))
  if (nrow(Corr_M) == 1) {
    Corr_M <- t(Corr_M)
  }
  nn <- dim(Sigma_M)[2]
  M <- (-1 + sqrt(1 + 8*nn))/2
  if (M > 1) {
    nm <- gsub("Sigma", "Corr", dimnames(Sigma_M)[[2]][-cumsum(1:M)])
  } else {
    nm <- gsub("Sigma", "Corr", dimnames(Sigma_M)[[2]])
  }
  dimnames(Corr_M) <- list(NULL, nm)
  new_res_mcmc <- as.mcmc(cbind(res_mcmc, Corr_M))

  return(new_res_mcmc)
}

#' @export
col.alpha <- function(acol, alpha = 0.2) {
    acol <- col2rgb(acol)
    acol <- rgb(acol[1]/255, acol[2]/255, acol[3]/255, alpha)
    acol
}

#' @export
plot_est <- function(res_sum, true_values, param, main = NULL) {
  idx <- grep(paste0(param, "["), dimnames(res_sum$quantiles)[[1]], fixed = TRUE)
  tv_tmp <- as.numeric(na.omit(as.numeric(true_values)))
  plot(res_sum$statistics[idx, 1], tv_tmp,
    ylim = c(min(as.numeric(res_sum$quantiles[idx, c(1, 5)]), true_values, na.rm = TRUE), max(as.numeric(res_sum$quantiles[idx, c(1, 5)]), true_values, na.rm = TRUE)),
    xlab = "posterior mean estimate",
    ylab = "true values",
    main = ifelse(is.null(main), param, main),
    type = "n")
  for (i in idx) {
    segments(res_sum$statistics[i, 1], res_sum$quantiles[i, 1],
      res_sum$statistics[i, 1], res_sum$quantiles[i, 5],
      col = "gray", lwd = 1)
  }
  abline(h = 0, lty = 2, lwd = .2, col = 1)
  abline(0, 1, lty = 2, lwd = .2, col = 1)
  points(res_sum$statistics[idx, 1], res_sum$statistics[idx, 1],
    pch = 20, col = "lightblue")
  points(res_sum$statistics[idx, 1], tv_tmp, pch = 20, col = "orange")
}

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
#' @author Sergio Venturini \email{sergio.venturini@sdabocconi.it}
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
