setClassUnion(name = "null_or_data_frame", members = c("data.frame", "NULL"))

#' An S4 class to represent a set of data for use in the \code{netcopula}
#' package.
#'
#' @slot type A length-one character vector representing the data type and must
#' be either \code{"binary"}, \code{"count"} or \code{"continuous"}. Currently
#' the package supports only the "binary" data type.
#' @slot n_study A length-one numeric vector representing the number of
#' studies.
#' @slot n_outcomes A length-one numeric vector representing the number of
#' studies.
#' @slot n_datapoints A length-one numeric vector representing the number of
#' data points in the sample (i.e., the overall number of arms for the trials
#' in the sample).
#' @slot n_treatments A length-one numeric vector representing the total
#' number of treatments involved in the analysis.
#' @slot study_id A four-columns data frame providing studies-specific data.
#' The four columns must be named \code{studyid}, \code{study}, \code{narms}
#' and \code{baseline}, with the second one that needs to be a progressive
#' number that match with those in the \code{study_data} slot. If NULL, a
#' default set of IDs will be created.
#' @slot study_data A data frame containing the outcome data for each study.
#' Its first five columns must be named \code{study}, \code{arm}, \code{trt},
#' \code{baseline} and \code{narms}. Values in the \code{study} column must
#' match with those in the same column of the \code{study_id} slot. Outcomes
#' must occupy the last \code{n_outcomes} columns. Finally, the columns right
#' before the outcomes must provide the sample sizes for each study.
#' @slot ref_trt A length-one integer vector providing the reference treatment
#' (which should not be confused with the baseline treatment for each study).
#'
#' @name nc_data-class
#' @rdname nc_data-class
#' @aliases nc_data
#'
#' @author Sergio Venturini \email{sergio.venturini@unito.it}
#'
#' @exportClass nc_data
#'
#' @examples
#' showClass("nc_data")
#' getSlots("nc_data")
setClass(Class = "nc_data",
  slots = list(
    type = "character",
    n_study = "numeric",
    n_outcomes = "numeric",
    n_datapoints = "numeric",
    n_treatments = "numeric",
    study_id = "null_or_data_frame",
    study_data = "data.frame",
    ref_trt = "numeric"
  ),
  validity = function(object) {
    if (object@type == "count" | object@type == "continuous") {
      stop("[nc_data validation]: Currently only the 'binary' data type is supported.")
    }
    if (!(object@type %in% c("binary", "count", "continuous"))) {
      stop("[nc_data validation]: Slot 'type' must be either 'binary', 'count' or 'continuous'.")
    }

    if (object@n_study < 2) {
      stop("[nc_data validation]: Slot 'n_study' must be larger than 1.")
    }
    if (object@n_study - as.integer(object@n_study) != 0) {
      stop("[nc_data validation]: Slot 'n_study' must be integer.")
    }

    if (object@n_outcomes < 1) {
      stop("[nc_data validation]: Slot 'n_outcomes' must be at least 1.")
    }
    if (object@n_outcomes - as.integer(object@n_outcomes) != 0) {
      stop("[nc_data validation]: Slot 'n_outcomes' must be integer.")
    }

    if (object@n_treatments < 2) {
      stop("[nc_data validation]: Slot 'n_treatments' must be at least 2.")
    }
    if (object@n_treatments - as.integer(object@n_treatments) != 0) {
      stop("[nc_data validation]: Slot 'n_treatments' must be integer.")
    }

    if (object@n_datapoints - as.integer(object@n_datapoints) != 0) {
      stop("[nc_data validation]: Slot 'n_datapoints' must be integer.")
    }

    if (nrow(object@study_id) != object@n_study) {
      stop(paste0("[nc_data validation]: Slot 'study_id' must have exactly n_study = ", object@n_study, " rows."))
    }
    if (ncol(object@study_id) != 4) {
      stop("[nc_data validation]: Slot 'study_id' must have exactly four columns reporting for each study the study ID and its progressive number, the number of arms in each study and the corresponding baseline treatment.")
    }
    if (any(colnames(object@study_id) != c("studyid", "study", "narms", "baseline"))) {
      stop("[nc_data validation]: Slot 'study_id' columns must be named 'studyid', 'study', 'narms' and 'baseline' respectively.")
    }
    if (any((sort(object@study_id[, "study"]) - (1:object@n_study)) != 0)) {
      stop(paste0("[nc_data validation]: Slot column 'study' in the 'study_id' data fame must be a progressive number from 1 to n_study = ", object@n_study, "."))
    }

    if (nrow(object@study_data) != object@n_datapoints) {
      stop(paste0("[nc_data validation]: Slot 'study_data' must have exactly n_datapoints = ", object@n_datapoints, " rows."))
    }
    if (ncol(object@study_data) != (6 + 2*object@n_outcomes)) {
      stop(paste0("[nc_data validation]: Slot 'study_data' must have exactly ", (6 + 2*object@n_outcomes), " columns (see ?nc_data)."))
    }
    if (any(colnames(object@study_data)[2:6] != c("study", "arm", "trt", "baseline", "narms"))) {
      stop("[nc_data validation]: Slot 'study_data' first five columns must be named 'study', 'arm', 'trt', 'baseline' and 'narms' respectively.")
    }
    if (any(colnames(object@study_data)[7:(6 + object@n_outcomes)] != paste0("n", 1:object@n_outcomes))) {
      stop("[nc_data validation]: Slot 'study_data' columns providing the study sample sizes must be progressively named 'n1', 'n2',...")
    }
    if (any(colnames(object@study_data)[(6 + object@n_outcomes + 1):(6 + 2*object@n_outcomes)] != paste0("y", 1:object@n_outcomes))) {
      stop("[nc_data validation]: Slot 'study_data' columns providing the outcome data must be progressively named 'y1', 'y2',...")
    }

    if (!is.integer(object@ref_trt)) {
      object@ref_trt <- as.integer(object@ref_trt)
    }
    if (object@ref_trt < 1) {
      stop("[nc_data validation]: Slot 'ref_trt' must be a positive integer.")
    }

    return(TRUE)
  }
)

#' @describeIn nc_data Create an instance of the \code{\link{nc_data}} class
#' using new/initialize.
#'
#' @param .Object Prototype object from the class \code{\link{nc_data}}.
#' @param type A length-one character vector representing the data type (see
#' the corresponding slot description for more details).
#' studies (see the corresponding slot description for more details).
#' @param n_study A length-one numeric vector representing the number of
#' studies (see the corresponding slot description for more details).
#' @param study_id A three-columns data frame providing the studies ID (see the
#' corresponding slot description for more details).
#' @param study_data A data frame containing the outcome data for each study
#' (see the corresponding slot description for more details).
#' @param ref_trt A length-one integer vector providing the reference treatment
#' (which should not be confused with the baseline treatment for each study).
#' @param \ldots Optional arguments to the function in its next call.
#'
#' @author Sergio Venturini \email{sergio.venturini@unito.it}
#'
#' @aliases initialize,nc_data-method
#' @aliases nc_data-initialize
#' 
#' @importFrom methods initialize
#' @exportMethod initialize
#' 
#' @export
setMethod("initialize", signature(.Object = "nc_data"),
  function(.Object, ..., type, n_study, study_id, study_data, ref_trt) {
      if (!(is.data.frame(study_id) | is.null(study_id))) {
        stop("[nc_data validation]: Slot 'study_id' must be either a data frame or NULL.") # not really a validation, rather an initialization check! ;)
      }

      if (is.null(study_id)) {
        id_arms <- match(unique(study_data[, "study"]), study_data[, "study"])
        study_id <- data.frame(studyid = 1:n_study, study = 1:n_study, narms = study_data[id_arms, "narms"], baseline = study_data[id_arms, "baseline"])
      }
      study_data <- cbind(studyid = study_id[match(study_data[, 1], study_id[, 2]), 1], study_data)

      if (type == "binary" | type == "count") {
        if (any(sapply(study_data[, -(1:6)], class) != "integer")) {
          study_data[, -(1:6)] <- lapply(study_data[, -(1:6)], function(x) as.integer(round(x, digits = 0)))
        }
      }

      if(is.null(ref_trt)) {
        ref_trt <- as.integer(1)
      }

      callNextMethod(.Object, ..., type = type, n_study = n_study, study_id = study_id, study_data = study_data, ref_trt = ref_trt) # see Gentleman (2009), p. 9
  }
)

#' @describeIn nc_data Show an instance of the \code{\link{nc_data}} class.
#'
#' @param object An object of class \code{\link{nc_data}}.
#'
#' @aliases show,nc_data-method
#' 
#' @author Sergio Venturini \email{sergio.venturini@unito.it}
#'
#' @aliases show,nc_data-method
#' @aliases nc_data-show
#' 
#' @importFrom methods show
#' @exportMethod show
setMethod("show", signature(object = "nc_data"),
  function(object) {
    display_limit <- 25
    cat("Object of class 'nc_data':\n")
    cat("Slots:\n")
    cat("@type:            ", object@type, "\n")
    cat("@n_study:         ", sprintf("%.0f", object@n_study), "\n")
    cat("@n_outcomes:      ", sprintf("%.0f", object@n_outcomes), "\n")
    cat("@n_datapoints:    ", sprintf("%.0f", object@n_datapoints), "\n")
    cat("@n_treatments:    ", sprintf("%.0f", object@n_treatments), "\n")
    cat("@ref_trt:         ", sprintf("%.0f", object@ref_trt), "\n")
    cat("@study_id:\n")
    if (nrow(object@study_id) > display_limit) {
        show(head(object@study_id, n = display_limit))
        cat("   (only first", display_limit, "rows of", nrow(object@study_id), "displayed)\n")
    } else {
        show(object@study_id)
    }
    cat("@study_data:\n")
    if (nrow(object@study_data) > display_limit) {
        show(head(object@study_data, n = display_limit))
        cat("   (only first", display_limit, "rows of", nrow(object@study_data), "displayed)\n")
    } else {
        show(object@study_data)
    }
  }
)

#' @describeIn nc_data Provide a summary of a \code{\link{nc_data}} class
#' instance.
#'
#' @author Sergio Venturini \email{sergio.venturini@unito.it}
#'
#' @aliases summary,nc_data-method
#' @aliases nc_data-summary
#' 
#' @exportMethod summary
setMethod("summary", signature(object = "nc_data"),
  function(object) {
    out_vec <- c("General Information:", "--------------------")
    out_char <- sprintf("%.0f", object@n_study)
    out_vec <- c(out_vec, paste("Number of studies:                       ", out_char))
    out_char <- sprintf("%.0f", object@n_treatments)
    out_vec <- c(out_vec, paste("Overall number of treatments:            ", out_char))
    out_char <- object@type
    out_vec <- c(out_vec, paste("Outcome type:                            ", out_char))
    out_char <- sprintf("%.0f", object@n_outcomes)
    out_vec <- c(out_vec, paste("Number of outcomes:                      ", out_char))
    out_char <- sprintf("%.0f", object@n_datapoints)
    out_vec <- c(out_vec, paste("Number of data points:                   ", out_char))
    out_vec <- c(out_vec, paste("Reference treatment:                     ", object@ref_trt))

    out_vec <- c(out_vec, "", "Studies Data Summary:", "---------------------")
    out_char <- sprintf("%.2f", object@n_datapoints/object@n_study)
    out_vec <- c(out_vec, paste("Average number of treatments per study:  ", out_char))
    out_col <- (ncol(object@study_data) - object@n_outcomes + 1):ncol(object@study_data)
    out_dup <- object@study_data[!duplicated(object@study_data[, "studyid"]), ]
    out_char <- sprintf("%.2f", sum(!is.na(out_dup[, out_col]))/object@n_study)
    out_vec <- c(out_vec, paste("Average number of outcomes per study:    ", out_char))

    if (object@type == "binary") {
      out_vec <- c(out_vec, "", "Number of events per outcome:")
      out_sum <- summary(object@study_data[, out_col])
      # temp_summary <- rbind(names(out_sum), out_sum)
      # temp_summary <- format(temp_summary, justify = "centre")
      # temp_summary_mat <- as.vector(apply(X = temp_summary, MARGIN = 1, FUN = function(X) Reduce(function(x, y) paste(x, y), X)))
      # out_vec <- c(out_vec, temp_summary_mat)
      temp_str <- capture.output(out_sum)
      temp_str <- gsub(pattern = "\t", replacement = "  ", x = temp_str)
      if (any(grepl(pattern = "NA", temp_str, fixed = TRUE))) {
        temp_str_index <- grep(pattern = "NA", temp_str, fixed = TRUE)
        temp_str <- temp_str[-temp_str_index]
      }
      out_vec <- c(out_vec, temp_str)

      out_vec <- c(out_vec, "", "Sample sizes per outcome:")
      out_sum <- summary(object@study_data[, (out_col - object@n_outcomes)])
      temp_str <- capture.output(out_sum)
      temp_str <- gsub(pattern = "\t", replacement = "  ", x = temp_str)
      if (any(grepl(pattern = "NA", temp_str, fixed = TRUE))) {
        temp_str_index <- grep(pattern = "NA", temp_str, fixed = TRUE)
        temp_str <- temp_str[-temp_str_index]
      }
      out_vec <- c(out_vec, temp_str)

      out_vec <- c(out_vec, "", "Raw rates (%) per outcome:")
      out_rates <- object@study_data[, out_col]/object@study_data[, (out_col - object@n_outcomes)]*100
      names(out_rates) <- paste0("rr", 1:object@n_outcomes)
      out_sum <- summary(out_rates)
      temp_str <- capture.output(out_sum)
      temp_str <- gsub(pattern = "\t", replacement = "  ", x = temp_str)
      if (any(grepl(pattern = "NA", temp_str, fixed = TRUE))) {
        temp_str_index <- grep(pattern = "NA", temp_str, fixed = TRUE)
        temp_str <- temp_str[-temp_str_index]
      }
      out_vec <- c(out_vec, temp_str)

      out_vec <- c(out_vec, "", "Treatment distribution:")
      out_sum <- t(as.matrix(table(object@study_data$trt)))
      colnames(out_sum) <- paste0("trt ", unique(object@study_data$trt))
      rownames(out_sum) <- "Number of studies"
      temp_str <- capture.output(out_sum)
      temp_str <- gsub(pattern = "\t", replacement = "  ", x = temp_str)
      if (any(grepl(pattern = "NA", temp_str, fixed = TRUE))) {
        temp_str_index <- grep(pattern = "NA", temp_str, fixed = TRUE)
        temp_str <- temp_str[-temp_str_index]
      }
      out_vec <- c(out_vec, temp_str)

      out_vec <- c(out_vec, "", "Study by treatment distribution:")
      out_vec <- c(out_vec, paste0("Treatment  ",
        ifelse(object@n_treatments < 10, "   ", "    "), "Studies"))
      out_tab <- as.vector(
        by(object@study_data$studyid, object@study_data$trt,
          function(x) {
            tb <- table(x)
            toString(names(tb[tb > 0]))})
        )
      for (j in 1:length(out_tab)) {
        out_vec <- c(out_vec, paste0("Treatment",
          ifelse(object@n_treatments > 9, ifelse(j < 10, "   ", "  "), " "), j,
          ":  ", out_tab[j]))#, ifelse(j < length(out_tab), "\n", "")))
      }
    }

    for (i in 1:length(out_vec)) {
      cat(paste0(out_vec[i], "\n"))
    }

    return(invisible(out_vec))
  }
)

#' An S4 class to represent a MCMC run for a model in the \code{netcopula}
#' package.
#'
#' @description
#'   An S4 class to represent the results of fitting the model in the
#'   \code{netcopula} package using a single Markov Chain Monte Carlo chain.
#'
#' @slot x An object of class \code{array}; posterior draws from
#'   the MCMC algorithm for the (adjusted) latent variables.
#' @slot x_unadj An object of class \code{array}; posterior draws from
#'   the MCMC algorithm for the (unadjusted) latent variables.
#' @slot delta An object of class \code{array}; posterior draws from
#'   the MCMC algorithm for the trial-specific log-odds ratios relative to the
#'   baseline.
#' @slot mu An object of class \code{array}; posterior draws from
#'   the MCMC algorithm for the study-specific baseline effects.
#' @slot Gamma An object of class \code{array}; posterior draws from
#'   the MCMC algorithm for the copula correlation matrix.
#' @slot r An object of class \code{array}; posterior draws from
#'   the MCMC algorithm for the copula correlation matrix elements (used only
#'   if Gamma.update is set to 'JointR').
#' @slot d An object of class \code{array}; posterior draws from
#'   the MCMC algorithm for the pooled effects (across trials) of each
#'   treatment relative to the baseline.
#' @slot Sigma An object of class \code{array}; posterior draws from
#'   the MCMC algorithm for the between-study covariance structure.
#' @slot accept An object of class \code{list}; final acceptance rates
#'   for the MCMC algorithm.
#'   simulation.
#' @slot control An object of class \code{list}; list of the control
#'   parameters (number of burnin and sample iterations, etc.). See
#'   \code{\link{nc_control}()} for more information.
#' @slot prior An object of class \code{list}; list of the prior
#'   hyperparameters. See \code{\link{nc_prior}()} for more information.
#' @slot dim An object of class \code{list}; list of dimensions for
#'   the estimated model.
#' @slot data An object of class \code{link{nc_data-class}}; data used in the
#'   analysis.
#' @slot call An object of class \code{call} providing the matched call.
#'
#' @name nc_mcmc-class
#' @rdname nc_mcmc-class
#' @aliases nc_mcmc
#'
#' @exportClass nc_mcmc
#'
#' @examples
#' getSlots("nc_mcmc")
setClass(Class = "nc_mcmc",
  slots = list(
    mu = "array",
    delta = "array",
    d = "array",
    Sigma = "array",
    Gamma = "array",
    r = "array",
    x = "array",
    x_unadj = "array",
    accept = "list",
    control = "list",
    prior = "list",
    dim = "list",
    data = "nc_data",
    call = "call"
  )
)

#' Show an instance of the \code{nc_mcmc} class.
#'
#' @param object An object of class \code{\link{nc_mcmc}}.
#'
#' @author Sergio Venturini \email{sergio.venturini@unito.it}
#'
#' @aliases show,nc_mcmc-method
#' @aliases nc_mcmc-show
#' 
#' @importFrom methods show
#' @exportMethod show
setMethod("show",
  "nc_mcmc",
  function(object) {
    cat("Copula-based multivariate network meta-analysis simulated chain\n")
    cat("Number of studies:", object@dim$n_study, "\n")
    cat("Number of outcomes:", object@dim$n_outcomes, "\n")
    cat("Number of treatments:", object@dim$n_treatments, "\n")
    cat("Reference treatment:", object@dim$ref_trt, "\n")
    cat("\n")
    cat("To get a summary of the object, use the 'summary()' function.")
  }
)

#' Provide a summary of a \code{nc_mcmc} class instance.
#'
#' @param object An object of class \code{\link{nc_mcmc}}.
#' @param include.burnin A length-one logical vector. If \code{TRUE} the
#'  burnin iterations (if available) are included in the summary.
#' @param regex_pars An optional \code{\link[=grep]{regular expression}} to use
#'  for parameter selection. It can be a character vector with multiple strings
#'  or NULL.
#' @param ... Further arguments to pass on (currently ignored).
#'
#' @return A list object containing the summary of the \code{nc_mcmc} object.
#'
#' @author Sergio Venturini \email{sergio.venturini@unito.it}
#'
#' @aliases summary,nc_mcmc-method
#' @aliases nc_mcmc-summary
#' 
#' @exportMethod summary
setMethod("summary",
  "nc_mcmc",
  function(object, include.burnin = FALSE, regex_pars = NULL, ...) {
    if (!is.null(regex_pars)) {
      for (i in 1:length(regex_pars)) {
        if (regex_pars[i] != "mu" && regex_pars[i] != "delta" && regex_pars[i] != "d" && 
          regex_pars[i] != "Sigma" && regex_pars[i] != "Gamma" && regex_pars[i] != "r" &&
          regex_pars[i] != "x" && regex_pars[i] != "x_unadj") {
          stop("the 'regex_pars' elements must be one of 'mu', 'delta', 'd', 'Sigma', 'Gamma', 'r', 'x' or 'x_unadj'.")
        }
      }
    }

    res.coda <- nc_mcmc_to_mcmc(object, include.burnin = include.burnin, verbose = FALSE)
    if (!is.null(regex_pars)) {
      x_idx <- numeric(0)
      for (i in 1:length(regex_pars)) {
        x_idx <- c(x_idx, grep(paste0(regex_pars[i], "["), colnames(res.coda), fixed = TRUE))
      }
    } else {
      x_idx <- grep("d[", colnames(res.coda), fixed = TRUE)
      x_idx <- c(x_idx, grep("mu[", colnames(res.coda), fixed = TRUE))
      x_idx <- c(x_idx, grep("Sigma[", colnames(res.coda), fixed = TRUE))
      x_idx <- c(x_idx, grep("Gamma[", colnames(res.coda), fixed = TRUE))
    }
    res.coda <- res.coda[, x_idx]

    out <- summary(res.coda)

    return(out)
  }
)

#' Provide a graphical summary of a \code{nc_mcmc} class instance.
#'
#' @param x An object of class \code{\link{nc_mcmc}}.
#' @param what A length-one character vector providing the plot type to produce.
#'   Admissible values are those provided by the \pkg{\link{bayesplot}} package,
#'   that is: \code{acf}, \code{areas}, \code{dens}, \code{hex}, \code{hist},
#'   \code{intervals}, \code{neff}, \code{pairs}, \code{parcoord}, \code{recover},
#'   \code{rhat}, \code{scatter}, \code{trace}, \code{violin} or \code{combo}.
#'   In particular, \code{combo} allows to mix different plot types. For more
#'   details see the documentation of the \pkg{\link{bayesplot}} package,
#'   starting from \code{\link[=MCMC-overview]{this overview page}}.
#' @param pars An optional character vector of parameter names. If neither 
#'   \code{pars} nor \code{regex_pars} is specified, the default is to use all
#'   parameters.
#' @param regex_pars An optional \code{\link[=grep]{regular expression}} to use for
#'   parameter selection. Can be specified instead of \code{pars} or in addition to
#'   \code{pars}.
#' @param include.burnin A length-one logical vector. If \code{TRUE} the
#'   burnin iterations (if available) are included in the summary.
#' @param combo A character vector providing the plot types to combine (see
#'   \code{\link[bayesplot]{mcmc_combo}}).
#' @param ... Further arguments to pass on.
#'
#' @return An invisible \code{\link{ggplot}} object providing the required graphical
#'  representation.
#'
#' @author Sergio Venturini \email{sergio.venturini@unito.it}
#'
#' @aliases plot,nc_mcmc-method
#' @aliases nc_mcmc-plot
#' 
#' @exportMethod plot
setMethod("plot",
  signature(x = "nc_mcmc"),
  function(x, what = character(), pars = character(), regex_pars = character(), include.burnin = FALSE,
    combo = NULL, ...) {
    stopifnot(is.character(pars),
              is.character(regex_pars),
              is.character(what))
    
    if (length(what) == 0)
      stop("specify the plot type with the 'what' argument.")

    if (!(what %in% unlist(all_plots_list, use.names = FALSE)))
      stop("the plot type specified is not available.")

    x_mcmc <- nc_mcmc_to_mcmc(x, include.burnin = include.burnin, verbose = FALSE)

    control <- x@control

    ow <- options("warn")
    # options(warn = -1) # suppress all warnings

    if (what %in% acf_plot_list) {
      if (what == "acf") {
        p <- bayesplot::mcmc_acf(x = x_mcmc, pars = pars, regex_pars = regex_pars, ...)
      } else if (what == "acf_bar") {
        p <- bayesplot::mcmc_acf_bar(x = x_mcmc, pars = pars, regex_pars = regex_pars, ...)
      }
    }

    if (what %in% areas_plot_list) {
      if (what == "areas") {
        p <- bayesplot::mcmc_areas(x = x_mcmc, pars = pars, regex_pars = regex_pars, ...)
      } else if (what == "areas_ridges") {
        p <- bayesplot::mcmc_areas_ridges(x = x_mcmc, pars = pars, regex_pars = regex_pars, ...)
      }
    }

    if (what %in% dens_plot_list) {
      if (what == "dens") {
        p <- bayesplot::mcmc_dens(x = x_mcmc, pars = pars, regex_pars = regex_pars, ...)
      } else if (what == "dens_overlay") {
        p <- bayesplot::mcmc_dens_overlay(x = x_mcmc, pars = pars, regex_pars = regex_pars, ...)
      } else if (what == "dens_chains") {
        p <- bayesplot::mcmc_dens_chains(x = x_mcmc, pars = pars, regex_pars = regex_pars, ...)
      }
    }

    if (what %in% hex_plot_list) {
      if (what == "hex") {
        p <- bayesplot::mcmc_hex(x = x_mcmc, pars = pars, regex_pars = regex_pars, ...)
      }
    }

    if (what %in% hist_plot_list) {
      if (what == "hist") {
        p <- bayesplot::mcmc_hist(x = x_mcmc, pars = pars, regex_pars = regex_pars, ...)
      } else if (what == "hist_by_chain") {
        p <- bayesplot::mcmc_hist_by_chain(x = x_mcmc, pars = pars, regex_pars = regex_pars, ...)
      }
    }

    if (what %in% intervals_plot_list) {
      if (what == "intervals") {
        p <- bayesplot::mcmc_intervals(x = x_mcmc, pars = pars, regex_pars = regex_pars, ...)
      }
    }

    if (what %in% neff_plot_list) {
      x_sub <- subset(x, pars = pars, regex_pars = regex_pars)
      nsample <- floor((control[["burnin"]] + control[["nsim"]])/control[["thin"]])
      neff <- coda::effectiveSize(x_sub)
      ratio <- neff/nsample
      if (what == "neff") {
        p <- bayesplot::mcmc_neff(ratio = ratio, ...)
      } else if (what == "neff_hist") {
        p <- bayesplot::mcmc_neff_hist(ratio = ratio, ...)
      }
    }

    if (what %in% pairs_plot_list) {
      if (what == "pairs") {
        p <- bayesplot::mcmc_pairs(x = x_mcmc, pars = pars, regex_pars = regex_pars, ...)
      }
    }

    if (what %in% parcoord_plot_list) {
      if (what == "parcoord") {
        p <- bayesplot::mcmc_parcoord(x = x_mcmc, pars = pars, regex_pars = regex_pars, ...)
      }
    }

    if (what %in% recover_plot_list) {
      x_sub <- subset(x, pars = pars, regex_pars = regex_pars)
      if (what == "recover_hist") {
        p <- bayesplot::mcmc_recover_hist(x = x_sub, ...)
      } else if (what == "recover_intervals") {
        p <- bayesplot::mcmc_recover_intervals(x = x_sub, ...)
      } else if (what == "recover_scatter") {
        p <- bayesplot::mcmc_recover_scatter(x = x_sub, ...)
      }
    }

    if (what %in% rhat_plot_list) {
      x_sub <- subset(x, pars = pars, regex_pars = regex_pars)
      rhat <- coda::gelman.diag(x_sub, multivariate = FALSE)$psrf[, 1]
      if (what == "rhat") {
        p <- bayesplot::mcmc_rhat(rhat = rhat, ...)
      } else if (what == "rhat_hist") {
        p <- bayesplot::mcmc_rhat_hist(rhat = rhat, ...)
      }
    }

    if (what %in% scatter_plot_list) {
      if (what == "scatter") {
        p <- bayesplot::mcmc_scatter(x = x_mcmc, pars = pars, regex_pars = regex_pars, ...)
      }
    }

    if (what %in% trace_plot_list) {
      if (what == "trace") {
        p <- bayesplot::mcmc_trace(x = x_mcmc, pars = pars, regex_pars = regex_pars, ...)
      } else if (what == "trace_highlight") {
        p <- bayesplot::mcmc_trace_highlight(x = x_mcmc, pars = pars, regex_pars = regex_pars, ...)
      }
    }

    if (what %in% violin_plot_list) {
      if (what == "violin") {
        p <- bayesplot::mcmc_violin(x = x_mcmc, pars = pars, regex_pars = regex_pars, ...)
      }
    }

    if (what == "combo") {
      if (!is.null(combo)) {
        p <- bayesplot::mcmc_combo(x = x_mcmc, pars = pars, regex_pars = regex_pars, combo = combo, ...)
      } else {
        stop("to produce an 'mcmc_combo' plot, the 'combo' option must be specified.")
      }
    }

    options(ow) # reset to previous, typically 'warn = 0'

    p
    invisible(p)
  }
)
