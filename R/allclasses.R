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
#' @exportClass nc_data
#'
#' @examples
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
#' @param \ldots Optional arguments to the function in its next call.
#'
#' @aliases initialize,nc_data-method
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
#' @export
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
#' @aliases summary,nc_data-method
#' 
#' @export
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
      out_tab <- as.vector(by(object@study_data$studyid, object@study_data$trt, function(x) toString(names(table(x)))))
      out_tab <- data.frame(Treatment = paste0("Treatment ", unique(object@study_data$trt)), Studies = out_tab)
      temp_str <- capture.output(out_tab)
      out_vec <- c(out_vec, temp_str)
    }

    outTab <- as.table(matrix(out_vec,ncol = 1))
    dimnames(outTab) <- list(rep("", length(out_vec)), "")

    return(outTab)
  }
)

# setMethod("plot", signature(x = "nc_data"),
#   function(x, y, ...) {
#     survey.Data <- x@surveyData
#     if (length(survey.Data@nAnimalVec > 0)) {
#       par(mfrow = c(2, 2))      
#       ## Mean number of animals per herd:
#       plot(x@herdSensVec, x@nAnimalsMeanVec/x@nHerdsVec, type = "l",
#               xlab = "Herd sensitivity", ylab = "Mean no. of animals per herd to be tested")
#       ## Number of herds to be tested:
#       plot(x@herdSensVec, x@nHerdsVec, type = "l",
#               xlab = "Herd sensitivity", ylab = "No. of herds to be tested")
#       ## Total number of animals to be tested:
#       plot(x@herdSensVec, x@nAnimalsMeanVec, type = "l",
#               xlab = "Herd sensitivity", ylab = "Expected total no. of animals to be tested")
#       ## Expected cost:
#       plot(x@herdSensVec, x@expectedCostVec, type = "l",
#               xlab = "Herd sensitivity", ylab = "Expected cost")           
#       ## Titel:
#       par(oma = c(2,1,3,1)) 
#       title("Analysis individual sampling", outer = TRUE)                  
#     } else {
#       cat("Object of class 'IndSamplingSummary' contains no data.\n")
#     }
#   }
# )

setOldClass("mcmc")

#' An S4 class to represent a MCMC run for a model in the \code{netcopula}
#' package.
#'
#' @slot x BLA BLA BLA...
#' @slot delta BLA BLA BLA...
#' @slot mu BLA BLA BLA...
#' @slot Gamma BLA BLA BLA...
#' @slot d BLA BLA BLA...
#' @slot Sigma BLA BLA BLA...
#' @slot accept BLA BLA BLA...
#' @slot dens BLA BLA BLA...
#' @slot control BLA BLA BLA...
#' @slot dim BLA BLA BLA...
#' @slot call BLA BLA BLA...
#'
#' @name nc_mcmc-class
#' @rdname nc_mcmc-class
#' @aliases nc_mcmc-class
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
    Sigma = "mcmc",
    Gamma = "array",
    r = "array",
    x = "array",
    x_unadj = "array",
    # a = "array",
    # b = "array",
    # D = "array",
    # S_q = "array",
    # Sigma_q_prop = "array",
    # rho_mu = "array",
    # cov_mu = "array",
    # ar_mu_vec = "array",
    # mu_rate = "array",
    # mu_prop = "array",
    # tp_mu = "array",
    # delta_prop = "array",
    accept = "list",
    dens = "list",
    control = "list",
    dim = "list",
    call = "call",
    data = "nc_data"
  ),
  # http://stackoverflow.com/questions/13841400/use-s3-virtual-class-as-slot-of-an-s4-class-got-error-got-class-s4-should-b
  prototype = prototype(
    x = structure(list(), class = "array"),
    x_unadj = structure(list(), class = "array"),
    # a = structure(list(), class = "array"),
    # b = structure(list(), class = "array"),
    delta = structure(list(), class = "array"),
    mu = structure(list(), class = "array"),
    Gamma = structure(list(), class = "array"),
    r = structure(list(), class = "array"),
    d = structure(list(), class = "array"),
    Sigma = structure(list(), class = "mcmc")
    # D = structure(list(), class = "array"),
    # S_q = structure(list(), class = "array"),
    # Sigma_q_prop = structure(list(), class = "array")
    # rho_mu = structure(list(), class = "array"),
    # cov_mu = structure(list(), class = "array"),
    # ar_mu_vec = structure(list(), class = "array"),
    # mu_rate = structure(list(), class = "array"),
    # mu_prop = structure(list(), class = "array"),
    # tp_mu = structure(list(), class = "array"),
    # delta_prop = structure(list(), class = "array")
  )
)

#' @describeIn nc_mcmc Provide a summary of a \code{\link{nc_mcmc}} class
#' instance.
#'
#' @aliases summary,nc_mcmc-method
#' 
#' @export
setMethod("summary", signature(object = "nc_mcmc"),
  function(object, digits_summary = 1, latent_var = FALSE, add_corr = FALSE, print = TRUE, ...) {
    moreargs <- list(...)
    if (is.null(object@control$nthin)) {
      nthin <- moreargs[[match("nthin", names(moreargs))]]
      data <- moreargs[[match("data", names(moreargs))]]
      ref_trt <- data@ref_trt
      res <- as.mcmc(object, latent_var = latent_var, nthin = nthin, data = data)
      moreargs <- moreargs[-c(match("nthin", names(moreargs)), match("data", names(moreargs)))]
    } else {
      res <- as.mcmc(object, latent_var = latent_var)
    }
    if (add_corr) {
      res <- add_corr(res)
    }

    parameter.names <- dimnames(res)[[2]]
    n.keep <- nrow(res)
    n.parameters <- length(parameter.names)
    n.chains <- 1
    sims.array <- array(NA, c(n.keep, n.chains, n.parameters))
    for (i in 1:n.chains) {
      sims.array[, i, ] <- as.matrix(res)
    }
    dimnames(sims.array) <- list(NULL, NULL, parameter.names)

    summary <- build_output(sims.array)

    if (print) {
      print(round(summary, digits_summary), ... = moreargs)
    }
    
    invisible(summary)
  }
)

#' @describeIn nc_mcmc Provide a graphical summary of a \code{\link{nc_mcmc}}
#' class instance.
#'
#' @aliases plot,nc_mcmc-method
#' 
#' @export
setMethod("plot", signature(x = "nc_mcmc"),
  function(x, y, latent_var = FALSE, add_corr = TRUE, to_plot, type = "trace", reorder = FALSE, ...) {
    to_plot_list <- c("Gamma", "Sigma_M", "Corr_M", "mu", "delta", "d", "x", "x_unadj", "S_q", "Sigma_q_prop")
    if (!(to_plot %in% to_plot_list)) {
      stop(paste0("no chain available for parameters '", to_plot,"'."))
    }
    if (!latent_var & (to_plot == "x")) {
      stop("to inspect the latent variables chain set the 'latent_var' argument to 'TRUE'.")
    }
    if (!add_corr & (to_plot == "Corr_M")) {
      stop("to inspect the common between-study outcome correlations chain set the 'add_corr' argument to 'TRUE'.")
    }

    res <- x
    res_mcmc <- as.mcmc(res, latent_var = latent_var)
    if (add_corr) {
      res_mcmc <- add_corr(res_mcmc)
    }

    if (type == "caterpillar") {
      out <- summary(res, digits_summary = 5, latent_var = latent_var, add_corr = add_corr, print = FALSE)
      to_keep <- substr(rownames(out), 1, nchar(to_plot) + 1) == paste0(to_plot, "[")
      out_tmp <- out[to_keep, , drop = FALSE]
      caterplot(res_mcmc, to_plot, reorder = reorder, quantiles = list(outer = c(0.025, 0.975), inner = c(0.1, 0.9)))
      if (reorder) {
        out_i <- sort.int(out_tmp[, 5], decreasing = TRUE, index.return = TRUE)$ix
      } else {
        out_i <- 1:nrow(out_tmp)
      }
      caterpoints(out_tmp[out_i, 1], pch = "x", col = "darkorange")
      title(paste0("parameters ", to_plot), cex.main = 1.25)
    } else if (type == "trace") {
      traplot(res_mcmc, parms = to_plot, mar = c(1.0, 1.5, .75, .15) + .1, col = NULL, lty = 1, plot.title = NULL, main = NULL, greek = FALSE, style = "plain", cex.main = .85, cex.axis = .7, tcl = -.3, xaxt = "n")
    }
  }
)
