#' Fit the model via the variational Bayes
#' 
#' Obtain the posterior distribution
#' 
#' @export

fit.vbsrm2 <- function(time = NULL, fault = NULL, type = NULL,
                       te = NULL, data = data.frame(), alpha = NA,
                       prior.omega = list(shape=0, rate=1),
                       prior.beta = list(shape=0, rate=1),
                       control = list(), ...) {
  data <- Rsrat::faultdata(time, fault, type, te, data)
  con <- vbsrm2.options()
  nmsC <- names(con)
  con[(namc <- names(control))] <- control
  if (length(noNms <- namc[!namc %in% nmsC]))
    warning("unknown names in control: ", paste(noNms, collapse = ", "))

  model <- vb2.gamma$new(alpha = alpha)
  model$set_prior(c(prior.omega$shape, prior.beta$shape), c(prior.omega$rate, prior.beta$rate))
  model$set_data(data)
  tres <- system.time(result <- model$vbem(con))
  c(alpha=model$alpha, result)
}

#' Options for vbsrm2
#'
#' Generate a list of option values.
#'
#' @return A list of options.
#' @export

vbsrm2.options <- function() {
  list(maxiter = 10000,
    reltol = sqrt(.Machine$double.eps),
    abstol = 1.0e+200,
    eps = sqrt(.Machine$double.eps),
    trace = FALSE,
    rmax = 1000)
}

