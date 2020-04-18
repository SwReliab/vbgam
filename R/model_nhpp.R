#' Class for NHPP-based software reliability model
#'
#' @docType class
#' @name NHPP
#' @return Object of \code{\link{R6Class}} with methods for NHPP-based software reliability model.
#' @format \code{\link{R6Class}} object.
#' @field name A character string for the name of model.
#' @field params A numeric vector for the model parameters.
#' @field df An integer for the degrees of freedom of the model.
#' @field data Data to esimate parameters.
#'
#' @section Methods:
#' \describe{
#'   \item{\code{print()}}{This method prints model parameters.}
#'   \item{\code{omega()}}{This method returns the number of total faults.}
#'   \item{\code{mvf(t)}}{This method returns the mean value function at time t.}
#'   \item{\code{dmvf(t)}}{This method returns the mean value function on discrete time domain.}
#'   \item{\code{inv_mvf(x)}}{This method returns the time at which the mean value function attains x.}
#'   \item{\code{intensity(t)}}{This method returns the intensity function at time t.}
#'   \item{\code{reliab(t, s)}}{This method returns the software reliability at time t from the orign s.}
#'   \item{\code{residual(t)}}{This method returns the expected residual number of faults at time t.}
#'   \item{\code{ffp(t)}}{This method returns the fault-free probability at time t.}
#'   \item{\code{imtbf(t)}}{This method returns the instantaneous MTBF at time t.}
#'   \item{\code{cmtbf(t)}}{This method returns the cumulative MTBF at time t.}
#'   \item{\code{median(s, p = 0.5)}}{This method returns the time at which the software reliability attains the proability p from the orign s.}
#'   \item{\code{init_params(data)}}{This method changes the model parameters based on a given data. This is used to set the initial value for the fitting algorithm.}
#'   \item{\code{set_params(params)}}{This method sets the model parameters.}
#'   \item{\code{em(params, data)}}{This method returns a list with an updated parameter vector (param),
#'          absolute difference of parameter vector (pdiff),
#'          log-likelihood function for a given parameter vector (llf),
#'          the number of total faults (total) via EM algorithm for a given data. \emph{divide} in GammaSRM is the number of integration points.}
#'   \item{\code{llf(data)}}{This method returns the log-likelihood function for a given data.}
#' }
#' @seealso \code{\link{srm}}
NULL

#' @rdname NHPP
vb2.gamma <- R6::R6Class("vb2.gamma",
  # private = list(
  #   prior = NA,
  #   posterior = NA
  # ),
  public = list(
    alpha = NA,
    prior = NA,
    posterior = NA,
    data = NA,
    # print = function(digits = max(3, getOption("digits") - 3), ...) {
    #   cat(gettextf("Model name: %s\n", self$name))
    #   print.default(format(self$params, digits = digits), print.gap = 2, quote = FALSE)
    # },
    initialize = function(alpha = NA) {
      self$alpha <- alpha
      self$prior <- list(omega=c(1,1), beta=c(1,1))
    },
    set_alpha = function(alpha) {
      self$alpha <- alpha
    },
    set_prior = function(omega, beta) {
      self$prior <- list(omega=omega, beta=beta)
    },
    set_data = function(data) {
      self$data <- data
    },
    set_posterior = function(posterior) {
      self$posterior <- posterior
    },
    set_vfe = function(vfe) {
      self$vfe <- vfe
    },
    get_residual = function() {
      n <- self$posterior$n - self$data$total
      p <- self$posterior$mix
      list(prob = data.frame(n=n, p=p),
        mean = sum(n * p),
        ffp = p[1])
    },
    vbem = function(options) {
      if (is.na(self$alpha)) {
        optim.fn <- function(alpha) {
          -vb2fit(alpha=exp(alpha), prior=self$prior, data=self$data, options=options)$vfe
        }
        result <- optim(par=1, fn=optim.fn, control=options)
        if (result$convergence == 0) {
          self$alpha <- exp(result$par)
        }
      }
      vb2fit(alpha=self$alpha, prior=self$prior, data=self$data, options=options)
    }
  )
)
