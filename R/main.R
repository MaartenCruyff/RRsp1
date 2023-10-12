#' One-sayers model for the ECWM
#'
#' @description Fits the ECWM, the one-sayers model and, if predictors are
#' included in one or both formulas, the
#' logistic regression version of the one-sayers model.
#'
#' @param f0 formula in the form \code{response ~ terms}, where
#' \code{response} is the numeric variable containing the observed responses 1
#' (one 'yes' answer and one 'no' answer) or 2 (two 'yes' or two 'no' answers),
#' and \code{terms} its linear predictor.
#' @param f1 formula in the form \code{~ terms} with the rhs containing
#' the linear predictor of the probability of 'self-protective one-saying'.
#' Defaults to the intercept-only model \code{~ 1}.
#' @param pvar variable containing the probabilities of answering
#' 'yes' to the innocuous (number sequence) question in the subsample.
#' @param data data frame with individual records containing the response
#' variable, the variable \code{pvar}, and the predictors (if any).
#' @return
#' List with the parameter estimates and g.o.f. statistics of the fitted models.
#' @examples
#' ## Intercepts-only models for survey A
#' sp_one(Q4 ~ 1, ~ 1, dplyr::slice(doping[1:4, ], rep(1:4, n)), p4)
#'
#' ## Model for all surveys with survey as predictor of the prevalence and one-saying
#' sp_one(Q4 ~ survey, ~ survey, dplyr::slice(doping, rep(1:20, n)), p4)
#' @importFrom stats get_all_vars model.matrix optim pchisq pt xtabs
#' @importFrom dplyr %>% mutate case_when select count across everything slice
#' @importFrom rlang sym !!
#' @importFrom numDeriv jacobian
#' @export

sp_one <- function(f0,  f1 = ~ 1, data, pvar)
{
  dv   <- f0[[2]]
  y    <- data[[dv]]
  data <- data.frame(data, pyes = data[[substitute(pvar)]], y = y)

  # recode conditional randomization probabilities so that max(pyes)
  # is the probability that the randomized response coincides with the
  # true response.
  pq   <- mutate(data,
                 p = case_when(pyes == min(pyes) & y == 2 ~ min(pyes),
                               pyes == max(pyes) & y == 2 ~ max(pyes),
                               pyes == min(pyes) & y == 1 ~ max(pyes),
                               pyes == max(pyes) & y == 1 ~ min(pyes)),
                 q = 1 - p) %>%
    dplyr::select(p, q) %>%
    as.matrix()

  # prevalence estimate of the standard ecwm without sp-one
  # =======================================================
  m0 <- optim(.2, logl0, pq = pq, method = "Brent", lower = 0, upper = 1, hessian = T)

  out0 <- stats0(data = data, m0 = m0, pyes = pyes, dv = dv)


  # Model matrices and starting value for the fitted model
  # ======================================================
  D0 <- model.matrix(f0, data = data)
  D1 <- model.matrix(f1, data)
  par0 <- rep(0, ncol(D0))
  par1 <- rep(0, ncol(D1))
  par0[1] <- -1
  par1[1] <- -2
  pars    <- c(par0, par1)

  # Fitting the intercepts-only model
  # =======================================================
  m1 <- optim(c(par0[1], par1[1]), logl_reg, pq = pq, y = y,
              D0 = model.matrix(~ 1, data = data),
              D1 = model.matrix(~ 1, data = data),
              hessian = T)

  out1 <- stats1(m0 = m0, m1 = m1, dv = dv, n = nrow(data))

  # Fitting the requested model, if not the intercepts-only model
  # =============================================================
  if (length(pars > 2))
  {
    # 2nd fit with BFGS to avoid NaN's in inverse of hessian
    m2 <- optim(pars, logl_reg, pq = pq, y = y, D0 = D0, D1 = D1, hessian = T)
    m2 <- optim(m2$par, logl_reg, pq = pq, y = y, D0 = D0, D1 = D1, method = "BFGS", hessian = T)

    out2 <- stats2(m1 = m1, m2 = m2, dv = dv, n = nrow(data), pars = pars,
                  names0 = colnames(D0), names1 = colnames(D1))
  }



  cat("Call: \n")
  print(match.call())
  cat("\n")
  cat("M0: ECWM")
  cat("\n")
  cat("============================================================")
  cat("\n")
  cat("Prevalence estimate: \n\n")
  print(out0$coefs)
  cat("\n")
  cat("Observed vs fitted \n\n")
  print(out0$fit)
  cat("\n")
  cat("Fit measures \n\n")
  print(out0$G2)
  cat("\n\n")
  cat("M1: One-sayers model \n")
  cat("============================================================")
  cat("\n")
  cat("Prevalence estimates: \n\n")
  print(out1$coefs)
  cat("\n")
  cat("Fit measures \n\n")
  print(out1$LR)
  cat("\n\n")


  if (length(pars) > 2)
  {
    cat("Logistic regression model \n")
    cat("============================================================")
    cat("\n")
    cat("Logistic parameter estimates \n\n")
    print(out2$coefs)
    cat("\n")
    cat("Fit measures \n\n")
    print(out2$LR)
    cat("\n")
    cat("============================================================")


  }

  invisible(list(M0 = out0,
                 M1 = out1,
                 if (length(pars) > 2) M2 = out2))

}


#' Power curves for one-sayers model
#'
#' @description Renders a plot with the power curves of \eqn{\hat\pi} and
#' \eqn{\hat\theta} given the null hypotheses \eqn{H_0:\pi=0} and \eqn{H_0:\theta=0}.
#'
#' @param pi prevalence \eqn{\pi} of the sensitive attribute.
#' @param theta prevalence \eqn{\theta} of one-saying.
#' @param p probability of answering the innocuous question with 'yes' in one
#' sub-sample and with 'no' in the other.
#' @param n vector with sample sizes for which the power is to be computed.
#' @param alpha significance level,
#' @examples
#' plot_power(pi = .2, theta = .1, p = .2)
#' @importFrom stats pnorm qnorm
#' @importFrom ggplot2 ggplot aes geom_line scale_color_discrete
#' scale_y_continuous theme_minimal
#' @importFrom tidyr pivot_longer
#' @export


plot_power <- function(pi, theta, p, n = seq(10, 1000, by = 10), alpha = 0.05) {

  power      <- matrix(0, length(n), 3, F, list(NULL, c("n", "pi", "theta")))
  power[, 1] <- n

  z <- qnorm(1 - alpha)

  for (j in 1:length(n))
  {
    se_H0_pi    <- se(n = n[j], p = p, pi = 0, theta = theta)
    se_H1_pi    <- se(n = n[j], p = p, pi = pi, theta = theta)
    se_H0_theta <- se(n = n[j], p = p, pi = pi, theta = 0)
    se_H1_theta <- se(n = n[j], p = p, pi = pi, theta = theta)
    power[j, 2] <- pnorm((pi    - z * se_H0_pi[1])    / se_H1_pi[1])
    power[j, 3] <- pnorm((theta - z * se_H0_theta[2]) / se_H1_theta[2])

  }
  pivot_longer(as.data.frame(power), 2:3, names_to = "par", values_to = "power") %>%
    ggplot(aes(n, power, col = par)) +
    geom_line() +
    scale_y_continuous(limits = c(0, 1),
                       breaks = c(0, 0.2, 0.4, 0.6, 0.8, 1)) +
    scale_color_discrete(labels = c(expression(pi),
                                    expression(theta))) +
    theme_minimal()
}

