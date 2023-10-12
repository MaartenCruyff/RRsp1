
### logistic function

logistic <- function(x) exp(x) / ( 1 + exp(x))


# log likelihood for the intercepts-only model

logl0 <- function(ptrue, pq)
{
  -sum(log(pq %*% c(ptrue[1], 1 - ptrue[1])))
}

### log likelihood for regression models

logl_reg <- function(pars, pq, y, D0, D1)
{
  pi    <- logistic(D0 %*% pars[1:ncol(D0)])
  ptrue <- cbind(c(pi), 1 - c(pi))
  theta <- c(logistic(D1 %*% pars[-(1:ncol(D0))]))

  spq          <- (1 - theta) * pq
  spq[y == 1]  <- spq[y == 1] + theta[y == 1]

  -sum(log(rowSums(spq * ptrue)))
}

### summary statistics model m0

stats0 <- function(data, m0, pyes, dv)
{
  ptrue <- m0$par
  fit   <- count(data, pyes, !!sym(dv))
  p     <- fit[1, 1]
  est1  <- matrix(c(1 - p, p, p, 1 - p), 2) %*% c(ptrue, 1 - ptrue) * sum(fit$n[1:2])
  est2  <- matrix(c(p, 1 - p, 1 - p, p), 2) %*% c(ptrue, 1 - ptrue) * sum(fit$n[3:4])
  fit   <- data.frame(fit, nfitted = c(est1, est2))

  G2    <- data.frame(loglike = -m0$val,
                      G2     = 2 * sum(fit$n * log(fit$n / fit$nfitted)),
                      df     = 1, row.names = "stats") %>%
    mutate(p  = 1 - pchisq(G2, 1),
           across(c(loglike, G2), ~round(.x, 2)),
           p  = round(p, 4))

  coefs <- data.frame(est = m0$par,
                      se  = sqrt(1/m0$hessian),
                      row.names = paste(dv)) %>%
    mutate(t = est / se,
           p = pt(-abs(t), nrow(data) - 1),
           across(everything(), ~ round(.x, 3)))

  list(G2 = G2, coefs = coefs)
}


### summary statistics model m1

stats1 <- function(m0, m1, dv, n)
{
  par   <- m1$par
  vcv   <- solve(m1$hessian)
  var   <- t(jacobian(logistic, par)) %*% vcv %*% jacobian(logistic, par)
  se    <- sqrt(diag(var))

  coefs <- data.frame(est = logistic(par),
                      se  = se,
                      t   = logistic(par) / se,
                      p   = pt(q = -abs(logistic(par) / se), df = n),
                      row.names = c(paste(dv), "sp-one")) %>%
    mutate(across(everything(), ~ round(.x, 3)))

  LR    <- data.frame (loglike = -m1$val,
                       LR      = 2 * (m0$val - m1$val),
                       df      = 1,
                       row.names = "stats") %>%
    mutate(p  = 1 - pchisq(LR, 1),
           across(c(loglike, LR), ~round(.x, 2)),
           p  = round(p, 4))
  list(LR = LR, coefs = coefs)
}

### summary statistics model m2

stats2 <- function(m1, m2, dv, n, pars, names0, names1)
{

  par <- m2$par
  se  <- sqrt(diag(solve(m2$hessian)))
  coefs  <- data.frame(est = par,
                       se  = se,
                       t   = par / se,
                       p   = pt(-abs(par / se), df  = n - length(par)),
                       row.names = c(paste(dv, names0), paste("sp-one", names1))) %>%
    mutate(across(everything(), ~round(.x, 3)))

  LR <- data.frame(loglike = -m2$val,
                    LR      = 2 * (m1$val - m2$val),
                    df      = length(pars) - 2,
                    row.names = "stats") %>%
    mutate(p  = 1 - pchisq(LR, length(pars) - 2),
           across(c(loglike, LR), ~round(.x, 2)),
           p  = round(p, 4))

  list(LR = LR, coefs = coefs)
}

## Function that returns the standard errors of the estimators of pi and theta

se <- function(n, p, theta, pi){


  # transition matrix ECWM
  P     <-  matrix(c(p,    1 - p,
                     1 - p,    p,
                     1 - p,    p,
                     p,    1 - p),
                   4, 2)

  # transition matrix one-sayers model
  Q            <- (1 - theta) * P
  Q[c(1, 3), ] <- Q[c(1, 3), ] + theta

  # conditional randomized response probabilities

  cp <- Q %*% c(pi, 1 - pi)

  # standard error pihat
  var1  <- cp[1] * cp[2] * (cp[4] / ((2 * p - 1) * (cp[2] + cp[4])^2))^2 / (n / 2)
  var2  <- cp[3] * cp[4] * (cp[2] / ((2 * p - 1) * (cp[2] + cp[4])^2))^2 / (n / 2)
  se_pi <- sqrt(var1 + var2)

  # standard error thetahat
  se_theta <- sqrt((cp[1] * cp[2] + cp[3] * cp[4]) / (n / 2))

  c(se_pi, se_theta)

}
