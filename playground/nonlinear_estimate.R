#' @title The Zhang (2008) model
#'
#' @param data input data.frame
#' @param counts observed counts
#' @param pop_size population size
#' @param auxiliary auxiliary variable (single)
#' @param family what to assume: gaussian, poisson, NB
#' @param control control parameters (TODO)
#'
#' @export
nonlinear_estimate <- function(data,
                               counts,
                               pop_size,
                               auxiliary,
                               family,
                               control) {


  loglik <- function(alpha, beta, phi, m, N, n, X, Z) {
    phi_exp <- exp(-phi)
    mu <- (N^(X %*% alpha)) * ((n/N)^(Z %*% beta))
    ll <- sum(m * log(mu) - (m + phi_exp) * log(mu + phi_exp) + phi_exp * log(phi_exp) -
          lgamma(phi_exp) - lgamma(m + 1) + lgamma(m + phi_exp))
    return(-ll)
  }

  grad_analytical <- function(alpha, beta, phi, m, N, n, X, Z) {
    phi_exp <- exp(-phi)
    mu <- (N^(X %*% alpha)) * ((n/N)^(Z %*% beta))
    dmu <- m/mu - (phi_exp + m)/(mu + phi_exp)
    dalpha <- dmu * mu * log(N)
    dbeta <- dmu * mu * log(n/N)
    dphi <- sum(digamma(phi_exp + m) - log(mu + phi_exp) + log(phi_exp) -
                  digamma(phi_exp) + 1 - (phi_exp + m)/(phi_exp + mu))
    dphi <- dphi * (-phi_exp) # neg log link
    grad <- c(t(X) %*% dalpha, t(Z) %*% dbeta, dphi)
    return(-grad)
  }

  hess_analytical <- function(alpha, beta, phi, m, N, n, X, Z) {
    phi_exp <- exp(-phi)
    mu <- (N^(X %*% alpha)) * ((n/N)^(Z %*% beta))
    dmu <- m/mu - (phi_exp + m)/(mu + phi_exp)
    dmu_2 <- (phi_exp + m)/(mu + phi_exp)^2 - m/mu^2
    dmudphi <- -(mu - m)/(mu + phi_exp)^2
    # neg log link
    dmudphi <- dmudphi * (-phi_exp)
    dphi_2 <- sum(trigamma(phi_exp + m) - trigamma(phi_exp) - 2/(phi_exp + mu) +
                    (phi_exp + m)/(phi_exp + mu)^2 + 1/phi_exp)
    # log link
    dphi_2 <- dphi_2 * phi_exp^2
    dphi_2 <- dphi_2 + sum(digamma(phi_exp + m) - log(mu + phi_exp) + log(phi_exp) -
                             digamma(phi_exp) + 1 - (phi_exp + m)/(phi_exp + mu)) * phi_exp
    dalpha_2 <- dmu_2 * (mu * log(N))^2 + dmu * mu * (log(N)^2)
    dbeta_2 <- dmu_2 * (mu * log(n/N))^2 + dmu * mu * (log(n/N)^2)
    dalphadbeta <- dmu_2 * mu * log(N) * mu * log(n/N) + dmu * mu * log(N) * log(n/N)
    dalphadphi <- dmudphi * mu * log(N)
    dbetadphi <- dmudphi * mu * log(n/N)

    hess <- rbind(
      cbind(t(X) %*% (dalpha_2 * X), t(X) %*% (dalphadbeta * Z), t(X) %*% dalphadphi),
      cbind(t(t(X) %*% (dalphadbeta * Z)), t(Z) %*% (dbeta_2 * Z), t(Z) %*% dbetadphi),
      cbind(t(t(X) %*% dalphadphi), t(t(Z) %*% dbetadphi), dphi_2)
    )
    return(-hess)
  }

}
