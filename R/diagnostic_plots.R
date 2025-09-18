#' @importFrom ggplot2 ggplot aes geom_point geom_smooth geom_hline geom_vline labs theme_classic theme element_rect stat_qq stat_qq_line coord_cartesian geom_line
#' @importFrom ggrepel  geom_text_repel
#' @importFrom gridExtra grid.arrange
#'
#' @title Plot diagnostics for a \code{hidden} object
#'
#' @description Diagnostic plotting method for objects of class \code{hidden}.
#' Depending on the estimation method used in \code{estimate_hidden_pop()} function, different sets of diagnostic
#' plots are displayed.
#'
#' @param x an object of class \code{hidden} returned by the main estimation function \code{estimate_hidden_pop()}
#' @param which optional integer indication which plot should be displayed. If \code{NULL}, the default set of plots is shown.
#' @param label logical. If \code{TRUE} adds observation labels (country names) to the plots (only if country information is in the object).
#'
#' @details
#' The function displays model diagnostics depending on the estimation method used.
#' The numbers in parentheses indicate the value to use in the \code{which} argument.
#' \itemize{
#'  \item \code{'ols'}: Residuals vs fitted (1), Q-Q plot (2), scale-location plot (3), residuals vs leverage (4), observed vs fitted (5).
#'  \item \code{'nls'}: Residuals vs fitted (1), Q-Q plot (2), observed vs fitted (3).
#'  \item \code{'glm'}: The same set of plots as in \code{ols}, but based on deviance/standardized deviance residuals.
#'  \item \code{'mle'}: √Observed vs √fitted (1), Q-Q plot for Anscombe residuals (2), absolute Anscombe residuals vs √Fitted (3), Anscombe residuals vs √Fitted (4).
#' }
#' If \code{which} is specified, only the corresponding plot is shown.
#'
#'
#' @examples
#' \dontrun{
#' # Load Polish irregular migration data
#' data(foreigners_pl)
#'
#' model_data <- subset(foreigners_pl, year == 2018 & half == 1)
#'
#' # Basic Zhang model estimation using GLM (recommended)
#' result_zhang <- estimate_hidden_pop(
#'   data = model_data,
#'   observed = ~ border,
#'   auxiliary = ~ police,
#'   reference_pop = ~ pesel,
#'   method = "mle",
#'   family = "poisson"
#' )
#'
#' plot(result_zhang)
#'
#' }
#'@export
plot.hidden <- function(x, which = NULL, label = FALSE){

  if (x$method == 'ols') {

    df <- data.frame(fitted = x$fitted,
                     observed = x$m,
                     residuals = x$residuals,
                     resid_stand = x$resid_stand,
                     cooks = x$cooks,
                     leverage = x$leverage)
    if (!is.null(x$countries)) {
      df$nationality <- x$countries
    }

    # 1 residuals vs fitted
    p1 <- ggplot(df, aes(x = fitted, y = residuals)) +
      geom_point(shape = 1) +
      geom_hline(yintercept = 0, linetype = 'dotted') +
      labs(title = 'Residuals vs Fitted values', x = 'Fitted values', y = 'Residuals') +
      theme_classic() +
      theme(panel.border = element_rect(fill = NA))
    if (label && !is.null(x$countries)) p1 <- p1 + geom_text_repel(data = df, aes(label = nationality))

    # 2 q-q residuals
    p2 <- ggplot(df, aes(sample = residuals)) +
      stat_qq_line(lwd = 0.3) +
      stat_qq(shape=1) +
      labs(title = 'Q-Q Residuals', x = 'Theoretical Quantiles', y = 'Standardized residuals') +
      theme_classic() +
      theme(panel.border = element_rect(fill = NA))

    # 3 scale location plot
    p3 <- ggplot(df, aes(x = fitted, y = sqrt(abs(resid_stand)))) +
      geom_point(shape = 1) +
      geom_smooth(se = FALSE, color = 'red', lwd = 0.5) +
      labs(title = 'Scale-Location', x = 'Fitted values', y = '√|Standardized residuals|') +
      theme_classic() +
      theme(panel.border = element_rect(fill = NA))
    if (label && !is.null(x$countries)) p3 <- p3 + geom_text_repel(data = df, aes(label = nationality))

    # 4 residuals vs leverage
    h_seq <- seq(min(df$leverage), max(df$leverage) + 0.025, length.out = 100)

    cook_line <- function(h, q) {
      sqrt(q * (1 - h) / h)
    }

    cook_df <- data.frame(
      leverage = rep(h_seq, 4),
      resid = c(cook_line(h_seq, 0.5),
                -cook_line(h_seq, 0.5),
                cook_line(h_seq, 1),
                -cook_line(h_seq, 1)),
      group = factor(rep(c('+0.5', '-0.5', '+1', '-1'), each = length(h_seq)))
    )

    p4 <- ggplot(df, aes(x = leverage, y = resid_stand)) +
      geom_smooth(se = FALSE, color = 'red', lwd = 0.5) +
      geom_point(shape = 1) +
      geom_hline(yintercept = 0, linetype = 'dotted') +
      geom_vline(xintercept = 0, linetype = 'dotted') +
      geom_line(data = cook_df, aes(x = leverage, y = resid, group = group),
                color = 'gray', linetype = 'dashed', lwd = 0.4) +
      labs(title = 'Residuals vs Leverage',
           x = 'Leverage', y = 'Standardized residuals') +
      coord_cartesian(xlim = c(-0.001, max(df$leverage) + 0.001), ylim = c(min(df$resid_stand) - 1, max(df$resid_stand) + 1)) +
      theme_classic() +
      theme(panel.border = element_rect(fill = NA))
    if (label && !is.null(x$countries)) p4 <- p4 + geom_text_repel(data = df, aes(label = nationality))

    # 5 observed vs fitted
    p5 <- ggplot(df, aes(x = fitted, y = observed)) +
      geom_smooth(method = 'lm', se = FALSE, color = 'black',lty = 'dashed', lwd = 0.5) +
      geom_point(shape = 1) +
      labs(title = 'Observed vs Fitted', x = 'Fitted', y = 'Observed') +
      theme_classic() +
      theme(panel.border = element_rect(fill = NA))
    if (label && !is.null(x$countries)) p5 <- p5 + geom_text_repel(data = df, aes(label = nationality))

    if (!is.null(which) && which %in% seq(1:5)){
      plot <- switch(which, '1' = p1, '2' = p2, '3' = p3, '4' = p4, '5' = p5)
      print(plot)
    } else{
      grid.arrange(p1, p2, p3, p4, ncol = 2, nrow = 2)
    }

  } else if (x$method == 'nls') {

    df <- data.frame(fitted = x$fitted,
                     observed = x$m,
                     residuals = x$residuals)
    if (!is.null(x$countries)) {
      df$nationality <- x$countries
    }

    # 1 residuals vs fitted
    p1 <- ggplot(df, aes(x = fitted, y = residuals)) +
      geom_point(shape = 1) +
      geom_hline(yintercept = 0, linetype = 'dotted') +
      labs(title = 'Residuals vs Fitted values', x = 'Fitted values', y = 'Residuals') +
      theme_classic() +
      theme(panel.border = element_rect(fill = NA))
    if (label && !is.null(x$countries)) p1 <- p1 + geom_text_repel(data = df, aes(label = nationality))

    # 2 q-q residuals
    p2 <- ggplot(df, aes(sample = residuals)) +
      stat_qq_line(lwd = 0.3) +
      stat_qq(shape=1) +
      labs(title = 'Q-Q Residuals', x = 'Theoretical Quantiles', y = 'Sample Quantiles') +
      theme_classic() +
      theme(panel.border = element_rect(fill = NA))

    # 3 observed vs fitted
    p3 <- ggplot(df, aes(x = fitted, y = observed)) +
      geom_smooth(method = 'lm', se = FALSE, color = 'black',lty = 'dashed', lwd = 0.5) +
      geom_point(shape = 1) +
      labs(title = 'Observed vs Fitted', x = 'Fitted', y = 'Observed') +
      theme_classic() +
      theme(panel.border = element_rect(fill = NA))
    if (label && !is.null(x$countries)) p3 <- p3 + geom_text_repel(data = df, aes(label = nationality))

    if (!is.null(which) && which %in% seq(1:3)){
      plot <- switch(which, '1' = p1, '2' = p2, '3' = p3)
      print(plot)
    } else{
      grid.arrange(p1, p2, ncol = 2)
    }

  } else if (x$method == 'glm - Poisson') {

    df <- data.frame(fitted = x$fitted,
                     observed = x$m,
                     residuals = x$residuals,     # deviance residuals
                     resid_stand = x$resid_stand, # standardized deviance residuals
                     cooks = x$cooks,
                     leverage = x$leverage)
                     # pearson = (results$m - results$fitted)/sqrt(results$fitted),
                     # anscombe = (3*results$m^(2/3)-3*results$fitted^(2/3))/(2 * results$fitted^(1/6))
    if (!is.null(x$countries)) {
      df$nationality <- x$countries
    }

    # 1 residuals vs fitted
    p1 <- ggplot(df, aes(x = fitted, y = residuals)) +
      geom_point(shape = 1) +
      geom_hline(yintercept = 0, linetype = 'dotted') +
      labs(title = 'Residuals vs Fitted values', x = 'Fitted values', y = 'Residuals') +
      theme_classic() +
      theme(panel.border = element_rect(fill = NA))
    if (label && !is.null(x$countries)) p1 <- p1 + geom_text_repel(data = df, aes(label = nationality))

    # 2 q-q residuals
    p2 <- ggplot(df, aes(sample = resid_stand)) +
      stat_qq_line(lwd = 0.3) +
      stat_qq(shape=1) +
      labs(title = 'Q-Q Residuals', x = 'Theoretical Quantiles', y = 'Standardized residuals') +
      theme_classic() +
      theme(panel.border = element_rect(fill = NA))

    # 3 scale location plot
    p3 <- ggplot(df, aes(x = fitted, y = sqrt(abs(resid_stand)))) +
      geom_point(shape = 1) +
      geom_smooth(se = FALSE, color = 'red', lwd = 0.5) +
      labs(title = 'Scale-Location', x = 'Fitted values', y = '√|Standardized residuals|') +
      theme_classic() +
      theme(panel.border = element_rect(fill = NA))
    if (label && !is.null(x$countries)) p3 <- p3 + geom_text_repel(data = df, aes(label = nationality))

    # 4 residuals vs leverage
    h_seq <- seq(min(df$leverage), max(df$leverage) + 0.025, length.out = 100)

    cook_line <- function(h, q) {
      sqrt(q * (1 - h) / h)
    }

    cook_df <- data.frame(
      leverage = rep(h_seq, 4),
      resid = c(cook_line(h_seq, 0.5),
                -cook_line(h_seq, 0.5),
                cook_line(h_seq, 1),
                -cook_line(h_seq, 1)),
      group = factor(rep(c('+0.5', '-0.5', '+1', '-1'), each = length(h_seq)))
    )

    p4 <- ggplot(df, aes(x = leverage, y = resid_stand)) +
      geom_smooth(se = FALSE, color = 'red', lwd = 0.5) +
      geom_point(shape = 1) +
      geom_hline(yintercept = 0, linetype = 'dotted') +
      geom_vline(xintercept = 0, linetype = 'dotted') +
      geom_line(data = cook_df, aes(x = leverage, y = resid, group = group),
                color = 'gray', linetype = 'dashed', lwd = 0.4) +
      labs(title = 'Residuals vs Leverage',
           x = 'Leverage', y = 'Standardized residuals') +
      coord_cartesian(xlim = c(-0.001, max(df$leverage) + 0.001), ylim = c(min(df$resid_stand) - 1, max(df$resid_stand) + 1)) +
      theme_classic() +
      theme(panel.border = element_rect(fill = NA))
    if (label && !is.null(x$countries)) p4 <- p4 + geom_text_repel(data = df, aes(label = nationality))

    # 5 observed vs fitted
    p5 <- ggplot(df, aes(x = fitted, y = observed)) +
      geom_smooth(method = 'lm', se = FALSE, color = 'black',lty = 'dashed', lwd = 0.5) +
      geom_point(shape = 1) +
      labs(title = 'Observed vs Fitted', x = 'Fitted', y = 'Observed') +
      theme_classic() +
      theme(panel.border = element_rect(fill = NA))
    if (label && !is.null(x$countries)) p5 <- p5 + geom_text_repel(data = df, aes(label = nationality))

    if (!is.null(which) && which %in% seq(1:5)){
      plot <- switch(which, '1' = p1, '2' = p2, '3' = p3, '4' = p4, '5' = p5)
      print(plot)
    } else{
      grid.arrange(p1, p2, p3, p4, ncol = 2, nrow = 2)
    }

  } else if (x$method == 'mle') {

    phi <- x$coefficients[length(x$coefficients)]
    df <- data.frame(fitted = x$fitted,
                     observed = x$m,
                     residuals = x$residuals,
                     anscombe = (3*phi*(1 + x$m/phi)^(2/3) - (1 + x$fitted/phi)^(2/3) + 3*(x$m^(2/3) - x$fitted^(2/3)))/ (2*(x$fitted + x$fitted^2/phi))^(1/6))
                     # pearson = (results$m-results$fitted)/sqrt(results$fitted - results$fitted^2/phi),
                     # deviance = sign(results$m - results$fitted)*(2*(results$m*log(results$m/results$fitted) - (results$m + 1/phi)*log((results$m + 1/phi)/(results$fitted + 1/phi))))^(1/2)
    if (!is.null(x$countries)) {
      df$nationality <- x$countries
    }

    # 1 -- observed vs fitted (square root)
    p1 <- ggplot(df, aes(x = sqrt(fitted), y = sqrt(observed))) +
      geom_smooth(method = 'lm', se = FALSE, color = 'black',lty = 'dashed', lwd = 0.5) +
      geom_point(shape = 1) +
      labs(title = 'Observed vs Fitted', x = '√Fitted', y = '√Observed') +
      theme_classic() +
      theme(panel.border = element_rect(fill = NA))
    if (label && !is.null(x$countries)) p1 <- p1 + geom_text_repel(data = df, aes(label = nationality))

    # 2 -- normal QQ plot for Anscombe residuals
    p2 <- ggplot(df, aes(sample = anscombe)) +
      stat_qq_line() +
      stat_qq(shape = 1) +
      labs(title = 'Q-Q Plot for Anscombe residuals', x = 'Theoretical Quantiles', y = 'Sample Quantiles') +
      theme_classic() +
      theme(panel.border = element_rect(fill = NA))

    # 3 -- comparison of square fitted and absolute Anscombe residuals
    p3 <- ggplot(df, aes(x = sqrt(fitted), y = abs(anscombe))) +
      geom_point(shape = 1) +
      geom_smooth(se = FALSE, color = 'red', lwd = 0.5) +
      # geom_hline(yintercept = 0, linetype = 'dashed', color = 'red') +
      labs(title = 'Absolute Anscombe residuals vs √Fitted', x = '√Fitted', y = '|Anscombe|') +
      theme_classic() +
      theme(panel.border = element_rect(fill = NA))
    if (label && !is.null(x$countries)) p3 <- p3 + geom_text_repel(data = df, aes(label = nationality))

    # 4 -- comparison of square fitted and Anscombe residuals
    p4 <- ggplot(df, aes(x = sqrt(fitted), y = anscombe)) +
      geom_hline(yintercept = 0, linetype = 'dotted') +
      geom_point(shape = 1) +
      geom_smooth(se = FALSE, color = 'red', lwd = 0.5) +
      labs(title = 'Anscombe residuals vs √Fitted', x = '√Fitted', y = 'Anscombe') +
      theme_classic() +
      theme(panel.border = element_rect(fill = NA))
    if (label && !is.null(x$countries)) p4 <- p4 + geom_text_repel(data = df, aes(label = nationality))

    if (!is.null(which) && which %in% seq(1:4)){
      plot <- switch(which, '1' = p1, '2' = p2, '3' = p3, '4' = p4)
      print(plot)
    } else{
      grid.arrange(p1, p2, p3, p4, ncol = 2, nrow = 2)
    }

  }

}
