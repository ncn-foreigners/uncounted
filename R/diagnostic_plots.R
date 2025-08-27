
#'@export
plot.hidden <- function(results, which = NULL){
  
  if (results$method == 'ols') {
    
    df <- data.frame(fitted = results$fitted, 
                     observed = results$m,
                     residuals = results$residuals,
                     resid_stand = results$resid_stand,
                     cooks = results$cooks,
                     leverage = results$leverage)
    
    # 1 residuals vs fitted
    p1 <- ggplot(df, aes(x = fitted, y = residuals)) +
      geom_point(shape = 1) + 
      geom_hline(yintercept = 0, linetype = 'dotted') + 
      labs(title = 'Residuals vs Fitted values', x = 'Fitted values', y = 'Residuals') + 
      theme_classic() +
      theme(panel.border = element_rect(fill = NA))
    
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
    
    # 5 observed vs fitted
    p5 <- ggplot(df, aes(x = fitted, y = observed)) + 
      geom_smooth(method = 'lm', se = FALSE, color = 'black',lty = 'dashed', lwd = 0.5) +
      geom_point(shape = 1) +
      labs(title = 'Observed vs Fitted', x = 'Fitted', y = 'Observed') + 
      theme_classic() +
      theme(panel.border = element_rect(fill = NA))
    
    if (!is.null(which) && which %in% seq(1:5)){
      plot <- switch(which, '1' = p1, '2' = p2, '3' = p3, '4' = p4, '5' = p5)
      print(plot)
    } else{
      grid.arrange(p1, p2, p3, p4, ncol = 2, nrow = 2)
    }
    
  } else if (results$method == 'nls') {
    
    df <- data.frame(fitted = results$fitted, 
                     observed = results$m,
                     residuals = results$residuals)
    
    # 1 residuals vs fitted
    p1 <- ggplot(df, aes(x = fitted, y = residuals)) +
      geom_point(shape = 1) + 
      geom_hline(yintercept = 0, linetype = 'dotted') + 
      labs(title = 'Residuals vs Fitted values', x = 'Fitted values', y = 'Residuals') + 
      theme_classic() +
      theme(panel.border = element_rect(fill = NA))
    
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
    
    if (!is.null(which) && which %in% seq(1:3)){
      plot <- switch(which, '1' = p1, '2' = p2, '3' = p3)
      print(plot)
    } else{
      grid.arrange(p1, p2, ncol = 2)
    }
    
  } else if (results$method == 'glm - Poisson') {
    
    df <- data.frame(fitted = results$fitted, 
                     observed = results$m,
                     residuals = results$residuals,     # deviance residuals 
                     resid_stand = results$resid_stand, # standardized deviance residuals 
                     cooks = results$cooks,
                     leverage = results$leverage)
                     # pearson = (results$m - results$fitted)/sqrt(results$fitted), 
                     # anscombe = (3*results$m^(2/3)-3*results$fitted^(2/3))/(2 * results$fitted^(1/6))
    
    # 1 residuals vs fitted
    p1 <- ggplot(df, aes(x = fitted, y = residuals)) +
      geom_point(shape = 1) + 
      geom_hline(yintercept = 0, linetype = 'dotted') + 
      labs(title = 'Residuals vs Fitted values', x = 'Fitted values', y = 'Residuals') + 
      theme_classic() +
      theme(panel.border = element_rect(fill = NA))
    
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
    
    # 5 observed vs fitted
    p5 <- ggplot(df, aes(x = fitted, y = observed)) + 
      geom_smooth(method = 'lm', se = FALSE, color = 'black',lty = 'dashed', lwd = 0.5) +
      geom_point(shape = 1) +
      labs(title = 'Observed vs Fitted', x = 'Fitted', y = 'Observed') + 
      theme_classic() +
      theme(panel.border = element_rect(fill = NA))
    
    if (!is.null(which) && which %in% seq(1:5)){
      plot <- switch(which, '1' = p1, '2' = p2, '3' = p3, '4' = p4, '5' = p5)
      print(plot)
    } else{
      grid.arrange(p1, p2, p3, p4, ncol = 2, nrow = 2)
    }
    
  } else if (results$method == 'mle') {
    
    phi <- results$coefficients$phi
    df <- data.frame(fitted = results$fitted, 
                     observed = results$m,
                     residuals = results$residuals,
                     anscombe = (3*phi*(1 + results$m/phi)^(2/3) - (1 + results$fitted/phi)^(2/3) + 3*(results$m^(2/3) - results$fitted^(2/3)))/ (2*(results$fitted + results$fitted^2/phi))^(1/6))
                     # pearson = (results$m-results$fitted)/sqrt(results$fitted - results$fitted^2/phi),
                     # deviance = sign(results$m - results$fitted)*(2*(results$m*log(results$m/results$fitted) - (results$m + 1/phi)*log((results$m + 1/phi)/(results$fitted + 1/phi))))^(1/2)
    
  
    # 1 -- observed vs fitted (square root)
    p1 <- ggplot(df, aes(x = sqrt(fitted), y = sqrt(observed))) + 
      geom_smooth(method = 'lm', se = FALSE, color = 'black',lty = 'dashed', lwd = 0.5) +
      geom_point(shape = 1) +
      labs(title = 'Observed vs Fitted', x = '√Fitted', y = '√Observed') + 
      theme_classic() +
      theme(panel.border = element_rect(fill = NA))
    
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

    # 4 -- comparison of square fitted and Anscombe residuals
    p4 <- ggplot(df, aes(x = sqrt(fitted), y = anscombe)) +
      geom_hline(yintercept = 0, linetype = 'dotted') +
      geom_point(shape = 1) +
      geom_smooth(se = FALSE, color = 'red', lwd = 0.5) +
      labs(title = 'Anscombe residuals vs √Fitted', x = '√Fitted', y = 'Anscombe') +
      theme_classic() +
      theme(panel.border = element_rect(fill = NA))

    if (!is.null(which) && which %in% seq(1:4)){
      plot <- switch(which, '1' = p1, '2' = p2, '3' = p3, '4' = p4)
      print(plot)
    } else{
      grid.arrange(p1, p2, p3, p4, ncol = 2, nrow = 2)
    }
    
  }
  
}
