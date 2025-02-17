uniform_quantile_scale <- function(x, lower_quantile = 0.33, upper_quantile = 0.99, center_quantile = NULL) {
  
  lower <- quantile(x, lower_quantile)
  upper <- quantile(x, upper_quantile)
  zeroq <- min(x)
  oneq <- max(x)
  if (lower<=(zeroq + 1e-6)) { 
    lower <- zeroq + 1e-6
    lower_quantile <- 1e-6
  }
  if (upper>= oneq - 1e-6) { 
    upper <- oneq - 1e-6
    upper_quantile <- 1 - 1e-6
  }
  if (!is.null(center_quantile)){
    
    center <- quantile(x, center_quantile)
    
    if (center<=lower + 1e-6) {
      center <- lower + 1e-6
      center_quantile <- lower_quantile + 1e-6
    } else if (center>=upper - 1e-6) {
      center <- upper - 1e-6
      center_quantile <- upper_quantile - 1e-6
    }
    scaled <- lower_quantile + ifelse(x < center,
                                      (x - lower) / (center - lower) *(center_quantile - lower_quantile),
                                      (center_quantile - lower_quantile)+ (x - center) / (upper - center) * (upper_quantile - center_quantile))
  } else {
    scaled <- lower_quantile + (x - lower) / (upper - lower) *(upper_quantile - lower_quantile)
  }
  # Clip values to [0, 1] range
  pmax(pmin(scaled, 1), 0)
}
