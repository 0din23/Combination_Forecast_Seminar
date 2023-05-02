expanding_uni_regression <- function(y, x, date, m,...){
  
  # Convert inputs to data frame
  res <- data.frame(
    "y" = as.numeric(y),
    "x" = as.numeric(x),
    "date" = date
  )
  
  # Calculate expanding window regression coefficients
  res <- res %>% 
    cbind(., res%>% 
            select(-date) %>% 
            rollapplyr(., width = seq_along(x), FUN = function(d) {
              d <- as.data.frame(d)
              if(nrow(d) <= m){
                return(c(NA,NA))
              } else{
                reg <- lm(as.numeric(d[-nrow(d),1]) ~ as.numeric(d[-nrow(d),2]))
                return(coef(reg)) 
              }
              
            }, by.column = FALSE)
    )
  
  # Rename new variabels
  colnames(res) <- c("y", "x", "date", "alpha", "beta")
  
  # Calculate predicted values and residuals
  res$y_hat <- res$alpha + res$beta * res$x
  res$epsilon <- res$y - res$y_hat
  
  return(res)
}
