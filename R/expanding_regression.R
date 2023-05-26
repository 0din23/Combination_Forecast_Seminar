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
################################################################################
# df <- data %>% 
#   select(date = yyyyq, Index_excess_forward, all_of(predictive_features)) %>% 
#   filter(between(date, as.yearqtr("1947 Q1"), as.yearqtr("2005 Q4")))
# label <- "Index_excess_forward"
# m <- 60

expanding_multi_regression <- function(df, label, m,...){
  
  f <- paste0(label, " ~ .")
  
  # Calculate expanding window regression coefficients
  res <- df %>% 
    cbind(., df %>% 
            select(-date) %>% 
            rollapplyr(., width = seq_along(df$Index_excess_forward), FUN = function(d) {
              d <- as.data.frame(d) %>% 
                mutate_all(as.numeric)
           
              if(nrow(d) <= m){
                return(rep(NA, ncol(df)-1))
              } else{
                reg <- lm(f, data = d)
                return(coef(reg)) 
              }
              
            }, by.column = FALSE)
    )
  
  # Rename new variabels
  colnames(res) <- c(colnames(df), "alpha", paste("ß_",colnames(df)[-c(1,2)]))
  
  # Calculate predicted values and residuals
  res$y_hat <- res$alpha + rowSums(as.data.frame(as.matrix(res[,paste("ß_",colnames(df)[-c(1,2)])]) * as.matrix(res[,colnames(df)[-c(1,2)]])), na.rm = T)
  res$epsilon <- res$y - res$y_hat
  return(res)
}
