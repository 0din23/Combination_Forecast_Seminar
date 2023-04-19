
data <- read.csv(PATH_DATA) %>% 
  mutate(
    yyyymm = yyyymm %>% as.character() %>% as.yearmon(., "%Y%m"),
    Index_forward = shift(Index,-1)
  )


expanding_uni_regression <- function(y, x, date, m){
  
  res <- data.frame(
    "y" = y %>% as.numeric(),
    "x" = x %>% as.numeric(),
    "date" = date
  )
  
  res <- res %>% 
    cbind(., res%>% 
            select(-date) %>% 
            rollapplyr(., width = seq(1,nrow(.),1), FUN = function(d) {
              d <- as.data.frame(d)
              if(nrow(d) <= m){
                return(c(NA,NA))
              } else{
                reg <- lm(as.numeric(d[-nrow(d),1]) ~ as.numeric(d[-nrow(d),2]))
                return(coef(reg)) 
              }
      
            }, by.column = FALSE)
    )
  
  colnames(res) <- c("y", "x", "date", "alpha", "beta")
  res$y_hat <- res$alpha + res$beta * res$x
  res$epsilon <- res$y - res$y_hat
  res$y_naive <- res$y %>% lag()
  res$epsilon_naive <- res$y - res$y_naive
  res <- res %>% na.omit()
  res$Net_SSE <- cumsum(res$epsilon_naive^2 -  res$epsilon^2)
  
  return(res)
}

df <- cbind(as.yearmon(data %>% ), res) %>% as.data.frame()
colnames(df)[1] <- "date"
res %>% 
  ggplot(.) +
  geom_line(aes(x=date, y=Net_SSE))



uni_df <- expanding_uni_regression(y = data$Index_forward,
                         x = data$E12,
                         date = data$yyyymm,
                         m = 60)

