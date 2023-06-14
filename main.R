# Sourcing 
source("R/dependencies.R")

################################################################################
# Data and Input #
################################################################################
# Read data from xlsx file and mutate the data frame
data <- openxlsx::read.xlsx(xlsxFile = PATH_DATA) %>% 
  mutate(
    yyyyq = yyyyq %>% as.character() %>% as.yearqtr(., "%Y%q"),
    Rfree = Rfree %>% as.numeric(),
    Index = Index %>% gsub(",", "", .) %>% as.numeric(),
    #Index_excess = (as.numeric(CRSP_SPvwx) - Rfree),
    Index_excess = log(1 + as.numeric(CRSP_SPvw)) - log(1 + as.numeric(Rfree)),
    Index_excess_forward = shift(Index_excess, -1),
  ) %>% 
  filter(!is.na(Index_excess)) %>% 
  filter(yyyyq >= "1947 Q1") %>% 
  mutate(
    Index_hist_mean = (cumsum(Index_excess)-as.numeric(Index_excess)) / (c(1:nrow(.))-1)
  ) %>% 
  arrange(yyyyq) %>% 
  select(yyyyq, Index, Index_excess, Index_excess_forward, "B/M" = "b/m",
         SVAR = svar, "NTIS" = ntis, LTY = lty,
         LTR = ltr, E12, D12, "TBL" = tbl, AAA, BAA, Index_hist_mean,
         "INFL" = infl, ik, CORPR = corpr, CRSP_SPvw, Rfree) %>% 
  mutate(
    "I/K" = ik %>% as.numeric() %>% lag(),
    "D12" = D12 %>% as.numeric(),
    "E12" = E12 %>% as.numeric(),
    "LTY" = LTY %>% as.numeric(),
    "LTR" = LTR %>% as.numeric(),
    "AAA" = AAA %>% as.numeric(),
    "BAA" = BAA %>% as.numeric(),
    "INFL" = INFL %>% as.numeric() %>% lag(),
    "CORPR" = CORPR %>% as.numeric(),
    "D/Y" = log(D12 / lag(Index)),
    "D/P" = log(D12 / Index),
    "E/P" = E12 / Index,
    "D/E" = log(D12 / E12),
    "TMS" = LTY - as.numeric(TBL),
    "DFY" = AAA - BAA,
    "DFR" = CORPR - LTR
  ) %>% 
  mutate(
    estimated_var_excess =runSD(Index_excess, 40),
    estimated_var = runSD(CRSP_SPvw, 40),
  ) %>%
  filter(yyyyq >= "1947 Q1")


os_backhold <- "1995 Q1"
os_start <- "2005 Q4"
end <- "2022 Q4"
# data <- filter(data,yyyyq >= "1955 Q1")
# data <- filter(data,yyyyq >= "1995 Q4")
# long_sample: yyyyq >= "1965 Q1"
# short_sample: yyyyq >= "2006 Q1"
# original_sample: (yyyyq >= "1965 Q1" & yyyyq <= "2005 Q4")

################################################################################
# Single Parameter Predictive Regression #
################################################################################

# Create Univariate Forecasts
predictive_features <- c("D/P", "D/Y", "E/P", "D/E", "SVAR", "B/M", "NTIS",
                         "TBL", "LTY", "LTR", "TMS", "DFY", "DFR", "INFL", "I/K")

univariate_forecast <- predictive_features %>%
  lapply(., function(x){
    print(paste0("Feature: ", x))
    res <- data[!is.na(as.numeric(data[[x]])),]
    res <- res[!is.na(as.numeric(data[["Index_excess_forward"]])),]
    
    res <- expanding_uni_regression(x = res[[x]],
                                    y = res[["Index_excess_forward"]],
                                    m = 40,
                                    date = res[["yyyyq"]]
    ) %>% 
      select(date, alpha, beta, y_hat, epsilon, feature_value = x, y) %>% 
      mutate("feature" = x)
  }) %>% 
  rbindlist() %>% 
  filter(!is.na(alpha)) %>%
  left_join(., data %>% select(date = yyyyq, Index_hist_mean), by = c("date")) %>% 
  filter(date >= os_backhold & date <= end) %>% 
  group_by(feature) %>% 
  mutate(
    epsilon_hist = Index_hist_mean - y
  ) %>% 
  ungroup(feature)

## campbell thompson forecast
pf_df <- data.frame(
  "feature" = predictive_features,
  "CT_sing"=c(1, 1, 1, 1, 1, 1, -1,
           -1, -1, 1, -1, -1, 1, -1, -1)
)

univariate_forecast <- univariate_forecast %>% 
  left_join(., pf_df, by = "feature") %>% 
  mutate(
    CT_beta = ifelse(beta*CT_sing>0, beta,0),
    CT_y_hat = ifelse(CT_beta > 0, alpha + CT_beta * feature_value, Index_hist_mean),
    CT_y_hat = ifelse(CT_y_hat < 0, 0, CT_y_hat)
  )



## Combination Forecast
combination_forecast <- univariate_forecast%>% 
  select(date, feature, y_hat) %>% 
  pivot_wider(data=., names_from = feature, values_from = y_hat) %>% 
  mutate(
    simple_mean = rowMeans(select(.,-date)),
    simple_median = apply(select(.,-date),1,median),
    trimmed_mean = apply(select(.,-date),1,function(x){
      x[max(x)] <- 0
      x[min(x)] <- 0
      return(sum(x)/(length(x)-2))
    })
  ) %>% 
  select(date, simple_mean, simple_median, trimmed_mean) %>% 
  pivot_longer(data=., cols = c(simple_mean, simple_median, trimmed_mean),
               names_to = "feature", values_to = "y_hat") %>%
  left_join(.,
            univariate_forecast %>% 
              filter(feature == "D/P") %>% 
              select(date, y, Index_hist_mean, epsilon_hist),
            by = "date") %>% 
  filter(date >= os_start & date <= end) %>% 
  group_by(feature) %>% 
  mutate(
    epsilon = y-y_hat,
    Net_SSE = cumsum(epsilon_hist^2 - epsilon^2)
  ) %>% 
  ungroup(feature)


m = 40
theta09 <- 0.9
combination_forecast_2 <- univariate_forecast %>% 
  group_by(feature) %>% 
  mutate(
    t = c(1:n()),
    theta_09 = theta09^(t-1-m),
    
    DMSPE_09_weight = cumsum(theta_09 * epsilon^2),
    DMSPE_1_weight = cumsum(epsilon^2),
    
  ) %>% 
  group_by(date) %>% 
  mutate(
    DMPSPE09 = DMSPE_09_weight / sum(DMSPE_09_weight),
    DMPSPE1 = DMSPE_1_weight / sum(DMSPE_1_weight),
  
    DMPSPE09 = ifelse(t > m, DMPSPE09, 0),
    DMPSPE1 = ifelse(t > m, DMPSPE1, 0) ,
  ) %>% 
  summarize(
    DMPSPE09 = sum(DMPSPE09 * y_hat),
    DMPSPE1 = sum(DMPSPE1 * y_hat),
    y = mean(y),
    Index_hist_mean = mean(Index_hist_mean),
    epsilon_hist = mean(epsilon_hist)
  ) %>% 
  pivot_longer(data=., cols = c("DMPSPE09", "DMPSPE1"), values_to = "y_hat",
               names_to = "feature") %>% 
  filter(date >= os_start & date <= end) %>% 
  group_by(feature) %>% 
  mutate(
    epsilon = y-y_hat,
    Net_SSE = cumsum(epsilon_hist^2 - epsilon^2)
  ) %>% 
  select(all_of(combination_forecast %>% colnames()))


## Combine all data
combination_forecast <- combination_forecast %>% 
  rbind(.,combination_forecast_2) 

univariate_forecast <- univariate_forecast %>% 
  filter(date >= os_start & date <= end) %>% 
  group_by(feature) %>% 
  mutate(
    epsilon_ct = y-CT_y_hat,
    Net_SSE = cumsum(epsilon_hist^2 - epsilon^2),
    Net_SSE_ct = cumsum(epsilon_hist^2 - epsilon_ct^2),
  ) %>% 
  ungroup(feature)



## White Noise Test 
set.seed(1234)
wn_test <- paste0("WN_",c(1:250)) %>%
  lapply(., function(x){
    print(paste0("Feature: ", x))
    
    res <- data[!is.na(as.numeric(data[["Index_excess_forward"]])),]
    
    res <- expanding_uni_regression(x = rnorm(n=nrow(res)),
                                    y = res[["Index_excess_forward"]],
                                    m = 40,
                                    date = res[["yyyyq"]]
    ) %>% 
      select(date, alpha, beta, y_hat, epsilon, feature_value = x, y) %>% 
      mutate("feature" = x)
  }) %>% 
  rbindlist() %>% 
  filter(!is.na(alpha)) %>%
  left_join(., data %>% select(date = yyyyq, Index_hist_mean), by = c("date")) %>% 
  filter(date >= os_backhold & date <= end) %>% 
  group_by(feature) %>% 
  mutate(
    epsilon_hist = Index_hist_mean - y
  ) %>% 
  ungroup(feature) %>% 
  filter(date >= os_start & date <= end) %>% 
  group_by(feature) %>% 
  mutate(
    Net_SSE = cumsum(epsilon_hist^2 - epsilon^2)
  ) %>% 
  ungroup(feature)

wn_realization <- wn_test %>% 
  group_by(feature) %>% 
  summarize(
    net_SSE = sum(epsilon^2)
  ) %>% 
  filter(
    net_SSE <= quantile(net_SSE, probs = 0.05)
  ) %>%
  arrange(net_SSE) %>% 
  tail(.,1) %>% 
  pull(feature)

################################################################################
# Plots #
################################################################################
# Univariate
univariate_forecast %>% 
  select(-c(CT_sing, CT_beta, CT_y_hat, epsilon_ct)) %>% 
  rbind(., wn_test %>% filter(feature == wn_realization) %>% 
          mutate(Net_SSE_ct = Net_SSE)) %>% 
  ggplot(.) + 
  geom_line(aes(x = date, y = Net_SSE, color = "Original")) + 
  geom_line(aes(x = date, y = Net_SSE_ct, color = "CT Restricted")) + 
  facet_wrap(~feature) + 
  #scale_y_continuous(limits = c(-0.1, .1), breaks = seq(-.1, .1, by = .05)) + 
  geom_hline(yintercept = 0) +  # Add horizontal line at y = 0
  theme_bw() +
  ggtitle(paste0("Evaluation Plots for: ",univariate_forecast$date[1], " to ", univariate_forecast$date %>% tail(.,1)))+
  theme(axis.text.x = element_text(angle = 70, vjust=0, hjust = 0))


combination_forecast %>% 
  ggplot(.) + 
  geom_line(aes(x = date, y = Net_SSE)) + 
  facet_wrap(~feature) + 
  #scale_y_continuous(limits = c(-0.1, .1), breaks = seq(-.1, .1, by = .05)) + 
  geom_hline(yintercept = 0) +  # Add horizontal line at y = 0
  theme_bw() +
  ggtitle(paste0("Evaluation Plots for: ",univariate_forecast$date[1], " to ", univariate_forecast$date %>% tail(.,1)))

################################################################################
# Tables #
################################################################################

necessary_columns <- c("date", "feature", "epsilon", "epsilon_hist", "Index_hist_mean", "y","y_hat")

eval_data <- univariate_forecast %>% 
  filter(date >= os_start & date <= end) %>% 
  select(all_of(necessary_columns)) %>% 
  rbind(., combination_forecast %>%
          filter(date >= os_start & date <= end) %>% 
          select(all_of(necessary_columns))) %>% 
  left_join(., data %>% select(yyyyq,estimated_var_excess, estimated_var, TBL, CRSP_SPvw), by = c("date"="yyyyq"))

# R2os
eval_data %>% 
  na.omit() %>% 
  group_by(feature) %>% 
  summarize(
    Ros = 100*(1 - (sum((y-y_hat)^2) / sum((y-Index_hist_mean)^2)))
  )

# compute p-value
#T = length(unique(eval_data$date))
#q0=20 #initial holdout period
eval_data = eval_data %>% mutate(f = epsilon_hist^2 - (epsilon^2-(Index_hist_mean-y_hat)^2))#%>%
            #group_by(feature) %>%
            #filter(.,between(row_number(), m+q0+1, T)) 
features = unique(eval_data$feature)
statistical_significance <- features %>%
  sapply(., function(x){
    df = eval_data %>% filter(feature == x) 
    f = df$f
    model = lm(f~1)
    freedom = df.residual(model)
    se <- summary(model)$coefficients[, "Std. Error"]
    t_value <- coef(model) / se
    p_value <- pt(t_value, freedom, lower.tail = FALSE)
  })
statistical_significance 


# utility gain
gamma <- 3

portfolio_eval <- eval_data %>% 
  group_by(feature) %>% 
  mutate(
    
    y_hat_simple = exp(y_hat)-1,
    
    omega_0 = (1/gamma) * (Index_hist_mean / estimated_var^2),
    omega_j = (1/gamma) * (y_hat_simple / estimated_var^2),
    
    # restriction
    omega_0 = ifelse(omega_0 > 1.5,1.5,omega_0),
    omega_0 = ifelse(omega_0 < 0,0,omega_0),
    
    omega_j = ifelse(omega_j > 1.5,1.5,omega_j),
    omega_j = ifelse(omega_j < 0,0,omega_j),
    
    
    port_bench = 0.6*as.numeric(CRSP_SPvw) + 0.4*((1+TBL)^(1/4)-1),
    portfolio_hist_mean = lag(omega_0)*as.numeric(CRSP_SPvw) + (1-lag(omega_0))*((1+TBL)^(1/4)-1),
    portfolio_forecast = lag(omega_j)*as.numeric(CRSP_SPvw) + (1-lag(omega_j))*((1+TBL)^(1/4)-1),
    
  ) %>% 
  group_by(feature) %>% 
  summarize(
    
    utility_60_40 = mean(port_bench,na.rm=T) - (gamma/2)*sd(port_bench,na.rm=T)^2,
    utility_0 = mean(portfolio_hist_mean,na.rm=T) - (gamma/2)*sd(portfolio_hist_mean,na.rm=T)^2,
    utility_j = mean(portfolio_forecast,na.rm=T) - (gamma/2)*sd(portfolio_forecast,na.rm=T)^2,
    
    utility_gain_bench = (utility_j- utility_60_40)*400,
    utility_gain_hist_mean = (utility_j - utility_0)*400,
    
    #under_perf_alloc_stock = mean((omega_j)* ifelse(portfolio_forecast < portfolio_hist_mean,1,0), na.rm=T),
    #under_perf_alloc_stock_diff = mean((omega_j-omega_0) * ifelse(portfolio_forecast < portfolio_hist_mean,1,0), na.rm=T),
    
    #outlier_cach_upside_hist_mean = mean(omega_0 * ifelse(lag(y) < quantile(y, probs = 0.8),1,0), na.rm=T),
    #outlier_cach_upside_forecast = mean(omega_j * ifelse(lag(y) < quantile(y, probs = 0.8),1,0), na.rm=T)
    
  ) %>%
  #select(feature, utility_gain_hist_mean) %>% 
  mutate(utility_gain_hist_mean = round(utility_gain_hist_mean,2))
  
portfolio_eval %>% 
  select(feature, utility_j, utility_gain_bench, utility_gain_hist_mean) %>% 
  pivot_longer(data=., cols = c("utility_j", "utility_gain_bench", "utility_gain_hist_mean")) %>% 
  ggplot(.) +
  geom_col(aes(y=value, x = feature, fill = name), position="dodge")







# Here are Niklas Plots, Please do not touch ####################################
gamma <- 3

portfolio_eval <- eval_data %>% 
  group_by(feature) %>% 
  mutate(
    
    y_hat_simple = exp(y_hat)-1,
    TBL = as.numeric(TBL),
    
    omega_0 = (1/gamma) * (Index_hist_mean / estimated_var^2),
    omega_j = (1/gamma) * (y_hat_simple / estimated_var^2),
    
    # restriction
    omega_0 = ifelse(omega_0 > 1.5,1.5,omega_0),
    omega_0 = ifelse(omega_0 < 0,0,omega_0),
    
    omega_j = ifelse(omega_j > 1.5,1.5,omega_j),
    omega_j = ifelse(omega_j < 0,0,omega_j),
    
    
    port_bench = 0.6*as.numeric(CRSP_SPvw) + 0.4*((1+TBL)^(1/4)-1),
    portfolio_hist_mean = lag(omega_0)*as.numeric(CRSP_SPvw) + (1-lag(omega_0))*((1+TBL)^(1/4)-1),
    portfolio_forecast = lag(omega_j)*as.numeric(CRSP_SPvw) + (1-lag(omega_j))*((1+TBL)^(1/4)-1),
    
  ) 
portfolio_eval%>% 
  group_by(feature) %>% 
  summarize(
    
    utility_60_40 = mean(port_bench,na.rm=T) - (gamma/2)*sd(port_bench,na.rm=T)^2,
    utility_0 = mean(portfolio_hist_mean,na.rm=T) - (gamma/2)*sd(portfolio_hist_mean,na.rm=T)^2,
    utility_j = mean(portfolio_forecast,na.rm=T) - (gamma/2)*sd(portfolio_forecast,na.rm=T)^2,
    
    utility_gain_bench = round((utility_j- utility_60_40)*400,2),
    utility_gain_hist_mean = round((utility_j - utility_0)*400,2)
  )

## Tabelle und Barchart for empirical results
barchart_features <- c("simple_mean", "trimmed_mean", "simple_median", "D/Y", "TMS",
                       "TMS", "INFL", "I/K", "SVAR", "TBL")
  
portfolio_eval %>% 
  group_by(feature) %>% 
  summarize(
    utility_60_40 = mean(port_bench,na.rm=T) - (gamma/2)*sd(port_bench,na.rm=T)^2,
    utility_0 = mean(portfolio_hist_mean,na.rm=T) - (gamma/2)*sd(portfolio_hist_mean,na.rm=T)^2,
    utility_j = mean(portfolio_forecast,na.rm=T) - (gamma/2)*sd(portfolio_forecast,na.rm=T)^2,
    
    utility_gain_bench = round((utility_j- utility_60_40)*400,2),
    utility_gain_hist_mean = round((utility_j - utility_0)*400,2),
    utility_j = round(utility_j*400,2)
  ) %>% 
  select(feature, "Utility gain vs. Benchmark" = utility_gain_bench,
         "Utility gain vs. hist. mean" = utility_gain_hist_mean) %>% 
  filter(feature %in% barchart_features) %>% 
  pivot_longer(data=., cols = c("Utility gain vs. Benchmark", "Utility gain vs. hist. mean")) %>% 
  ggplot(.) +
  geom_col(aes(y=value, x = feature, fill = name), position="dodge") +
  theme_tq() + ylab("Model") + xlab("Utility Gain") +
  theme(axis.text.x = element_text(angle = 0, vjust=0.2, hjust = 0))



## Equitylines to show performance
portfolio_eval %>% 
  filter(feature %in% barchart_features) %>%
  na.omit() %>% 
  mutate(forecast_cumret = cumprod(1+portfolio_forecast)-1,
         bench_cumret = cumprod(1+port_bench)-1) %>% 
  ggplot(.) +
  geom_line(aes(x=date, y=forecast_cumret, color = feature)) +
  geom_line(aes(x=date, y=bench_cumret, color = "Bench"))

portfolio_eval %>% 
  group_by(feature) %>% 
  na.omit() %>% 
  summarize(
    return_ann = 100*(tail(cumprod(portfolio_forecast+1),1)^(4/n())-1),
    Vol = sd(portfolio_forecast, na.rm=T)*200,
    Sharpe_simple = return_ann / Vol
  )

## Weight line chart
portfolio_eval %>% 
  filter(feature == "TMS") %>% 
  ggplot(.) +
  geom_line(aes(x=date, y=omega_j, color = feature))


## Contribution to performance
portfolio_eval %>% 
  filter(feature == "TMS") %>% 
  ggplot(.) +
  geom_point(aes(x=lag(omega_j), y=as.numeric(CRSP_SPvw)))


portfolio_eval %>% 
  na.omit() %>% 
  group_by(feature) %>% 
  summarize(
    ret_cor = mean((lag(omega_j)*as.numeric(CRSP_SPvw))^2, na.rm=T) -
      mean(lag(omega_j)*as.numeric(CRSP_SPvw), na.rm=T)^2,
    
    tbl_cor = mean((lag(omega_j)*as.numeric(TBL))^2, na.rm=T) -
      mean(lag(omega_j)*as.numeric(TBL), na.rm=T)^2,
    
    ret_cor = round(ret_cor, 3),
    tbl_cor = round(tbl_cor, 3)
  )


eval_data %>% 
  mutate(
    omega_0 = (1/gamma) * (Index_hist_mean / estimated_var_excess^2),
    omega_j = (1/gamma) * (y_hat / estimated_var_excess^2),
    
    # restriction
    omega_j = ifelse(omega_j > 1.5,1.5,omega_j),
    omega_j = ifelse(omega_j < 0,0,omega_j),
    
    port_bench = lag(0.6*y) + 0.4*((1+TBL)^(1/4)-1),
    portfolio_hist_mean = lag(omega_0*y) + (1-omega_0)*((1+TBL)^(1/4)-1),
    portfolio_forecast = lag(omega_j*y) + (1-omega_j)*((1+TBL)^(1/4)-1),
    
  ) %>% 
  na.omit() %>% 
  group_by(feature) %>% 
  mutate(
    port_bench_cum = cumprod(1+port_bench),
    portfolio_hist_mean_cum = cumprod(1+ portfolio_hist_mean),
    portfolio_forecast_cum = cumprod(1+ portfolio_forecast)
  ) %>% 
  #filter(date >= "2005 Q4") %>% 
  ggplot(.) +
  geom_line(aes(x=date, y = portfolio_hist_mean_cum-1, color = "hist. mean")) +
  geom_line(aes(x=date, y = portfolio_forecast_cum-1, color = "Forecast")) +
  geom_line(aes(x=date, y = port_bench_cum-1, color = "Benchmark")) +
  facet_wrap(~feature) +
  theme_bw() + ylab("cum. Return") +
  theme(axis.text.x = element_text(angle = 70, vjust=0, hjust = 0))


# eval_data %>% 
#   ggplot(.)+
#   #geom_histogram(aes(x = epsilon ),bins = 200) +
#   geom_density(aes(x=epsilon))+
#   geom_vline(xintercept = mean(eval_data$epsilon))+
#   facet_wrap(~feature)


eval_data %>%
  filter(feature == "TBL") %>% 
  ggplot(.)+
  #geom_histogram(aes(x = epsilon ),bins = 200) +
  geom_density(aes(x=epsilon_hist))


eval_data %>% 
  filter(feature != "SVAR") %>% 
  ggplot(.) +
  geom_line(aes(x=date, y = y_hat)) +
  facet_wrap(~feature)
  

################################################################################
# EXPLORATION #
################################################################################
data %>% 
  filter(yyyyq >= "1965 Q1") %>% 
  select(all_of(predictive_features)) %>% 
  mutate_all(as.numeric) %>% 
  cor() %>% 
  round(., 3) %>% xtable::xtable()
