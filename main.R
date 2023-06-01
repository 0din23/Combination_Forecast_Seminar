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
    Index_excess = (as.numeric(CRSP_SPvwx) - Rfree),
    Index_excess_forward = shift(as.numeric(CRSP_SPvwx) - Rfree, -1),
  ) %>% 
  filter(!is.na(Index_excess)) %>% 
  mutate(
    Index_hist_mean = rollapplyr(Index_excess, width = seq_along(Index_excess), mean)
  ) %>% 
  arrange(yyyyq) %>% 
  select(yyyyq, Index, Index_excess, Index_excess_forward, "B/M" = "b/m",
         SVAR = svar, "I/K" = ik, "NTIS" = ntis, LTY = lty,
         LTR = ltr, E12, D12, "TBL" = Rfree, AAA, BAA, Index_hist_mean,
         "INFL" = infl, CORPR = corpr) %>% 
  mutate(
    "D12" = D12 %>% as.numeric(),
    "E12" = E12 %>% as.numeric(),
    "LTY" = LTY %>% as.numeric(),
    "LTR" = LTR %>% as.numeric(),
    "AAA" = AAA %>% as.numeric(),
    "BAA" = BAA %>% as.numeric(),
    "INFL" = INFL %>% as.numeric(),
    "CORPR" = CORPR %>% as.numeric(),
    "D/Y" = log(D12 / lag(Index)),
    "D/P" = log(D12 / Index),
    "E/P" = E12 / Index,
    "D/E" = log(D12 / E12),
    "TMS" = LTY - TBL,
    "DFY" = AAA - BAA,
    "DFR" = CORPR - LTR
  ) %>% 
  mutate(
    estimated_var_excess =runSD(Index_excess, 40),
    estimated_var = runSD(Index_excess + TBL, 40),
  )
# 
data <- filter(data,yyyyq >= "1947 Q1" & yyyyq <= "2005 Q4")
# data <- filter(data,yyyyq >= "1955 Q1")
# data <- filter(data,yyyyq >= "1995 Q4")
# long_sample: yyyyq >= "1965 Q1"
# short_sample: yyyyq >= "2006 Q1"
# original_sample: (yyyyq >= "1965 Q1" & yyyyq <= "2005 Q4")

################################################################################
# Single Parameter Predictive Regression #
################################################################################
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
  filter(date >= "1955 Q1" & date <= "2005 Q4") %>% 
  # filter(date >= "1965 Q1") %>% 
  # filter(date >= "2005 Q4") %>% 
  group_by(feature) %>% 
  mutate(
    epsilon_hist = Index_hist_mean - y, 
    Net_SSE = cumsum( epsilon_hist^2 - epsilon^2)
  ) %>% 
  ungroup(feature)


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
  filter(date >= "1965 Q1" & date <= "2005 Q4") %>% 
  group_by(feature) %>% 
  mutate(
    epsilon = y-y_hat,
    Net_SSE = cumsum(epsilon_hist^2 - epsilon^2)
  ) %>% 
  ungroup(feature)


m = 40
theta09 <- 0.9
theta1 <- 1

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
  filter(date >= "1965 Q1" & date <= "2005 Q4") %>% 
  group_by(feature) %>% 
  mutate(
    epsilon = y-y_hat,
    Net_SSE = cumsum(epsilon_hist^2 - epsilon^2)
  ) %>% 
  select(all_of(combination_forecast %>% colnames()))
  
combination_forecast <- combination_forecast %>% 
  rbind(.,combination_forecast_2)

################################################################################
# Plots #
################################################################################
# Univariate
univariate_forecast %>% 
  filter(date >= "1965 Q1" & date <= "2005 Q4") %>% 
  ggplot(.) + 
  geom_line(aes(x = date, y = Net_SSE)) + 
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
  filter(date >= "1965 Q1" & date <= "2005 Q4") %>% 
  select(all_of(necessary_columns)) %>% 
  rbind(., combination_forecast %>%
          filter(date >= "1965 Q1" & date <= "2005 Q4") %>% 
          select(all_of(necessary_columns))) %>% 
  left_join(., data %>% select(yyyyq,estimated_var_excess, estimated_var, TBL), by = c("date"="yyyyq"))

# R2os
eval_data %>% 
  group_by(feature) %>% 
  summarize(
    Ros = 1 - (sum(epsilon^2) / sum(epsilon_hist^2))
  )


# utility gain
gamma <- 3

eval_data %>% 
  mutate(
    omega_0 = (1/gamma) * (Index_hist_mean / estimated_var_excess^2),
    omega_j = (1/gamma) * (y_hat / estimated_var_excess^2),
    
    port_bench = lag(0.6*y) + 0.4*((1+TBL)^(1/4)-1),
    portfolio_hist_mean = lag((omega_0*y)) + (1-omega_0)*((1+TBL)^(1/4)-1),
    portfolio_forecast = lag((omega_j*y)) + (1-omega_j)*((1+TBL)^(1/4)-1),
  ) %>% 
  group_by(feature) %>% 
  summarize(
    utility_60_40 = mean(port_bench,na.rm=T) - (1/2)*sd(port_bench,na.rm=T)^2,
    utility_0 = mean(portfolio_hist_mean,na.rm=T) - (1/2)*sd(portfolio_hist_mean,na.rm=T)^2,
    utility_j = mean(portfolio_forecast,na.rm=T) - (1/2)*sd(portfolio_forecast,na.rm=T)^2,
    
    utility_gain_bench = utility_j- utility_60_40
    utility_gain_hist_mean = utility_j - utility_0
  )



eval_data %>% 
  mutate(
    omega_0 = (1/gamma) * (Index_hist_mean / estimated_var_excess^2),
    omega_j = (1/gamma) * (y_hat / estimated_var_excess^2),
    
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
  #filter(!(feature %in% c("LTY", "TBL"))) %>% 
  ggplot(.) +
  geom_line(aes(x=date, y = portfolio_hist_mean_cum, color = "hist. mean")) +
  geom_line(aes(x=date, y = portfolio_forecast_cum, color = "Forecast")) +
  geom_line(aes(x=date, y = port_bench_cum, color = "Benchmark")) +
  facet_wrap(~feature)
D



