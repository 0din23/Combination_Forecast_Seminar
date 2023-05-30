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
  mutate_at(colnames(.)[-1], as.numeric) %>% 
  filter(!is.na(Index_excess)) %>% 
  mutate(
    Index_hist_mean = rollapplyr(Index_excess, width = seq_along(Index_excess), mean) %>% lag()
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
  mutate_at(colnames(.)[-1], as.numeric)

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
                                    m = 60,
                                    date = res[["yyyyq"]]
                             ) %>% 
      select(date, alpha, beta, y_hat, epsilon, feature_value = x, y) %>% 
      mutate("feature" = x)
    return(res)
  }) %>% 
  rbindlist() %>% 
  left_join(., data %>% select(date = yyyyq, Index_hist_mean), by = c("date")) %>% 
  filter(!is.na(alpha), between(date, as.yearqtr("1947 Q1"), as.yearqtr("2005 Q4"))) %>%
  group_by(feature) %>% 
  mutate(
    epsilon_hist = Index_hist_mean - y, 
    Net_SSE = cumsum(epsilon_hist^2 - epsilon^2)
  ) %>% 
  ungroup(feature)
  

# Simple mean combination forecast
t = 165
m = 1
s = c(m:(t-1))
combination_forecast <- univariate_forecast %>% 
  select(date,feature,y,y_hat,epsilon_hist) %>%
  group_by(date) %>%
  mutate(mean_weight = 1/n(),
         mean_combination = sum(mean_weight * y_hat),
         median_combination = median(y_hat),
         trimmed_mean_weight = 1/(n()-2),
         trimmed_mean_combination = sum(trimmed_mean_weight[c(-which.max(y_hat),-which.min(y_hat))]*y_hat[c(-which.max(y_hat),-which.min(y_hat))])) %>%
  group_by(feature) %>%
  filter(.,row_number() == s) %>%
  mutate(discount0.9 = 0.9^(t-1-s),
         phi0.9 = discount0.9 * (y-y_hat)^2,
         discount1 = 1^(t-1-s),
         phi1 = discount1 * (y-y_hat)^2) %>%
  group_by(date) %>%
  mutate(dmspe_weight0.9 = phi0.9^(-1)/sum(phi0.9^(-1)),    ####???
         dmspe_combination0.9 = sum(dmspe_weight0.9 * y_hat),
         dmspe_weight1 = phi1^(-1)/sum(phi1^(-1)),          ####???
         dmspe_combination1 = sum(dmspe_weight1 * y_hat)) %>%
  group_by(date) %>%
  mutate(error1= cumsum((unique(epsilon_hist))^2-(unique(mean_combination-y))^2),
         error2= cumsum((unique(epsilon_hist))^2-(unique(median_combination-y))^2),
         error3= cumsum((unique(epsilon_hist))^2-(unique(trimmed_mean_combination-y))^2),
         error4= cumsum((unique(epsilon_hist))^2-(unique(dmspe_combination0.9-y))^2),
         error5= cumsum((unique(epsilon_hist))^2-(unique(dmspe_combination1-y))^2)) %>% 
  select(date,error1,error2,error3,error4,error5) %>%
  distinct(date, .keep_all = TRUE) 

multivariate_forecast <- data %>% 
  select(date = yyyyq, Index_excess_forward, all_of(predictive_features)) %>%
  filter(between(date, as.yearqtr("1947 Q1"), as.yearqtr("2005 Q4"))) %>% 
  expanding_multi_regression(df = ., label = "Index_excess_forward", m = 60) %>% 
  left_join(., data %>% select(yyyyq, Index_hist_mean), by = c("date" = "yyyyq")) %>% 
  filter(!is.na(epsilon)) %>% 
  mutate(
    epsilon_hist = Index_hist_mean - Index_excess_forward, 
    Net_SSE = cumsum(epsilon_hist^2 - epsilon^2)
  )

################################################################################
# Plots #
################################################################################
# Univariate
univariate_forecast %>% 
    as.data.frame() %>% 
    ggplot(.) + 
    geom_line(aes(x = date, y = Net_SSE)) + 
    facet_wrap(~feature)
    
multivariate_forecast %>% 
  ggplot(.) +
  geom_line(aes(x = date, y = Net_SSE))
# Combination Forecast
# combination_forecast %>% 
#   mutate(
#     Mean= cumsum(error1),
#     Median = cumsum(error2),
#     Trimmed_mean = cumsum(error3),
#     DMSPE1.0 = cumsum(error5),
#     DMSPE0.9 = cumsum(error4)
#   ) %>%
#   select(date, Mean, Median, Trimmed_mean, DMSPE1.0, DMSPE0.9) %>% 
#   pivot_longer(data=., cols = colnames(.)[colnames(.)!="date"],
#                names_to = "feature", values_to = "Net_SSE") %>% 
#   ggplot(.) + 
#   geom_line(aes(x = date, y = Net_SSE)) + 
#   facet_wrap(~feature)

data.frame(date = combination_forecast$date, Mean = cumsum(combination_forecast$error1),Median = cumsum(combination_forecast$error2),
           Trimmed_mean = cumsum(combination_forecast$error3),
           DMSPE1.0 = cumsum(combination_forecast$error5),
           DMSPE0.9 = cumsum(combination_forecast$error4))%>%
  select(date, Mean, Median, Trimmed_mean, DMSPE1.0, DMSPE0.9) %>% 
  pivot_longer(data=., cols = colnames(.)[colnames(.)!="date"],
               names_to = "feature", values_to = "Net_SSE") %>% 
  ggplot(.) + 
  geom_line(aes(x = date, y = Net_SSE)) + 
  facet_wrap(~feature)
