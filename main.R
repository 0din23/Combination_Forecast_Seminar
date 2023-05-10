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
  )

# long_sample: yyyyq >= "1965 Q1"
# short_sample: yyyyq >= "2006 Q1"
# original_sample: (yyyyq >= "1965 Q1" & yyyyq <= "2005 Q4")

################################################################################
# Single Parameter Predictive Regression #
################################################################################
predictive_features <- c("D/P", "D/Y", "E/P", "D/E", "SVAR", "B/M", "NTIS",
                         "TBL", "LTY", "LTR", "TMS", "DFY", "DFR", "INFL", "I/K")

original_sample_reg_df <- predictive_features %>%
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
  }) %>% 
  rbindlist() %>% 
  left_join(., data %>% select(date = yyyyq, Index_hist_mean), by = c("date")) %>% 
  filter(!is.na(alpha), between(date, as.yearqtr("1965 Q1"), as.yearqtr("2005 Q4"))) %>%
  group_by(feature) %>% 
  mutate(
    epsilon_hist = Index_hist_mean - y, 
    Net_SSE = cumsum(epsilon_hist^2 - epsilon^2)
  ) %>% 
  ungroup(feature)
  
################################################################################
# Plots #
################################################################################
original_sample_reg_df %>% 
    as.data.frame() %>% 
    ggplot(.) + 
    geom_line(aes(x = date, y = Net_SSE)) + 
    facet_wrap(~feature)
    


