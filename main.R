# Sourcing 
source("R/dependencies.R")

################################################################################
# Data and Input #
################################################################################
data <- read.csv(PATH_DATA) %>% 
  mutate(
    yyyyq = yyyyq %>% as.character() %>% as.yearqtr(., "%Y%q"),
    Index_excess = RETURN(Index) - Rfree,
    Index_excess_forward = shift(RETURN(Index) - Rfree,-1),
    Index_hist_mean = rollapplyr(RETURN(Index)- Rfree, width = seq(1,nrow(.),1), FUN = mean) 
  ) %>%
  .[-1,] %>% 
  mutate(
    Index_hist_mean = rollapplyr(Index_excess, width = seq(1,nrow(.),1), FUN = mean)
  )

data <- list(
  "complete"=data,
  "long_sample"= data %>% filter(yyyyq >= "1965 Q1"),
  "short_sample"=data %>% filter(yyyyq >= "2006 Q1"),
  "original_sample"= data %>% filter((yyyyq >= "1965 Q1" & yyyyq <= "2005 Q4"))
)


################################################################################
# Single Parameter Predictive Regression #
################################################################################
# "D3" "E3" is somehow to short
predictive_features <- c("D12", "E12", "b.m", "tbl", "AAA", "BAA", "lty", "cay",
                         "ntis", "infl", "ltr", "corpr", "svar", "csp",  "ik")
original_sample_reg_df <- predictive_features %>%
  lapply(., function(x){
    print(paste0("Feature: ", x))
    res <- expanding_uni_regression(x = data$original_sample[[x]],
                                    y = data$original_sample[["Index_excess_forward"]],
                                    m = 60, date = data$original_sample[["yyyyq"]]) %>% 
      select(date, alpha, beta, y_hat, epsilon, feature_value = x, y) %>% 
      mutate(
        "feature" = x
      )
    return(res)
  }) %>% 
  rbindlist() %>% 
  left_join(., data$original_sample %>% 
              select(date = yyyyq, Index_hist_mean),
            by = c("date")) %>% 
  mutate(epsilon_hist = y - Index_hist_mean) %>% 
  #na.omit() %>% 
  group_by(feature) %>% 
  mutate(
    Net_SSE = cumsum(epsilon^2 - epsilon_hist^2)
  ) %>% 
  ungroup(feature)
  pivot_wider(data=., names_from = feature, values_from = c("alpha", "beta", "y_hat", "epsilon"))



################################################################################
# Plots #
################################################################################
original_sample_reg_df %>% 
    as.data.frame() %>% 
    ggplot(.) + 
    geom_line(aes(x = date, y = epsilon^2 - epsilon_hist^2)) + 
    facet_wrap(~feature)
    


