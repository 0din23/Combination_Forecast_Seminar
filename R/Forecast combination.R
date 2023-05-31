###combination forecasting model

library(ggpubr)


####determine the weights

t = 165
m = 1
s = c(m:(t-1))
weight <- original_sample_reg_df %>%
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

#########%>%
# weight <- weight %>% mutate( Mean= cumsum(error1),
#         Median = cumsum(error2),
#        Trimmed_mean = cumsum(error3),
#        DMSPE1.0 = cumsum(error5),
#        DMSPE0.9 = cumsum(error4))
  
 
  
  
  
  
  
  
  
  
###############plot###########
  
  
dff = data.frame(date = weight$date, Mean = cumsum(weight$error1),Median = cumsum(weight$error2),
                 Trimmed_mean = cumsum(weight$error3),
                 DMSPE1.0 = cumsum(weight$error5),
                 DMSPE0.9 = cumsum(weight$error4)) 

p1 = dff %>% ggplot(.) + 
  geom_line(aes(x = date, y = Mean)) 
  
p2 = dff %>% ggplot(.) + 
  geom_line(aes(x = date, y = Median)) 

p3 = dff %>% ggplot(.) + 
  geom_line(aes(x = date, y = Trimmed_mean)) 

p4 = dff %>% ggplot(.) + 
  geom_line(aes(x = date, y = DMSPE1.0)) 

p5 = dff %>% ggplot(.) + 
  geom_line(aes(x = date, y = DMSPE0.9)) 
  
ggarrange(p1, p2,p3,p4,p5, ncol = 3, nrow = 2)  





##################forecast_evaluation##############
q0 = 40  ####ten years initial holdout period
q = 50   ####q=T-m-q0
m = 40   ####first m observations ,in_sample portion
forecast_evaluation <- original_sample_reg_df %>% 
  group_by(feature) %>%
  select(feature,y,y_hat,Index_hist_mean) %>%
  filter(., between(row_number(), m+q0+1, m+q0+q)) %>%
  group_by(feature) %>%
  mutate(OOSR2 = 1- sum((y-y_hat)^2)/sum((y-Index_hist_mean)^2))%>% 
  ungroup(feature)





