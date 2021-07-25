
library(glmnet)
library(tidyverse)
library(broom)



options(echo=TRUE) # if you want see commands in output file
args <- commandArgs(trailingOnly = TRUE)

print(args)
#rm(args)
training_data<-read_csv(args[1]) 
training_data<- training_data %>% rowwise() %>% mutate (y=sum(c_across(is.numeric)))
y <- training_data %>% pull(y)
X <- training_data %>% select (-y)
# fdev = 1.0e-5, devmax = 0.999
lasso_path = glmnet(x = X, y = y, nlambda = 1000, lambda.min.ratio = 0.00000000001, alpha =1,lower.limits=0,devmax = 0.999999,fdev = 1.0e-10)
#lasso = glmnet(x = X, y = y, nlambda = 100, alpha =1,lower.limits=0,relax = TRUE, gamma =0)
lasso_path_coefs_sum_tidy = tidy(lasso_path$df)
lasso_path_coefs_matrix_tidy = tidy(lasso_path$beta)
lasso_path_tidy = tidy(lasso_path)

obtained_coefs_cnt <- rowid_to_column(lasso_path_coefs_sum_tidy) %>% rename(step = rowid, obtained_non_zero_cnt = x)  
n_loci = ncol(X)
#0.005,
desired_pcts = as.numeric(unlist((str_split(args[4],"_"))))#c(0.01,0.05,0.1,0.15))

output <- lasso_path_tidy %>% inner_join(obtained_coefs_cnt) %>% group_by(obtained_non_zero_cnt) %>% filter (lambda == min(lambda) )%>% ungroup()  %>%
  rename(loci = term) %>% mutate(obtained_non_zero_pct = obtained_non_zero_cnt/n_loci) %>% crossing(desired_pcts) %>% mutate(diff = (obtained_non_zero_pct-desired_pcts)) %>% filter (diff>=0) %>% group_by(desired_pcts) %>% filter (diff == min(diff)) %>% ungroup() %>% arrange(desired_pcts, loci)
write_csv(output,paste(str_trim(args[2]),"/r_lasso.csv",sep=""))

if (as.logical(as.integer(args[3])))
{relevant_lambdas<- output %>% distinct (lambda) %>% pull(lambda)
lasso_path_relaxed <-  glmnet(x = X, y = y, lambda = relevant_lambdas, alpha =1,lower.limits=0, relax = TRUE)
lasso_path_relaxed_tidy = tidy(lasso_path_relaxed)
relaxed_coefs_sum_tidy = tidy(lasso_path_relaxed$df)
relaxed_obtained_coefs_cnt <- rowid_to_column(relaxed_coefs_sum_tidy) %>% rename(step = rowid, obtained_non_zero_cnt = x)  
relaxed_output <- lasso_path_relaxed_tidy %>% inner_join(relaxed_obtained_coefs_cnt) %>% group_by(obtained_non_zero_cnt) %>% filter (lambda == min(lambda) )%>% ungroup()  %>%
  rename(loci = term) %>% mutate(obtained_non_zero_pct = obtained_non_zero_cnt/n_loci) %>% crossing(desired_pcts) %>% mutate(diff = (obtained_non_zero_pct-desired_pcts)) %>% filter (diff>=0) %>% group_by(desired_pcts) %>% filter (diff == min(diff)) %>% ungroup() %>% arrange(desired_pcts, loci)
write_csv(relaxed_output,paste(str_trim(args[2]),"/r_lasso_relaxed.csv", sep=""))}





