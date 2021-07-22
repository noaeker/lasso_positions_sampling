
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
lasso_path = glmnet(x = X, y = y, nlambda = 1000, lambda.min.ratio = 0.0000001, alpha =1,lower.limits=0)
#lasso = glmnet(x = X, y = y, nlambda = 100, alpha =1,lower.limits=0,relax = TRUE, gamma =0)
lasso_path_coefs_sum_tidy = tidy(lasso_path$df)
lasso_path_coefs_matrix_tidy = tidy(lasso_path$beta)
lasso_path_tidy = tidy(lasso_path)

coefs_cnt <- rowid_to_column(lasso_path_coefs_sum_tidy) %>% rename(step = rowid, non_zero_cnt = x)  
n_loci = ncol(X)
pcts = c(0.005,0.01,0.05,0.1,0.15, max(coefs_cnt$non_zero_cnt)/n_loci)
counts = pcts*n_loci
relevant_non_zero_cnts<-coefs_cnt %>% crossing(pcts) %>%  mutate(counts = as.integer(pcts*n_loci),diff = (non_zero_cnt-counts)) %>% filter (diff>=0) %>% group_by(counts, pcts) %>% filter (diff ==min(diff)) %>% pull(non_zero_cnt) 
#lasso_path = glmnet(x = X, y = y, nlambda = 100, alpha =1,lower.limits=0, relax = TRUE)
output <- lasso_path_tidy %>% inner_join(coefs_cnt) %>% filter (non_zero_cnt %in% relevant_non_zero_cnts) %>% rename(loci = term)
write_csv(output,paste(str_trim(args[2]),"/r_lasso.csv",sep=""))
if (as.logical(as.integer(args[3])))
{relevant_lambdas<- output %>% distinct (lambda) %>% pull(lambda)
lasso_path_relaxed <-  glmnet(x = X, y = y, lambda = relevant_lambdas, alpha =1,lower.limits=0, relax = TRUE)
lasso_path_relaxed_tidy = tidy(lasso_path_relaxed)
relaxed_coefs_sum_tidy = tidy(lasso_path_relaxed$df)
relaxed_coefs_cnt <- rowid_to_column(relaxed_coefs_sum_tidy) %>% rename(step = rowid, non_zero_cnt = x)  
relaxed_output <- lasso_path_relaxed_tidy  %>% inner_join(relaxed_coefs_cnt) %>% rename(loci = term)

write_csv(relaxed_output,paste(str_trim(args[2]),"/r_lasso_relaxed.csv", sep=""))}





