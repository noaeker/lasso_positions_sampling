library(glmnet)
library(tidyverse)
library(broom)



options(echo=TRUE) # if you want see commands in output file
args <- commandArgs(trailingOnly = TRUE)
print(args)
rm(args)
training_data<-read_csv(args[1]) 
training_data<- training_data %>% rowwise() %>% mutate (y=sum(c_across(is.numeric)))
y <- training_data %>% pull(y)
X <- training_data %>% select (-y)
lasso_path = glmnet(x = X, y = y, nlambda = 100, alpha =1,lower.limits=0, relax = args[2])
#lasso = glmnet(x = X, y = y, nlambda = 100, alpha =1,lower.limits=0,relax = TRUE, gamma =0)
lasso_path_coefs_sum_tidy = tidy(lasso_path$df)
lasso_path_coefs_matrix_tidy = tidy(lasso_path$beta)
lasso_path_tidy = tidy(lasso_path)

coefs_cnt <- rowid_to_column(lasso_path_coefs_sum_tidy) %>% rename(step = rowid, non_zero_cnt = x)  
n_loci = ncol(X)
pcts = c(0.01,0.03,0.05,0.1,0.15)
counts = pcts*n_loci
relevant_non_zero_cnts<-coefs_cnt %>% crossing(counts) %>% mutate(diff = abs(counts-non_zero_cnt)) %>% group_by(counts) %>% filter (diff ==min(diff)) %>% pull(non_zero_cnt) 
#lasso_path = glmnet(x = X, y = y, nlambda = 100, alpha =1,lower.limits=0, relax = TRUE)
output <- lasso_path_tidy %>% inner_join(coefs_cnt) %>% filter (non_zero_cnt %in% relevant_non_zero_cnts) %>% rename(loci = term)
write_csv(output,"test.csv")
if (args[3])
{relevant_lambdas<- output %>% distinct (lambda) %>% pull(lambda)
print(relevant_lambdas)
relaxed_fit <-  glmnet(x = X, y = y, lambda = relevant_lambdas, alpha =1,lower.limits=0, relax = TRUE)
lasso_path_relaxed_tidy = tidy(relaxed_fit)
relaxed_output<-lasso_path_relaxed_tidy  %>% rename(loci = term)
write_csv(output,"test_relaxed.csv")}