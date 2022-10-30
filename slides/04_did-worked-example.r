#####################################################
# Difference-in-Differences for Categorical Outcomes
#####################################################

#-- Setup 
library(tidyverse)
library(expm)
library(glue)
library(fastDummies)
library(here)
library(broom)
library(gt)
library(Matrix)
library(nnet)

category_names <- 
  c("ESI","NG","PUB","UNIN")                        
names(category_names) <-  
  c("Employer Sponsored","Privately Purchased","Public","Uninsured")

df_ex <- read_rds("./simulated-data-example-data.rds")

# This matrix will eventually summarize marginal effect estimtes across approaches
marginal_att_summary <- matrix(nrow = length(category_names), ncol = 6)
rownames(marginal_att_summary) <- category_names
colnames(marginal_att_summary) <- c("1. PT-Nonparametric", "2. PT-LinearReg", "3. PT-MutinomLogit",
                                    "4. TR-Nonparametric", "5. TR-LinearReg", "6. TR-MultinomLogit")

#########################################################
# Part A. Estimation Under a Parallel Trends Assumption
#########################################################

#---- 1. Nonparametric 2x2

tbl_2x2 <- 
  df_ex %>% 
  fastDummies::dummy_cols(select_columns = "insurance_type") %>% 
  rename_at(vars(starts_with("insurance_type_")), ~gsub("insurance_type_","",.)) %>% 
  group_by(z,post) %>% 
  summarise_at(vars(all_of(category_names)),mean) %>% 
  gather(occupancy, value,-z,-post) %>% 
  mutate(group_time = glue("z{z}_t{post}")) %>% 
  ungroup() %>% 
  select(-z,-post) %>% 
  spread(group_time,value) %>% 
  mutate(att = (z1_t1 - z1_t0) - (z0_t1 - z0_t0)) %>% 
  select(occupancy,att,everything())

marginal_att_summary[,"1. PT-Nonparametric"] <- tbl_2x2$att

#---- 2. Linear Regression

df_ <- 
  df_ex %>% 
  # Dichotomize the categorical outcome variable insurance_type
  dummy_cols(select_columns = "insurance_type") 

#### Run a series of J linear probability models. Note, however, the sum-to-zero
#### constraint on the collection of J DID estimates means we really
#### only need to run J-1 models to get the full set of DID estimates. 
#### We include all so we can show the regression results across all 
#### four outcome categories. 

fit_meff_lpm_ESI  <- lm(insurance_type_ESI ~ post + z + post_z + X, data = df_)

fit_meff_lpm_NG   <- lm(insurance_type_NG  ~ post + z + post_z + X, data = df_)

fit_meff_lpm_PUB  <- lm(insurance_type_PUB ~ post + z + post_z + X, data = df_)

fit_meff_lpm_UNIN <- lm(insurance_type_UNIN ~ post + z + post_z + X, data = df_)


# Collect the regression models in a single list for processing output
fit_meff_lpm <- list(ESI = fit_meff_lpm_ESI, 
                     NG = fit_meff_lpm_NG,
                     PUB = fit_meff_lpm_PUB,
                     UNIN = fit_meff_lpm_UNIN)

tbl_lpm <- 
  fit_meff_lpm %>% 
  map(~(
    .x %>% 
    tidy() 
  )) %>% 
  bind_rows(.id = "outcome") %>% 
  mutate_at(vars(outcome,term),function(x) gsub("_|_\\$",".",x)) %>% 
  select(outcome,term,estimate,std.error) %>% 
  mutate_at(vars(estimate,std.error),function(x) round(x,4)) %>% 
  mutate(term = gsub("\\(|\\)","",term)) %>% 
  mutate(estimate = glue("{estimate} ({std.error})")) %>% 
  select(-std.error) %>% 
  spread(outcome,estimate) 

tbl_lpm %>% 
  gt() %>% 
  tab_header(title = "Marginal Effect Model Coefficients (std. error) Based on Linear Probability Model")  %>% 
  cols_label(term = "")

marginal_att_summary[,("2. PT-LinearReg")] <- 
  fit_meff_lpm %>% 
  map(~(
    .x %>% 
      tidy() 
  )) %>% 
  bind_rows(.id = "outcome") %>% 
  mutate_at(vars(outcome,term),function(x) gsub("_|_\\$",".",x)) %>% 
  filter(term=="post.z") %>% 
  pull(estimate)


#---- 3. Multinomial Logit

# Setup the data frame
df_ <- 
  df_ex %>% 
  mutate(insurance_type = factor(insurance_type, levels = c("ESI","NG","PUB","UNIN")))

fit_meff_mltnm <- multinom(insurance_type ~ post + z + post_z + X, data = df_)

tbl_multinom <- 
  tidy(fit_meff_mltnm ) %>% 
  rename(outcome = y.level) %>% 
  mutate_at(vars(outcome,term),function(x) gsub("_|_\\$",".",x)) %>% 
  select(outcome,term,estimate,std.error) %>% 
  mutate_at(vars(estimate,std.error),function(x) round(x,4)) %>% 
  mutate(term = gsub("\\(|\\)","",term)) %>% 
  mutate(estimate = glue("{estimate} ({std.error})")) %>% 
  select(-std.error) %>% 
  spread(outcome,estimate) 

tbl_multinom %>% 
  gt() %>% 
  tab_header(title = "Marginal Effect Model Coefficients Based on Multinomial Logit Model")  %>% 
  cols_label(term = "")

# Obtain DID estimates of the ATT via the method of recycled predictions. 
recpr_meff_mltnm <-
  expand.grid(z = c(0,1) , post =c(0,1),X=0)  %>% 
  as_tibble() %>% 
  mutate(pr = map2(z,post,
     ~(as_tibble(predict(fit_meff_mltnm,
                         newdata = (df_ %>% filter(z==1) %>%  
                                    mutate(z = .x , post = .y, post_z = .x * .y)),
                 "probs"))))) %>% 
  unnest(cols = c(pr)) %>% 
  group_by(z,post) %>% 
  summarise_at(vars(category_names),mean) %>% 
  gather(category,value,-z,-post) %>% 
  mutate(z = glue("z{z}"),
         post = glue("t{post}")) %>% 
  unite(group,z,post) %>% 
  spread(group,value)  %>% 
  mutate(att = (z1_t1 - z1_t0) - (z0_t1 - z0_t0))

tbl_mltnm_recp <- 
  recpr_meff_mltnm %>% 
  select(coverage = category,ATT = att,everything())

tbl_mltnm_recp  %>% 
  gt() %>% 
  cols_label(z0_t0 = "Pre-Intervention",
             z0_t1 = "Post-Intervention",
             z1_t0 = "Pre-Intervention",
             z1_t1 = "Post-Intervention") %>% 
  tab_spanner(columns = c(z0_t0,z0_t1),label = "Comparison Group") %>% 
  tab_spanner(columns = c(z1_t0,z1_t1),label = "Intevention Group") %>% 
  tab_spanner(columns = c(ATT), label = "DID Estimate") %>% 
  cols_label(coverage = "") %>% 
  fmt_number(columns = c(z0_t0,z0_t1,z1_t0,z1_t1,ATT), n_sigfig=3) %>% 
  tab_header(title = "Marginal Effect Estimates Based on Recycled Predictions from a Multinomial Logit Model")

marginal_att_summary[,c("3. PT-MutinomLogit")] <- tbl_mltnm_recp$ATT

#####################################################################
# Part B. Estimation Under a Multiplicative (Transitions) Assumption
#####################################################################

#---- Tabulate the baseline distribution of the outcome in each group.
p_0 <- 
  df_ex %>% 
  filter(post==0) %>% 
  dummy_cols(select_columns = "insurance_type") %>% 
  group_by(z) %>% 
  summarise_at(vars(starts_with("insurance_type_")),mean) %>% 
  rename_at(vars(starts_with("insurance_type_")),function(x) gsub("insurance_type_","",x))

p_0z1 <- p_0 %>% filter(z==1) %>% select(c("ESI","NG","PUB","UNIN")) %>% as.matrix() 
p_0z1

p_0z0 <- p_0 %>% filter(z==0) %>% select(c("ESI","NG","PUB","UNIN")) %>% as.matrix() 
p_0z0

#---- 3. Nonparametric: Tabulation of Occupancy changes

# Tabulate the pre-to-post transition matrix in each group. 
R_np <- 
  df_ex %>% 
  filter(post==1 ) %>% 
  count(insurance_type,pre,z) %>% 
  group_by(pre,z) %>% 
  mutate(n=n/sum(n)) %>% 
  spread(insurance_type,n) %>% 
  arrange(z,pre)

R_0_np <- 
  R_np %>% 
  filter(z==0) %>% 
  data.frame() %>% 
  select(-z) %>% 
  column_to_rownames(var = "pre") %>% 
  as.matrix()
R_0_np <- R_0_np[c("ESI","NG","PUB","UNIN"),c("ESI","NG","PUB","UNIN")]
R_0_np

R_1_np <- 
  R_np %>% 
  filter(z==1) %>% 
  data.frame() %>% 
  select(-z) %>% 
  column_to_rownames(var = "pre") %>% 
  as.matrix()
R_1_np <- R_1_np[c("ESI","NG","PUB","UNIN"),c("ESI","NG","PUB","UNIN")]
R_1_np

R_DD_np = R_1_np - R_0_np
R_DD_np

# Baseline occupancy in the treated group
p_ <- p_0 %>% filter(z==1)  %>% select(-z) %>% gather(coverage,value)  %>% pull(value)

# Marginal treatment effect estimate is p_0 %*% R_DD
pi_ <- p_ %*% R_DD_np[category_names,category_names] %>% data.frame() %>%  gather(coverage,pi)

tbl_pi_nonpara <-   
  p_0 %>% 
  filter(z==1) %>% 
  select(-z) %>% 
  gather(coverage,p) %>% 
  left_join(
    R_DD_np %>% 
    data.frame() %>% 
    rownames_to_column(var = "coverage"),"coverage") %>% 
  left_join(pi_,"coverage") 
  
tbl_pi_nonpara %>% 
  select(coverage,pi,p,ESI,NG,PUB,UNIN) %>% 
  gt() %>% 
  tab_header(title = "Summary of Marginal Effect Point Estimates Using Nonparametric Transitions Estimator")  %>% 
  cols_label(coverage = "",pi ="ATT",p="p") %>% 
  fmt_number(columns = c(p,ESI,NG,PUB,UNIN,pi),n_sigfig=3) %>% 
  tab_spanner(columns = c(ESI,NG,PUB,UNIN),label = "R_DD")

marginal_att_summary[,c( "4. TR-Nonparametric")] <- tbl_pi_nonpara$pi

#---- 5. Parametric: Linear Regression

ref_cat<- "UNIN"

# Define J-1 outcomes and structure the data with pre-period indicators
  Y_ <- 
    df_ex %>% 
    dummy_cols(select_columns = "insurance_type") %>% 
    rename_at(vars(starts_with("insurance_type_")),~gsub("insurance_type_","",.)) %>% 
    filter(post == 1) %>% 
    select_at(all_of(unname(category_names)[-grep(ref_cat,category_names)]))  %>% 
    as.matrix()
  
  # Collection of pre-period outcome indicators
  y_pre <- 
    df_ex %>% 
    dummy_cols(select_columns = "insurance_type") %>% 
    rename_at(vars(starts_with("insurance_type_")),~gsub("insurance_type_","",.)) %>% 
    filter(post == 0) %>% 
    select_at(vars(idnumber,all_of(category_names)))
  
  # Interactions of pre-period outcome indicators with treatment variable. 
  y_pre_z <- 
    df_ex %>% 
    filter(post == 0 ) %>% 
    select(-post_z) %>% 
    dummy_cols(select_columns = "insurance_type") %>% 
    rename_at(vars(starts_with("insurance_type_")),~gsub("insurance_type_","",.))  %>% 
    mutate_at(all_of(unname(category_names)),~as.integer((. * .data$z))) %>% 
    rename_at(vars(all_of(unname(category_names))),~(paste0(.,"_z"))) %>% 
    select(idnumber,ends_with("_z"))
  
  # Include covariate X
  xx <- 
    df_ex %>% 
    filter(post ==0) %>% 
    select(idnumber,X)
  
  # Combine them all together to ensure we maintain proper matching of elements by study unit. 
  df_ <- 
    y_pre %>% 
    inner_join(y_pre_z,"idnumber")  %>% 
    inner_join(xx,"idnumber") %>% 
    select(-idnumber)
  
  fit_att_transitions <-
    lm(Y_ ~ . - 1 , data = df_) 

tidy(fit_att_transitions) %>% 
  mutate_at(vars(response,term),function(x) gsub("_|_\\$",".",x)) %>% 
  select(response,term,estimate,std.error) %>% 
  mutate_at(vars(estimate,std.error),function(x) round(x,4)) %>% 
  mutate(term = gsub("\\(|\\)","",term)) %>% 
  mutate(term = gsub("\\(|\\)","",term)) %>% 
  mutate(response = gsub("insurance.type.","",response)) %>% 
  mutate(estimate = glue("{round(estimate,3)} ({round(std.error,3)})")) %>% 
  select(-std.error) %>% 
  spread(response,estimate) %>% 
  gt() %>% 
  tab_header("Linear Transition Regression Model Coefficients (SE)") %>% 
  cols_label(term = "") %>% 
  tab_spanner(columns = c(ESI,NG,PUB),label = "Outcome")

# Prepare a data frame for the recycled predictions.
xx_ <- 
    category_names %>% 
    map_df(~({
      data.frame(z = c(1,0)) %>% 
        mutate(pre = .x) %>% 
        mutate(X = 0)
    })) %>% 
    fastDummies::dummy_cols(select_columns = c("pre")) %>% 
    rename_at(vars(starts_with("pre_")),~gsub("pre_","",.)) %>% 
    mutate_at(vars(all_of(category_names)),list(z = ~(.data$z * .))) 

predicted_values_m_transitions <- 
  cbind(xx_[,c("z","pre")],predict(fit_att_transitions,newdata =xx_)) %>% 
  left_join(p_0 %>% gather(pre,p,-z),c("pre","z")) %>% 
  data.frame() %>% 
  mutate(pre = factor(pre,levels = c(category_names))) %>% 
  arrange(z,pre) %>% 
  select(z,pre,p,everything())
predicted_reference <- 
  1-rowSums(predicted_values_m_transitions[,category_names[-which(category_names==ref_cat)]])
predicted_values_m_transitions[[ref_cat]] <-  predicted_reference

R_0_linreg <- 
  predicted_values_m_transitions %>% 
  ungroup() %>% 
  filter(z == 0) %>% 
  select(-z,-p) %>% 
  column_to_rownames(var = "pre") %>% 
  as.matrix()
# Make sure all the categories line up in the same order ...
R_0_linreg <- R_0_linreg[unname(category_names),unname(category_names)]

R_1_linreg <- 
  predicted_values_m_transitions %>% 
  ungroup() %>% 
  filter(z == 1) %>% 
  select(-z,-p) %>% 
  column_to_rownames(var = "pre") %>% 
  as.matrix()
# Make sure all the categories line up in the same order ...
R_1_linreg <- R_1_linreg[unname(category_names),unname(category_names)]

R_DD_linreg <- R_1_linreg - R_0_linreg

pi_parametric_transitions <- p_ %*% R_DD_linreg[category_names,category_names] %>% data.frame() %>%  gather(coverage,pi)

tbl_pi_regression <-   
  p_0 %>% 
  filter(z==1) %>% 
  select(-z) %>% 
  gather(coverage,p) %>% 
  left_join(
    R_DD_linreg %>% 
    data.frame() %>% 
    rownames_to_column(var = "coverage"),"coverage") %>% 
  left_join(pi_parametric_transitions,"coverage") 
  
tbl_pi_regression %>% 
  select(coverage,pi,p,ESI,NG,PUB,UNIN) %>% 
  gt() %>% 
  tab_header(title = "Summary of Marginal Effect Point Estimates Using Regression-Based Transitions Estimator")  %>% 
  cols_label(coverage = "",pi ="ATT",p="p") %>% 
  fmt_number(columns = c(p,ESI,NG,PUB,UNIN,pi),n_sigfig=3) %>% 
  tab_spanner(columns = c(ESI,NG,PUB,UNIN),label = "R_DD")

marginal_att_summary[,c("5. TR-LinearReg")] <- tbl_pi_regression$pi

#---- 6. Parametric: Multinomial Logit

Y_ <- 
  df_ex %>% 
  filter(post ==1) %>% 
  select(Y = insurance_type) %>% 
  mutate(Y = factor(Y, levels = category_names)) %>% 
  as.matrix()

fit_att_transitions_multinom <-
  multinom(Y_ ~ . - 1 , data = df_)

tbl_multinom2 <- 
  tidy(fit_att_transitions_multinom ) %>% 
  rename(outcome = y.level) %>% 
  mutate_at(vars(outcome,term),function(x) gsub("_|_\\$",".",x)) %>% 
  select(outcome,term,estimate,std.error) %>% 
  mutate_at(vars(estimate,std.error),function(x) round(x,4)) %>% 
  mutate(term = gsub("\\(|\\)","",term)) %>% 
  mutate(estimate = glue("{estimate} ({std.error})")) %>% 
  select(-std.error) %>% 
  spread(outcome,estimate) 

tbl_multinom2 %>% 
  gt() %>% 
  tab_header(title = "Tansition Model Coefficients (std. error) Based on Multinomial Logit Model")  %>% 
  cols_label(term = "")

predicted_values_m_transitions_multinom <- 
  cbind(xx_[,c("z","pre")],predict(fit_att_transitions_multinom,newdata =xx_,"probs")) %>% 
  left_join(p_0 %>% gather(pre,p,-z),c("pre","z")) %>% 
  data.frame() %>% 
  mutate(pre = factor(pre,levels = c(category_names))) %>% 
  arrange(z,pre) %>% 
  select(z,pre,p,everything())
predicted_reference_multinom <- 
  1-rowSums(predicted_values_m_transitions_multinom[,category_names[-which(category_names==ref_cat)]])
predicted_values_m_transitions_multinom[[ref_cat]] <-  predicted_reference_multinom

hat_R_0_multinom <- 
  predicted_values_m_transitions_multinom %>% 
  ungroup() %>% 
  filter(z == 0) %>% 
  select(-z,-p) %>% 
  column_to_rownames(var = "pre") %>% 
  as.matrix()
# Make sure all the categories line up in the same order ...
hat_R_0_multinom <- hat_R_0_multinom[category_names,category_names]

R_1_multinom <- 
  predicted_values_m_transitions_multinom %>% 
  ungroup() %>% 
  filter(z == 1) %>% 
  select(-z,-p) %>% 
  column_to_rownames(var = "pre") %>% 
  as.matrix()
# Make sure all the categories line up in the same order ...
R_1_multinom <- R_1_multinom[category_names,category_names]

R_DD_multinom <- R_1_multinom - hat_R_0_multinom

pi_multinom_transitions <- p_ %*% R_DD_multinom[category_names,category_names] %>% data.frame() %>%  gather(coverage,pi)

tbl_pi_multinom <-   
  p_0 %>% 
  filter(z==1) %>% 
  select(-z) %>% 
  gather(coverage,p) %>% 
  left_join(
    R_DD_multinom %>% 
    data.frame() %>% 
    rownames_to_column(var = "coverage"),"coverage") %>% 
  left_join(pi_multinom_transitions,"coverage") 
  
tbl_pi_multinom %>% 
  select(coverage,pi,p,ESI,NG,PUB,UNIN) %>% 
  gt() %>% 
  tab_header(title = "Summary of Marginal Effect Point Estimates Using Multinomial Regression-Based Transitions Estimator")  %>% 
  cols_label(coverage = "",pi ="ATT",p="p") %>% 
  fmt_number(columns = c(p,ESI,NG,PUB,UNIN,pi),n_sigfig=3) %>% 
  tab_spanner(columns = c(ESI,NG,PUB,UNIN),label = "R_DD")

marginal_att_summary[,c("6. TR-MultinomLogit")] <- tbl_pi_multinom$pi

#####################################################################
# Part C. Estimation of R_DD Under a Parallel Trends Assumption 
#####################################################################

df_ <- 
  df_ex %>% 
  # Dichotomize the categorical outcome variable insurance_type
  dummy_cols(select_columns = "insurance_type") %>% 
  rename_at(vars(starts_with("insurance_type_")), ~gsub("insurance_type_","",.))

R_DD_pt <- matrix(0,nrow = length(category_names), ncol = length(category_names), 
                 dimnames = list(category_names, category_names))

for (rr in c("ESI","NG","PUB","UNIN")) {
  df_tmp <- df_ %>% filter(pre==rr)
  
  R_DD_pt[rr,"ESI"]  <- coef(lm(ESI ~ post + z + post_z , data = df_tmp))["post_z"] 
  R_DD_pt[rr,"NG"]    <- coef(lm(NG  ~ post + z + post_z , data = df_tmp))["post_z"] 
  R_DD_pt[rr,"PUB"]   <- coef(lm(PUB ~ post + z + post_z , data = df_tmp))["post_z"] 

}
# Estimates for the reference category are obtained by leveraging the sum-to-zero constraint
# on the full set of DID estimates: 
R_DD_pt[,"UNIN"] <- -Matrix::rowSums(R_DD_pt)

# To crosswalk back to the marginal effects estimate under a parallel trends assumption we 
# need to add in the "intervention effect" second term (see paper)

p_0z1 %*% R_DD_pt + (p_0z1 %*% R_0_np - p_0z1) - (p_0z0 %*% R_0_np - p_0z0)

# Once we do we arrive at the same marginal effect as in Part A. of this document
tbl_2x2$att


