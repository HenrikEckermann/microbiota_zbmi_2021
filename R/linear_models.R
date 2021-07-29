#####################################################################
##########               README                            ##########
#####################################################################

# The analyses was run on the HPC compute cluster at the Donders Institute 
# in Nijmegen. This script can be run on weaker system as well (tested 
# with 16GB and 2.7 Intel Quad). The memory of the system will be checked
# in this script and based on that either 5 or 100 imputation are used. 
# They lead to the same conclusions, but we report results based on 100
# imputations. Those will take ca. 2 days to fit on a HPC cluster.

options(repr.plot.width=30, repr.plot.height=10)

#####################################################################
##########               helper functions                  ##########
#####################################################################

# Bayesian helper functions for after model fitting 
get_lower <- function(x, level = 0.025) {quantile(x, level)}
get_upper <- function(x, level = 0.975) {quantile(x, level)}

print_estimate <- function(model, effect){
  tab <- posterior_samples(model, effect) %>%
    select(all_of(effect)) %>%
    summarise_all(
      list(beta = median, lower = get_lower, upper = get_upper)
    ) %>% mutate(across(everything(), function(x) format(round(x, 3), nsmall = 3)))
  return(glue("$\beta$ = {tab$beta}, 95% CI [{tab$lower}, {tab$upper}]"))
}

# to create tables 
brms_tbl <- function(model, exclude_id = TRUE) {
  tb <- posterior_samples(model) %>%
    pivot_longer(everything(), names_to = "Parameter", values_to = "Value") %>%
    group_by(Parameter) %>%
    summarise(
      Estimate = median(Value),
      SD = sd(Value),
      lower = quantile(Value, 0.025),
      upper = quantile(Value, 0.975),
      .groups = "drop"
    ) %>%
    mutate(across(where(is.numeric), function(x) format(round(x, 2), nsmall = 2))) %>%
    mutate("95% CI" = glue("{lower} : {upper}")) %>%
    select(Parameter, Estimate, SD, "95% CI")
    
  if (exclude_id) {
    tb <- tb  %>%
      filter(!grepl("^r_id", Parameter), Parameter != "lp__")
  }
  return(tb)
}

library(tidyverse)
library(brms)
library(mice)
library(glue)
library(patchwork)

# you must first run import_data.R to create the objects we load here:
load(here::here("rdata/mb_import_H.Rds"))

#####################################################################
##########                      setup MI                   ##########
#####################################################################

# to be able to inspect models locally I set a lower m that will give same 
# point estimate but potentially lower accuracy of credible interval. For 
# the final model the high m will be used which is fitted in the cloud.
# this happens automatically based on availabel memory (see readme)

mem <- memuse::Sys.meminfo()$totalram %>% as.numeric()
m <- ifelse(mem <= 17179869184, 5, 100)



ids <- unique(df_clr$id)
# add info for imputation 
wla <- readxl::read_xlsx(here::here("data/weight_length_all.xlsx")) %>%
  select(id = ID, everything(), -sex) %>%
  mutate(across(where(is.numeric), function(x) ifelse(x == 0, NA, x))) %>%
  arrange(id)

# add weight, length and age at measurement for imputation 
wla_test <- bind_cols(
  bind_rows(
      select(wla, id, length = Length_2, weight = Weight_2, age_wla = Age_2 ) %>% mutate(time = 1), 
      select(wla, id, length = Length_3, weight = Weight_3, age_wla = Age_3 ) %>% mutate(time = 2),
      select(wla, id, length = Length_4, weight = Weight_4, age_wla = Age_4 ) %>% mutate(time = 3), 
      select(wla, id, length = Length_6, weight = Weight_6, age_wla = Age_6 )%>% mutate(time = 4), 
      select(wla, id, length = Length_8, weight = Weight_8, age_wla = Age_8) %>% mutate(time = 5) 
  ),
  bind_rows(
      select(wla, length_tminus = Length_1, weight_tminus = Weight_1, age_wla_tminus = Age_1), 
      select(wla, length_tminus = Length_2, weight_tminus = Weight_2, age_wla_tminus = Age_2),
      select(wla, length_tminus = Length_3, weight_tminus = Weight_3, age_wla_tminus = Age_3), 
      select(wla, length_tminus = Length_5, weight_tminus = Weight_5, age_wla_tminus = Age_5), 
      select(wla, length_tminus = Length_7, weight_tminus = Weight_7, age_wla_tminus = Age_7)) 
  ) %>% filter(id %in% ids) %>%
        arrange(id, time)
df <- full_join(df, wla_test, by = c("id", "time")) %>%
      full_join(df_clr, by = c("id", "time")) 

# fix taxa name for modeling 
colnames(df) <- str_replace(colnames(df), "[pg]__", "")
colnames(df) <- str_replace(colnames(df), "\\[", "")
colnames(df) <- str_replace(colnames(df), "\\]", "")


df <- select(df, everything(), scfa_clr = scfa_tminus)
# take out redundant variable before MI
df_imp <- mice::mice(
  select(df, -sl_fb_ratio, -bw_grams, , -l_fb_ratio),
  m = m,
  maxiter = 40,
  method = "pmm"
)


# The goal is not to find perfect priors for any specific model (as I fit 
# many different models per time point and per covariate structure).
# Rather the goal is to find priors that allow for any possible assocation
# while only being slightly restricting parameter space. Given the amount of
# data we have, the chosen priors will play barely a role as long as I do not 
# prevent the model from searching in realistic parameter space. FOr such as 
# scenario, prior predictive simulations have been performed frequently with 
# data, e.g. here: https://github.com/HenrikEckermann/microbiota_executive_functioning2021/blob/main/R/sim.R # below we only specify
# prior for the slope parameters, the SD of the "random effects" and the error 
# of the model. That means that brms default priors are used for the intercept 
# (student_t) and correlation between slopes and intercepts (lkj).
# All these priors are easily "overwhelmed" by data.
priors <- c(
  prior(normal(0, 0.5), class = "b"),
  prior(exponential(1), class = "sd"),
  prior(exponential(1), class = "sigma")
)


#####################################################################
##########                      imputed models             ##########
#####################################################################




# fbratio vs BMI 
mfb1_mi <- brm_multiple(
  family = gaussian(),
  formula = bmi ~ bmi_tminus + s_fb_ratio + s_bw_grams + (1 + s_fb_ratio | id),
  prior = priors,
  data = df_imp,
  iter = 1e4,
  control = list(adapt_delta = 0.99),
  chains = 1,
  file = here::here(glue("rdata/mfb1_{m}"))
)
pp_check(mfb1_mi)
summary(mfb1_mi)



mfb1_if_mi <- mice::filter(df_imp, df$time %in% c(1, 2, 3)) %>% 
  brm_multiple(
    family = gaussian(),
    formula = bmi ~ bmi_tminus + s_fb_ratio + s_bw_grams + (1 + s_fb_ratio | id),
    prior = priors,
    data = .,
    iter = 1e4,
    control = list(adapt_delta = 0.99),
    file = here::here(glue("rdata/mfb1_if_{m}"))
  )
pp_check(mfb1_if_mi)
summary(mfb1_if_mi)


mfb1_ch_mi <- mice::filter(df_imp, df$time %in% c(4, 5)) %>% 
  brm_multiple(
    family = gaussian(),
    formula = bmi ~ bmi_tminus + s_fb_ratio + s_bw_grams + (1 + s_fb_ratio | id),
    prior = priors,
    data = .,
    iter = 1e4,
    control = list(adapt_delta = 0.99),
    file = here::here(glue("rdata/mfb1_ch_{m}"))
  )
pp_check(mfb1_ch_mi)
summary(mfb1_ch_mi)


mscfa1_mi <- brm_multiple(
  family = gaussian(),
  formula = bmi ~ bmi_tminus + scfa_clr + s_bw_grams + (1 + scfa_clr | id),
  prior = priors,
  data = df_imp,
  iter = 1e4,
  control = list(adapt_delta = 0.99),
  file = here::here(glue("rdata/mscfa1_{m}"))
)
pp_check(mscfa1_mi)
summary(mscfa1_mi)



mscfa1_if_mi <- mice::filter(df_imp, df$time %in% c(1, 2, 3)) %>% 
  brm_multiple(
    family = gaussian(),
    formula = bmi ~ bmi_tminus + scfa_clr + s_bw_grams + (1 + scfa_clr | id),
    prior = priors,
    data = .,
    iter = 1e4,
    control = list(adapt_delta = 0.99),
    file = here::here(glue("rdata/mscfa1_if_{m}"))
  )

pp_check(mscfa1_if_mi)
summary(mscfa1_if_mi)

mscfa1_ch_mi <- mice::filter(df_imp, df$time %in% c(4, 5)) %>% 
  brm_multiple(
    family = gaussian(),
    formula = bmi ~ bmi_tminus + scfa_clr + s_bw_grams + (1 + scfa_clr | id),
    prior = priors,
    data = .,
    iter = 1e4,
    control = list(adapt_delta = 0.99),
    file = here::here(glue("rdata/mscfa1_ch_{m}"))
  )

pp_check(mscfa1_ch_mi)
summary(mscfa1_ch_mi)



firm_mi <- brm_multiple(
  family = gaussian(),
  formula = bmi ~ bmi_tminus + Firmicutes + s_bw_grams + (1 + Firmicutes | id),
  prior = priors,
  data = df_imp,
  iter = 1e4,
  control = list(adapt_delta = 0.99),
  file = here::here(glue("rdata/firm_{m}"))
)
pp_check(firm_mi)
summary(firm_mi)

firm_if_mi <- brm_multiple(
  family = gaussian(),
  formula = bmi ~ bmi_tminus + Firmicutes + s_bw_grams + (1 + Firmicutes | id),
  prior = priors,
  data = mice::filter(df_imp, df$time %in% c(1, 2, 3)),
  iter = 1e4,
  control = list(adapt_delta = 0.99),
  file = here::here(glue("rdata/firm_if_{m}"))
)
pp_check(firm_if_mi)
summary(firm_if_mi)


priors <- c(
  prior(normal(0, 0.5), class = "b"),
  prior(exponential(15), class = "sd"),
  prior(exponential(1), class = "sigma")
)
firm_ch_mi <- brm_multiple(
  family = gaussian(),
  formula = bmi ~ bmi_tminus + Firmicutes + s_bw_grams + (1 + Firmicutes | id),
  prior = priors,
  data = mice::filter(df_imp, df$time %in% c(5, 6)),
  iter = 1e4,
  control = list(adapt_delta = 0.99),
  file = here::here(glue("rdata/firm_ch_{m}"))
)
pp_check(firm_ch)
summary(firm_ch_mi)

priors <- c(
  prior(normal(0, 0.5), class = "b"),
  prior(exponential(1), class = "sd"),
  prior(exponential(1), class = "sigma")
)
bact_mi <- brm_multiple(
  family = gaussian(),
  formula = bmi ~ bmi_tminus + Bacteroidetes + s_bw_grams + (1 + Bacteroidetes | id),
  prior = priors,
  data = df_imp,
  iter = 1e4,
  control = list(adapt_delta = 0.99),
  file = here::here(glue("rdata/bact_{m}"))
)
pp_check(bact_mi)
summary(bact_mi)


bact_if_mi <- brm_multiple(
  family = gaussian(),
  formula = bmi ~ bmi_tminus + Bacteroidetes + s_bw_grams + (1 + Bacteroidetes | id),
  prior = priors,
  data = mice::filter(df_imp, df$time %in% c(1, 2, 3)),
  iter = 1e4,
  control = list(adapt_delta = 0.99),
  file = here::here(glue("rdata/bact_if_{m}"))
)
pp_check(bact_if_mi)
summary(bact_if_mi)


priors <- c(
  prior(normal(0, 0.5), class = "b"),
  prior(exponential(15), class = "sd"),
  prior(exponential(1), class = "sigma")
)
bact_ch_mi <- brm_multiple(
  family = gaussian(),
  formula = bmi ~ bmi_tminus + Bacteroidetes + s_bw_grams + (1 + Bacteroidetes | id),
  prior = priors,
  data = mice::filter(df_imp, df$time %in% c(5, 6)),
  iter = 1e4,
  control = list(adapt_delta = 0.99),
  file = here::here(glue("rdata/bact_ch_{m}"))
)

pp_check(bact_ch_mi)
summary(bact_ch_mi)


priors <- c(
  prior(normal(0, 0.5), class = "b"),
  prior(exponential(15), class = "sd"),
  prior(exponential(1), class = "sigma")
)
fbsep_mi <- brm_multiple(
  family = gaussian(),
  formula = bmi ~ bmi_tminus + Bacteroidetes + Firmicutes + s_bw_grams + (1 + Bacteroidetes + Firmicutes | id),
  prior = priors,
  data = df_imp,
  iter = 1e4,
  control = list(adapt_delta = 0.99),
  file = here::here(glue("rdata/fbsep_{m}"))
)
pp_check(fbsep_mi)
summary(fbsep_mi)



fbsep_if_mi <- brm_multiple(
  family = gaussian(),
  formula = bmi ~ bmi_tminus + Bacteroidetes + Firmicutes + s_bw_grams + (1 + Bacteroidetes + Firmicutes | id),
  prior = priors,
  data = mice::filter(df_imp, df$time %in% c(1, 2, 3)),
  iter = 1e4,
  control = list(adapt_delta = 0.99),
  file = here::here(glue("rdata/fbsep_if_{m}"))
)
pp_check(fbsep_if_mi)
summary(fbsep_if_mi)



priors <- c(
  prior(normal(0, 0.5), class = "b"),
  prior(exponential(15), class = "sd"),
  prior(exponential(1), class = "sigma")
)
fbsep_ch_mi <- brm_multiple(
  family = gaussian(),
  formula = bmi ~ bmi_tminus + Bacteroidetes + Firmicutes + s_bw_grams + (1 + Bacteroidetes + Firmicutes | id),
  prior = priors,
  data = mice::filter(df_imp, df$time %in% c(5, 6)),
  iter = 1e4,
  control = list(adapt_delta = 0.99),
  file = here::here(glue("rdata/fbsep_ch_{m}"))
)
pp_check(fbsep_ch_mi)
summary(fbsep_ch_mi)




# I put a little narrower prior on the beta coefficients 
priors2 <- c(
  prior(normal(0, 0.25), class = "b"),
  prior(exponential(1), class = "sd"),
  prior(exponential(1), class = "sigma")
)



clr_imp <- brm_multiple(
  family = gaussian(),
  formula = bmi ~ bmi_tminus + Akkermansia + Alistipes + Anaerostipes + Bacteroides + Bifidobacterium + Blautia + Coprococcus_1 + Coprococcus_2 + Coprococcus_3 + Dialister + Faecalibacterium + Holdemanella + Phascolarctobacterium + Prevotella_2 + Prevotella_7 + Prevotella_9 + Roseburia + Subdoligranulum + s_bw_grams + Eubacterium_hallii_group + (1 | id),
  prior = priors2,
  data = df_imp,
  control = list(adapt_delta = 0.99),
  file = here::here(glue("rdata/clr_{m}")),
  iter = 1e4
)
pp_check(clr_imp)
 

# infancy 
clr_if_imp <- brm_multiple(
  family = gaussian(),
  formula = bmi ~ bmi_tminus + Akkermansia + Alistipes + Anaerostipes + Bacteroides + Bifidobacterium + Blautia + Coprococcus_1 + Coprococcus_2 + Coprococcus_3 + Dialister + Faecalibacterium + Holdemanella + Phascolarctobacterium + Prevotella_2 + Prevotella_7 + Prevotella_9 + Roseburia + Subdoligranulum + s_bw_grams + Eubacterium_hallii_group + (1 | id),
  prior = priors2,
  data = mice::filter(df_imp, df$time %in% c(1, 2, 3)),
  control = list(adapt_delta = 0.99),
  file = here::here(glue("rdata/clr_if_{m}")),
  iter = 1e4
)
pp_check(clr_if_imp)
summary(clr_if_imp)




clr_ch_imp <- brm_multiple(
  family = gaussian(),
  formula = bmi ~ bmi_tminus + Akkermansia + Alistipes + Anaerostipes + Bacteroides + Bifidobacterium + Blautia + Coprococcus_1 + Coprococcus_2 + Coprococcus_3 + Dialister + Faecalibacterium + Holdemanella + Phascolarctobacterium + Prevotella_2 + Prevotella_7 + Prevotella_9 + Roseburia + Subdoligranulum + s_bw_grams + Eubacterium_hallii_group + (1 | id),
  prior = priors2,
  data = mice::filter(df_imp, df$time %in% c(4, 5)),
  control = list(adapt_delta = 0.99),
  file = here::here(glue("rdata/clr_ch_{m}")),
  iter = 1e4
)
pp_check(clr_ch_imp)

#####################################################################
##########                 non-imputed                     ##########
#####################################################################


fb <- brm(
  family = gaussian(),
  formula = bmi ~ bmi_tminus + s_fb_ratio + s_bw_grams + (1 + s_fb_ratio | id),
  prior = priors,
  data = select(df, bmi, bmi_tminus, s_fb_ratio, s_bw_grams, id),
  file = here::here("rdata/fb")
)
pp_check(fb)
summary(fb)


fb_if <- brm(
  family = gaussian(),
  formula = bmi ~ bmi_tminus + s_fb_ratio + s_bw_grams + (1 + s_fb_ratio | id),
  prior = priors,
  data = filter(df, time %in% c(1, 2, 3)),
  file = here::here("rdata/fb_if")
)
pp_check(fb_if)
summary(fb_if)



fb_ch <- brm(
  family = gaussian(),
  formula = bmi ~ bmi_tminus + s_fb_ratio + s_bw_grams + (1 + s_fb_ratio | id),
  prior = priors,
  data = filter(df, time %in% c(4, 5)),
  file = here::here("rdata/fb_ch"),
  #control = list(adapt_delta = 0.9999)
)
pp_check(fb_ch)
summary(fb_ch)


scfa_clr <- brm(
  family = gaussian(),
  formula = bmi ~ bmi_tminus + scfa_clr + s_bw_grams + (1 + scfa_clr | id),
  prior = priors,
  data = df,
  control = list(adapt_delta = 0.99),
  file = here::here("rdata/scfa_clr")
)
pp_check(scfa_clr)
summary(scfa_clr)


scfa_clr_if <- brm(
  family = gaussian(),
  formula = bmi ~ bmi_tminus + scfa_clr + s_bw_grams + (1 + scfa_clr | id),
  prior = priors,
  data = df %>% filter(time %in% c(1, 2, 3)),
  #control = list(adapt_delta = 0.99),
  file = here::here("rdata/scfa_clr_if")
)
pp_check(scfa_clr_if)
summary(scfa_clr_if)


scfa_clr_ch <- brm(
  family = gaussian(),
  formula = bmi ~ bmi_tminus + scfa_clr + s_bw_grams + (1 + scfa_clr | id),
  prior = priors,
  data = df %>% filter(time %in% c(4, 5)),
  control = list(adapt_delta = 0.99),
  file = here::here("rdata/scfa_clr_ch")
)
pp_check(scfa_clr_ch)
summary(scfa_clr_ch)


firm <- brm(
  family = gaussian(),
  formula = bmi ~ bmi_tminus + Firmicutes + s_bw_grams + (1 + Firmicutes | id),
  prior = priors,
  data = df,
  iter = 2e3,
  control = list(adapt_delta = 0.99),
  file = here::here(glue("rdata/firm"))
)
pp_check(firm)
summary(firm)

firm_if <- brm(
  family = gaussian(),
  formula = bmi ~ bmi_tminus + Firmicutes + s_bw_grams + (1 + Firmicutes | id),
  prior = priors,
  data = filter(df, df$time %in% c(1, 2, 3)),
  iter = 2e3,
  control = list(adapt_delta = 0.99),
  file = here::here(glue("rdata/firm_if"))
)
pp_check(firm_if)
summary(firm_if)


priors <- c(
  prior(normal(0, 0.5), class = "b"),
  prior(exponential(15), class = "sd"),
  prior(exponential(1), class = "sigma")
)
firm_ch <- brm(
  family = gaussian(),
  formula = bmi ~ bmi_tminus + Firmicutes + s_bw_grams + (1 + Firmicutes | id),
  prior = priors,
  data = filter(df, df$time %in% c(5, 6)),
  iter = 2e3,
  control = list(adapt_delta = 0.99),
  file = here::here(glue("rdata/firm_ch"))
)
pp_check(firm_ch)
summary(firm_ch)

priors <- c(
  prior(normal(0, 0.5), class = "b"),
  prior(exponential(1), class = "sd"),
  prior(exponential(1), class = "sigma")
)
bact <- brm(
  family = gaussian(),
  formula = bmi ~ bmi_tminus + Bacteroidetes + s_bw_grams + (1 + Bacteroidetes | id),
  prior = priors,
  data = df,
  iter = 2e3,
  control = list(adapt_delta = 0.99),
  file = here::here(glue("rdata/bact"))
)
pp_check(bact)
summary(bact)


bact_if <- brm(
  family = gaussian(),
  formula = bmi ~ bmi_tminus + Bacteroidetes + s_bw_grams + (1 + Bacteroidetes | id),
  prior = priors,
  data = filter(df, df$time %in% c(1, 2, 3)),
  iter = 2e3,
  control = list(adapt_delta = 0.99),
  file = here::here(glue("rdata/bact_if"))
)
pp_check(bact_if)
summary(bact_if)


priors <- c(
  prior(normal(0, 0.5), class = "b"),
  prior(exponential(15), class = "sd"),
  prior(exponential(1), class = "sigma")
)
bact_ch <- brm(
  family = gaussian(),
  formula = bmi ~ bmi_tminus + Bacteroidetes + s_bw_grams + (1 + Bacteroidetes | id),
  prior = priors,
  data = filter(df, df$time %in% c(5, 6)),
  iter = 2e3,
  control = list(adapt_delta = 0.99),
  file = here::here(glue("rdata/bact_ch"))
)

pp_check(bact_ch)
summary(bact_ch_mi)


priors <- c(
  prior(normal(0, 0.5), class = "b"),
  prior(exponential(15), class = "sd"),
  prior(exponential(1), class = "sigma")
)
fbsep <- brm(
  family = gaussian(),
  formula = bmi ~ bmi_tminus + Bacteroidetes + Firmicutes + s_bw_grams + (1 + Bacteroidetes + Firmicutes | id),
  prior = priors,
  data = df,
  iter = 2e3,
  control = list(adapt_delta = 0.99),
  file = here::here(glue("rdata/fbsep"))
)
pp_check(fbsep)
summary(fbsep)



fbsep_if <- brm(
  family = gaussian(),
  formula = bmi ~ bmi_tminus + Bacteroidetes + Firmicutes + s_bw_grams + (1 + Bacteroidetes + Firmicutes | id),
  prior = priors,
  data = filter(df, df$time %in% c(1, 2, 3)),
  iter = 2e3,
  control = list(adapt_delta = 0.99),
  file = here::here(glue("rdata/fbsep_if"))
)
pp_check(fbsep_if)
summary(fbsep_if)



priors <- c(
  prior(normal(0, 0.5), class = "b"),
  prior(exponential(15), class = "sd"),
  prior(exponential(1), class = "sigma")
)
fbsep_ch <- brm(
  family = gaussian(),
  formula = bmi ~ bmi_tminus + Bacteroidetes + Firmicutes + s_bw_grams + (1 + Bacteroidetes + Firmicutes | id),
  prior = priors,
  data = filter(df, df$time %in% c(5, 6)),
  iter = 2e3,
  control = list(adapt_delta = 0.99),
  file = here::here(glue("rdata/fbsep_ch"))
)
pp_check(fbsep_ch)
summary(fbsep_ch)






priors2 <- c(
  prior(normal(0, 0.25), class = "b"),
  prior(exponential(1), class = "sd"),
  prior(exponential(1), class = "sigma")
)

# Prevotella and eubacterium wwere  excluded. CIs were very wide indicating 
# multicollinearity. it is possible that it covaries strongly with any of the
# other included genera and therefore the model could not estimate correctly.
# confirm:
# select(df, all_of(colnames(df)[24:length(colnames(df))])) %>%
#   cor(method = "pearson", use = "pairwise.complete.obs") %>%
#   as_tibble(rownames = NA) %>%
#   mutate(across(everything(), round, 2))

clr <- brm(
  family = gaussian(),
  formula = bmi ~ bmi_tminus + Akkermansia + Alistipes + Anaerostipes + Bacteroides + Bifidobacterium + Blautia + Coprococcus_1 + Coprococcus_2 + Coprococcus_3 + Dialister + Faecalibacterium + Holdemanella + Phascolarctobacterium + Prevotella_2 + Prevotella_7 + Prevotella_9 + Roseburia + Subdoligranulum + s_bw_grams + Eubacterium_hallii_group + (1 | id),
  prior = priors2,
  data = df,
  control = list(adapt_delta = 0.99),
  file = here::here("rdata/clr")
)
pp_check(clr)
summary(clr)




# infancy 
clr_if <- brm(
  family = gaussian(),
  formula = bmi ~ bmi_tminus + Akkermansia + Alistipes + Anaerostipes + Bacteroides + Bifidobacterium + Blautia + Coprococcus_1 + Coprococcus_2 + Coprococcus_3 + Dialister + Faecalibacterium + Holdemanella + Phascolarctobacterium + Prevotella_2 + Prevotella_7 + Prevotella_9 + Roseburia + Subdoligranulum + s_bw_grams + Eubacterium_hallii_group + (1 | id),
  prior = priors2,
  data = filter(df, time %in% c(1, 2, 3)),
  control = list(adapt_delta = 0.99),
  file = here::here("rdata/clr_if")
)
pp_check(clr_if)
summary(clr_if)




clr_ch <- brm(
  family = gaussian(),
  formula = bmi ~ bmi_tminus + Akkermansia + Alistipes + Anaerostipes + Bacteroides + Bifidobacterium + Blautia + Coprococcus_1 + Coprococcus_2 + Coprococcus_3 + Dialister + Faecalibacterium + Holdemanella + Phascolarctobacterium + Prevotella_2 + Prevotella_7 + Prevotella_9 + Roseburia + Subdoligranulum + s_bw_grams + Eubacterium_hallii_group + (1 | id),
  prior = priors2,
  data = filter(df, time %in% c(4, 5)),
  control = list(adapt_delta = 0.99),
  file = here::here("rdata/clr_ch"),
  iter = 2e3
)
pp_check(clr_ch)
summary(clr_ch)











scfa_plots <- map2(list(clr_imp, clr_if_imp, clr_ch_imp), c(1, 2, 3), function(model, no){
  df <- posterior_samples(model) %>% select(
    b_Akkermansia,
    b_Alistipes,
    b_Anaerostipes,
    b_Bacteroides,
    b_Bifidobacterium,
    b_Blautia,
    b_Coprococcus_1,
    b_Coprococcus_2,
    b_Coprococcus_3,
    b_Dialister,
    b_Faecalibacterium,
    b_Holdemanella,
    b_Phascolarctobacterium,
    b_Prevotella_2,
    b_Prevotella_7,
    b_Prevotella_9,
    b_Roseburia,
    b_Subdoligranulum,
    b_Eubacterium_hallii_group) %>%
    pivot_longer(everything(), names_to = "Genus") %>%
    mutate(Genus = str_replace(Genus, "b_", "")) %>%
    group_by(Genus) %>%
    summarise(Estimate = median(value), lower = get_lower(value), upper = get_upper(value), .groups = "drop")
  p <- df %>%
    ggplot(aes(Genus, Estimate)) +
      geom_errorbar(aes(ymin = lower, ymax = upper)) +
      geom_point() +
      geom_hline(aes(yintercept = 0), linetype = "dashed") +
      theme_bw(base_size = 30) +
      ylim(c(-0.45, 0.3)) 
    
    if (no > 1) {
      p <- p + theme(axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank())
    }
    if (no != 2) {
      p <- p + theme(axis.title.x=element_blank())
    }
    
    if (no == 1) {
      p <- p + scale_x_discrete(labels = c(
                              expression(bolditalic("Akkermansia")),
                              expression(bolditalic("Alistipes")),
                              expression(italic("Anaerostipes")),
                              expression(italic("Bacteroides")),
                              expression(italic("Bifidobacterium")),
                              expression(italic("Blautia")),
                              expression(italic("Coprococcus_1")),
                              expression(italic("Coprococcus_2")),
                              expression(italic("Coprococcus_3")),
                              expression(italic("Dialister")),
                              expression(italic("Eubacterium_hallii_group")),
                              expression(italic("Faecalibacterium")),
                              expression(italic("Holdemanella")),
                              expression(italic("Phascolarctobacterium")),
                              expression(italic("Prevotella_2")),
                              expression(italic("Prevotella_7")),
                              expression(italic("Prevotella_9")),
                              expression(italic("Roseburia")),
                              expression(bolditalic("Subdoligranulum")))) 
      df2 <- filter(df, Genus %in% c("Alistipes", "Subdoligranulum"))
      p <- p + 
            geom_errorbar(
              data = df2, 
              aes(x = Genus, y = Estimate, ymin = lower, ymax = upper), 
              size = 2) +
            geom_point(data = df2, aes(x = Genus, y = Estimate), size = 5) 
        
    }
    
    if (no == 2) {
      
      df2 <- filter(df, Genus %in% c("Alistipes", "Subdoligranulum", "Akkermansia"))
      p <- p + 
            geom_errorbar(
              data = df2, 
              aes(x = Genus, y = Estimate, ymin = lower, ymax = upper), 
              size = 2) +
            geom_point(data = df2, aes(x = Genus, y = Estimate), size = 5)
        
    }
    
    if (no == 3) {
      
      df2 <- filter(df, Genus %in% c("Alistipes", "Subdoligranulum"))
      p <- p + 
            geom_errorbar(
              data = df2, 
              aes(x = Genus, y = Estimate, ymin = lower, ymax = upper), 
              size = 2) +
            geom_point(data = df2, aes(x = Genus, y = Estimate), size = 5)
        
    }
    
    
    
    p + coord_flip()
})
 


scfa_plots
save(scfa_plots, file = here::here("rdata/publish.Rds"))















#####################################################################
##########        create tables for publication            ##########
#####################################################################

reported <- list(
  fb,
  
  fbsep_mi,
  fbsep_if_mi,
  fbsep_ch_mi,
  
  mscfa1_mi,
  mscfa1_if_mi,
  mscfa1_ch_mi,
  clr_imp,
  clr_if_imp,
  clr_ch_imp
)
addons <- c(
  "child samples", 
  
  "all samples",
  "infant samples", 
  "child samples",
  
  "all samples",
  "infant samples",
  "child samples",
  
  "all_samples", 
  "infant samples",
  "child samples")
  
table_nr <- 1:10

# print header followed by table for each model 
file.create(here::here("article/table_supplement.Rmd"))
write_lines(
  '---\noutput:\n  html_document:\n    toc: true\n---', file = here::here("article/table_supplement.Rmd"), append = FALSE)
# }

pmap(list(reported, addons, table_nr), function(model, addon, nr) {
  write_lines(
    x = glue("\n\n### Supplementary Table {nr}\n \n"),
    file = here::here("article/table_supplement.Rmd"), append = TRUE)
  
  write_lines(
    x = knitr::kable(
      brms_tbl(model, exclude_id = TRUE), # ID must be excluded to guarantee anonymity
      caption = glue("Coefficients for the model using {addon} and formula: {as.character(model$formula)[1]}\n \n")
      ),
    file = here::here("article/table_supplement.Rmd"), append = TRUE)
})



# NOTE: random_forests.R has to be run to complete this table.
if(!file.exists(here::here("rdata/rf_hp_tbl.Rds"))) {
  print("To get the complete table, you need to first complete random_forests.R")
 } else {
  load(here::here("rdata/rf_hp_tbl.Rds"))
  write_lines(
    x = glue("\n\n### Supplementary Table 11"),
    file = here::here("article/table_supplement.Rmd"), append = TRUE)

  write_lines(
    x = knitr::kable(rf_hp_tbl, caption = "Random Forest Hyperparameters Per Model"),
    file = here::here("article/table_supplement.Rmd"), append = TRUE)
  rmarkdown::render(here::here("article/table_supplement.Rmd"))
}








