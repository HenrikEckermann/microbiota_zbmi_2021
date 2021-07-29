library(tidyverse)
# load helper scrips
source(here::here("R/ml_helper.R"))
source(here::here("R/mb_helper.R"))


#####################################################################
##########           import and prepare data               ##########
#####################################################################

load(here::here("rdata/mb_import_H.Rds"))
 


#####################################################################
##########               Genus Level RF                    ##########
#####################################################################



otu_rf <- otu_rf %>% mutate(subject = str_sub(Sample, 2, 4))
df_tm <- readxl::read_xlsx(here::here("data/tm_df_1.xlsx")) %>%
  select(subject, contains("bmi")) %>%
  mutate(subject = as.character(subject))

# for each time point I fit a model that associates the cross sectional data 
# plus one model where we predict the later BMI time point 
times <- c("1m", "3m", "4m", "6y", "10y")
bmi_pairs <- list(
   c("bmiAgeZ_1mo", "bmiAgeZ_3mo"),
   c("bmiAgeZ_3mo", "bmiAgeZ_4mo"),
   c("bmiAgeZ_4mo", "bmiAgeZ_25y"),
   c("bmiAgeZ_6y", "bmiAgeZ_7y"),
   c("bmiAgeZ_10y", "bmiAgeZ_125y")
)

rf_genus <- map2(bmi_pairs, times, function(Y, time_var) {
  
  # create model df
  d <- otu_rf %>% filter(time == time_var) %>% 
    full_join(
      select(df_tm, subject, bmi_cs = all_of(Y[1]), bmi_fut = all_of(Y[2])),
      by = "subject"
    ) %>%
    na.omit() 

  X <- select(d, matches("^\\d+")) %>% colnames()
  y_cs <- "bmi_cs"
  y_fut <- "bmi_fut"




  if (!file.exists(here::here(glue("rdata/rf0_{time_var}_{y_cs}_cs_genus.Rds")))) {
    model_and_data_cs <- future_map(1:10, function(rep) {
      rf_cv(
        data = d,
        features = X,
        y = y_cs,
        k = 4,
        ntree = 5000
      )
    }, .options = furrr_options(seed = TRUE)) %>% flatten()
    
    save(model_and_data_cs, file = here::here(glue("rdata/rf0_{time_var}_{y_cs}_cs_genus.Rds")))
    
  } else {
    load(here::here(glue("rdata/rf0_{time_var}_{y_cs}_cs_genus.Rds")))
  }

  # tuning 
  if (!file.exists(here::here(glue("rdata/tr0_{time_var}_{y_cs}_cs_genus.Rds")))) {
    pars_cs <- tune_rf(
      data = d,
      features = X,
      y = y_cs,
      regression = TRUE,
      ntree = 5000,
      tune.parameters = c("mtry", "sample.fraction")
    )
    
    save(pars_cs, file = here::here(glue("rdata/tr0_{time_var}_{y_cs}_cs_genus.Rds")))
    
  } else {
    load(here::here(glue("rdata/tr0_{time_var}_{y_cs}_cs_genus.Rds")))
  }
  
  # model with tuned pars 
  if (!file.exists(here::here(glue("rdata/rf1_{time_var}_{y_cs}_cs_genus.Rds")))) {
    model_and_data2_cs <- future_map(1:10, function(rep) {
      rf_cv(
        data = d,
        features = X,
        y = y_cs,
        k = 4,
        ntree = 5000,
        regression = TRUE,
        mtry = pars_cs$recommended.pars[1, "mtry"],
        sample.fraction = pars_cs$recommended.pars[1, "sample.fraction"] 
      )}, .options = furrr_options(seed = TRUE)) %>% flatten()
      
      save(model_and_data2_cs, file = here::here(glue("rdata/rf1_{time_var}_{y_cs}_cs_genus.Rds")))
      
  } else {
    load(here::here(glue("rdata/rf1_{time_var}_{y_cs}_cs_genus.Rds")))
  }
  
  # now for future bmi 
  if (!file.exists(here::here(glue("rdata/rf0_{time_var}_{y_fut}_fut_genus.Rds")))) {
    model_and_data_fut <- future_map(1:10, function(rep) {
      rf_cv(
        data = d,
        features = X,
        y = y_fut,
        k = 4,
        ntree = 5000
      )
    }, .options = furrr_options(seed = TRUE)) %>% flatten()
    
    save(model_and_data_fut, file = here::here(glue("rdata/rf0_{time_var}_{y_fut}_fut_genus.Rds")))
    
  } else {
    load(here::here(glue("rdata/rf0_{time_var}_{y_fut}_fut_genus.Rds")))
  }


  # tuning 
  if (!file.exists(here::here(glue("rdata/tr0_{time_var}_{y_fut}_fut_genus.Rds")))) {
    pars_fut <- tune_rf(
      data = d,
      features = X,
      y = y_fut,
      regression = TRUE,
      ntree = 5000,
      tune.parameters = c("mtry", "sample.fraction")
    )
    
    save(pars_fut, file = here::here(glue("rdata/tr0_{time_var}_{y_fut}_fut_genus.Rds")))
    
  } else {
    load(here::here(glue("rdata/tr0_{time_var}_{y_fut}_fut_genus.Rds")))
  }
  
  # model with tuned pars 
  if (!file.exists(here::here(glue("rdata/rf1_{time_var}_{y_fut}_fut_genus.Rds")))) {
    model_and_data2_fut <- future_map(1:10, function(rep) {
      rf_cv(
        data = d,
        features = X,
        y = y_fut,
        k = 4,
        ntree = 5000,
        regression = TRUE,
        mtry = pars_fut$recommended.pars[1, "mtry"],
        sample.fraction = pars_fut$recommended.pars[1, "sample.fraction"]
      )}, .options = furrr_options(seed = TRUE)) %>% flatten()
      
      save(model_and_data2_fut, file = here::here(glue("rdata/rf1_{time_var}_{y_fut}_fut_genus.Rds")))
      
  } else {
    load(here::here(glue("rdata/rf1_{time_var}_{y_fut}_fut_genus.Rds")))
  }
  
  
  list(
    pars_cs = pars_cs$recommended.pars,
    pars_fut = pars_fut$recommended.pars,
    models_cs = model_and_data_cs,
    models_fut = model_and_data_fut,
    models_cs2 = model_and_data2_cs,
    models_fut2 = model_and_data2_fut
  )
})

# create table with the hyperparameters for supplement  
rf_hp_tbl <- pmap_dfr(list(bmi_pairs, times, rf_genus), function(Y, time_var, mlist) {
  bind_rows(
    mlist[[1]] %>% mutate(time = time_var, y = Y[1]),
    mlist[[2]] %>% mutate(time = time_var, y = Y[2])
  ) %>%
  mutate(y = str_replace(y, "bmiAgeZ_", "")) %>%
  select(zBMI = y, Microbiota = time, Mtry = mtry, `Sample Fraction` = sample.fraction)
})
rf_hp_tbl
save(rf_hp_tbl, file = here::here("rdata/rf_hp_tbl.Rds"))










# obtain null distribution ---------------------------
nperm <- 1000
nulldist <- map_dfr(1:nperm, function(nulliter) {
  # for cross sectional and future obtain 1 null value per nperm
  if (!file.exists(here::here(glue("rdata/perm{nulliter}.Rds")))) {
    rf_genus_null <- map2_dfr(bmi_pairs, times, function(Y, time_var) {

        # create model df
        d <- otu_rf %>% filter(time == time_var) %>% 
          full_join(
            select(df_tm, subject, bmi_cs = all_of(Y[1]), bmi_fut = all_of(Y[2])),
            by = "subject"
          ) %>%
          na.omit() %>%
          mutate(
            bmi_cs = sample(bmi_cs, replace = FALSE), 
            bmi_fut = sample(bmi_fut, replace = FALSE))

        X <- select(d, matches("^\\d+")) %>% colnames()
        y_cs <- "bmi_cs"
        y_fut <- "bmi_fut"
        

          model_and_data_cs <- rf_cv(
            data = d,
            features = X,
            y = y_cs,
            k = 4,
            ntree = 5000
          )


          pars_cs <- tune_rf(
            data = d,
            features = X,
            y = y_cs,
            regression = TRUE,
            ntree = 5000,
            tune.parameters = c("mtry", "sample.fraction")
          )

          model_and_data2_cs <- rf_cv(
            data = d,
            features = X,
            y = y_cs,
            k = 4,
            ntree = 5000,
            regression = TRUE,
            mtry = pars_cs$recommended.pars[1, "mtry"],
            sample.fraction = pars_cs$recommended.pars[1, "sample.fraction"] 
            )

          model_and_data_fut <- rf_cv(
            data = d,
            features = X,
            y = y_fut,
            k = 4,
            ntree = 5000
          )
          

          pars_fut <- tune_rf(
            data = d,
            features = X,
            y = y_fut,
            regression = TRUE,
            ntree = 5000,
            tune.parameters = c("mtry", "sample.fraction")
          )
          

          model_and_data2_fut <- rf_cv(
            data = d,
            features = X,
            y = y_fut,
            k = 4,
            ntree = 5000,
            regression = TRUE,
            mtry = pars_fut$recommended.pars[1, "mtry"],
            sample.fraction = pars_fut$recommended.pars[1, "sample.fraction"]
          )
            
        
        oob_cs <- get_oob(model_and_data2_cs) %>% .$median
        oob_fut <- get_oob(model_and_data2_fut) %>% .$median
        pearson_cs <- get_pearson(model_and_data2_cs, y_cs) %>% .$median
        pearson_fut <- get_pearson(model_and_data2_fut, y_fut) %>% .$median
        
        
        out <- list(
          iter = nulliter,
          time = time_var,
          "oob_cs" = oob_cs,
          "oob_fut" = oob_fut,
          "pearson_cs" = pearson_cs,
          "pearson_fut" = pearson_fut      
        )
      out
    })
    save(rf_genus_null, file = here::here(glue("rdata/perm{nulliter}.Rds")))
  } else {
    load(here::here(glue("rdata/perm{nulliter}.Rds")))
  }
  rf_genus_null
})

# Now for each median correlation/oob estimate, we need to see how likely
# it would occur in the null-distribution to get the p value 

 # interpret models ---------------------------

genus_metric <- map2_dfr(rf_genus, 1:5, function(ls, i) {
  time_var = times[i]
  ordering <- c("cs0", "cs1", "fut0", "fut1")

  oob1 <- get_oob(ls[["models_cs"]])
  oob2 <- get_oob(ls[["models_cs2"]])
  oob3 <- get_oob(ls[["models_fut"]])
  oob4 <- get_oob(ls[["models_fut2"]])
    
    
  pearson1 <- get_pearson(ls[["models_cs"]], "bmi_cs")
  pearson2 <- get_pearson(ls[["models_cs2"]], "bmi_cs")
  pearson3 <- get_pearson(ls[["models_fut"]], "bmi_fut")
  pearson4 <- get_pearson(ls[["models_fut2"]], "bmi_fut")
  
  df1 <- bind_rows(oob1, oob2, oob3, oob4) %>%
    mutate(
      order = ordering, 
      time = time_var, 
      y = rep(c("bmi_cs", "bmi_fut"), each = 2),
      statistic = "OOB"
    ) %>%
    select(statistic, everything())
  df2 <- bind_rows(pearson1, pearson2, pearson3, pearson4) %>%
    mutate(
      order = ordering, 
      time = time_var, 
      y = rep(c("bmi_cs", "bmi_fut"), each = 2),
      statistic = "pearson"
    ) %>% 
    select(statistic, everything())

  bind_rows(df1, df2)
})


null_time <- nulldist %>% group_by(time) %>% nest()
pvalues <- genus_metric %>% filter(order %in% c("cs1", "fut1")) %>%
  select(statistic, median, order, time) %>%
  pmap(function(statistic, median, order, time){
    ind <- which(null_time == time)
    var <- ifelse(order == "cs1", glue("{str_to_lower(statistic)}_cs"), glue("{str_to_lower(statistic)}_fut"))
    if (statistic == "OOB") {
      p <- mean(median > null_time$data[[ind]][[var]])
    } else {
      p <- mean(median < null_time$data[[ind]][[var]])
    }
    
  })
pvalues <- na.omit(as.numeric(pvalues))





rf_table <- genus_metric %>% filter(order %in% c("cs1", "fut1")) %>%
  mutate(zBMI = ifelse(time == "1m" & order == "cs1", "1m", ifelse(
    time == "1m" & order == "fut1", "3m", ifelse(
      time == "3m" & order == "cs1", "3m", ifelse(
        time == "3m" & order == "fut1", "4m", ifelse(
          time == "4m" & order == "cs1", "4m", ifelse(
            time == "4m" & order == "fut1", "2.5y", ifelse(
              time == "6y" & order == "cs1", "6y", ifelse(
                time == "6y" & order == "fut1", "7y", ifelse(
                  time == "10y" & order == "cs1", "10y", ifelse(
                    time == "10y" & order == "fut1", "12.5y", NA))))))))))) %>%
  mutate(p = pvalues, q = qvalue::qvalue(p, lambda=0)$qvalues) %>%
  mutate(across(where(is.numeric), round, 3)) %>%
  mutate(
    time = factor(time, levels = c("1m", "3m", "4m", "6y", "10y")),
    zBMI = factor(zBMI, levels =c("1m", "3m", "4m", "2.5y", "6y", "7y", "10y", "12.5y")),
    statistic = ifelse(statistic != "OOB", str_to_title(statistic), statistic)) %>%
  arrange(time, zBMI) %>%
  select("Time Microbiota" = time, "Time zBMI" = zBMI, "Accuracy Measure" = statistic, Median = median, "P-value" = p, "Q-value" = q) 




#####################################################################
##########               variable importances              ##########
#####################################################################



if (!file.exists(here::here("rdata/imps.Rds"))) {
  imps <- pmap(
    list(
      time_var = c("4m", "10y", "10y"),
      y = c("bmi_cs", "bmi_cs", "bmi_fut"),
      bmi_var = c("bmiAgeZ_4mo", "bmiAgeZ_10y", "bmiAgeZ_125y")
    ),
    function(time_var, y, bmi_var) {
    
    # create model df
    d <- otu_rf %>% filter(time == time_var) %>% 
      full_join(
        select(df_tm, subject, y = all_of(bmi_var)),
        by = "subject"
      ) %>%
      na.omit() 
    colnames(d) <- make.names(colnames(d))
    X <- select(d, matches("^X\\d+")) %>% colnames()
    
    if (y == "bmi_cs") {
      if (!file.exists(here::here(glue("rdata/tr0_{time_var}_{y}_cs_genus.Rds")))) {
        pars_cs <- tune_rf(
          data = d,
          features = X,
          y = "y",
          regression = TRUE,
          ntree = 5000,
          tune.parameters = c("mtry", "sample.fraction")
        )
    
        save(pars_cs, file = here::here(glue("rdata/tr0_{time_var}_{y}_cs_genus.Rds")))
    
      } else {
        load(here::here(glue("rdata/tr0_{time_var}_{y}_cs_genus.Rds")))
      }
      pars <- pars_cs
    } else {
      if (!file.exists(here::here(glue("rdata/tr0_{time_var}_{y}_fut_genus.Rds")))) {
        pars_fut <- tune_rf(
          data = d,
          features = X,
          y = "y",
          regression = TRUE,
          ntree = 5000,
          tune.parameters = c("mtry", "sample.fraction")
        )
    
        save(pars_fut, file = here::here(glue("rdata/tr0_{time_var}_{y}_fut_genus.Rds")))
    
      } else {
        load(here::here(glue("rdata/tr0_{time_var}_{y}_fut_genus.Rds")))
      }
      pars <- pars_fut
    }
    model <- ranger(
      x = select(d, all_of(X)),
      y = d[["y"]],
      num.trees = 5000,
      importance = "permutation",
      mtry = pars$recommended.pars[1, "mtry"],
      sample.fraction = pars$recommended.pars[1, "sample.fraction"]
    )

    importance_pvalues(
      model, 
      method = "altmann",
      num.permutations = 1000,
      data = select(d, all_of(c(X, "y"))),
      formula = y ~ .
    )
  })
  save(imps, here::here("rdata/imps.Rds"))
 } else {
  load(here::here("rdata/imps.Rds"))
}



imp_tabs <- map(imps, function(imp) {
  imp_df <- imp %>%
    as.data.frame() %>%
    rownames_to_column("features") %>%
    mutate(features = str_replace(features, "X", ""))
  df_rel %>% select(features = OTU, Phylum, Class, Order, Family, Genus) %>%
    count(features, Phylum, Class, Order, Family, Genus) %>%
    filter(features %in% imp_df$features) %>%
    full_join(imp_df, by = "features") %>%
    arrange(pvalue) %>%
    mutate(across(where(is.numeric), round, 4)) %>%
    #mutate(q = qvalue::qvalue(pvalue)$qvalues) %>%
    head(10) %>%
    #arrange(desc(importance)) %>%
    mutate(
      Genus = ifelse(
        Genus == "g__", Family, ifelse(
          Genus == "g__<empty>", Order, ifelse(
            Genus == "uncultured_bacterium", Family, Genus))),
      SCFA = ifelse(Genus %in% genus_scfa, "Yes", "No")
    ) %>% 
    select(Phylum, Taxon = Genus, SCFA, Importance = importance, "P-value" = pvalue) %>%
    mutate(
      #Genus = str_replace(Genus, "g__", ""),
      Phylum = str_replace(Phylum, "p__", ""))
})


imp_tabs
save(imp_tabs, rf_table, file = here::here("rdata/rf_tables.Rds"))

