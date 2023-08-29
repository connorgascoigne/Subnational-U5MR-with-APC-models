rm(list = ls())
## ENTER COUNTRY OF INTEREST -----------------------------------------------
# Please capitalize the first letter of the country name and replace " " in the country name to "_" if there is.
country <- 'Kenya'

## Libraries -----------------------------------------------

library(SUMMER)
library(INLA)
library(tidyverse)

## Retrieve directories and country info -----------------------------------------------
code.path <- rstudioapi::getActiveDocumentContext()$path
code.path.splitted <- strsplit(code.path, "/")[[1]]

home.dir <- paste(code.path.splitted[1: (length(code.path.splitted)-2)], collapse = "/")
data.dir <- paste0(home.dir,'/Data/',country) # set the directory to store the data
res.dir <- paste0(home.dir,'/Results/',country) # set the directory to store the results (e.g. fitted R objects, figures, tables in .csv etc.)
info.name <- paste0(country, "_general_info.Rdata")
load(file = paste0(home.dir,'/Info/',info.name, sep='')) # load the country info

source(file = paste0(home.dir, '/Code/functions.R'))

# Load polygon files  ------------------------------------------------------

setwd(data.dir)

poly.adm0 <- sf::st_read(dsn = poly.path, layer = as.character(poly.layer.adm0)) # load the national shape file
poly.adm1 <- sf::st_read(dsn = poly.path, layer = as.character(poly.layer.adm1)) # load the shape file of admin-1 regions

load(paste0(poly.path,'/', country, '_Amat.rda'))
load(paste0(poly.path,'/', country, '_Amat_Names.rda'))

## Load admin names -----------------------------------------------
setwd(data.dir)

load(paste0(poly.path,'/', country, '_Amat.rda'))
load(paste0(poly.path,'/', country, '_Amat_Names.rda'))

## Load data -----------------------------------------------

load(paste0(country,'_cluster_dat_APC.rda'), envir = .GlobalEnv)

survey_years <- unique(mod.dat$survey)
beg.year <- 2006
end.year.cv <- 2012
end.proj.year.cv <- end.year.cv + 1

births.2014 <- 
  mod.dat %>% 
  dplyr::filter(survey == survey_years[3], years %in% beg.year:end.proj.year.cv) %>% 
  dplyr::mutate(years = years %>% as.character() %>% as.numeric(),
                v005 = v005/1e6,
                died = Y,
                total = total %>% as.numeric())

## Load UR proportions -----------------------------------------------

if(dir.exists(paths = paste0(res.dir,'/UR/'))){
  setwd(paste0(res.dir,'/UR'))
  weight.strata.natl.u5 <- readRDS(paste0('U5_fraction/','natl_u5_urban_weights.rds'))
  weight.strata.natl.u5$rural <- 1-weight.strata.natl.u5$urban
  weight.strata.adm1.u5 <- readRDS(paste0('U5_fraction/','admin1_u5_urban_weights.rds'))
}

# Fit CV models  -----------------------------------------------

setwd(paste0(res.dir))

n.sim <- 1000
n.regions <- nrow(admin1.names)
apc.adm1.u5.cv <- ap.adm1.u5.cv <- ac.adm1.u5.cv <- array(NA, dim = c(n.regions, 2 + n.sim))
for(i in 1:n.regions){
  
  # i <- 1
  
  births.2014.cv <-
    births.2014 %>% 
    # cross validation on the final year included in the
    dplyr::mutate(Y = dplyr::if_else(admin1.char == admin1.names$Internal[i] & years == end.proj.year.cv, NA, Y))
  
  ### Admin1 -----------------------------------------------
  
  #### APC -----------------------------------------------
  
  tryCatch({
    apc.adm1.u5.fit <- 
      smoothAPC(data = births.2014.cv,
                Amat = admin1.mat,
                mod = 'apc',
                admin.level = 'Admin1',
                year.label = c(beg.year:end.proj.year.cv),
                slope.drop = 'c',
                stratified = TRUE,
                overdisp.mean = -7.5,
                overdisp.prec = 0.39) %>% 
      suppressWarnings() %>% 
      suppressMessages()
    
    apc.adm1.u5.res <- 
      getSmoothedAPC(inla.mod = apc.adm1.u5.fit,
                     n.sim = n.sim,
                     strata.weights = weight.strata.adm1.u5,
                     save.u5m.draws = TRUE) %>% 
      suppressWarnings() %>% 
      suppressMessages()
    
    apc.adm1.u5.cv[i,] <- 
      apc.adm1.u5.res$u5m.draws %>% 
      dplyr::filter(region == admin1.names$Internal[i], years == end.proj.year.cv) %>%
      dplyr::mutate(across(starts_with('theta:'), ~ .x * proportion)) %>%
      dplyr::summarise(across(starts_with('theta:'), sum),
                       .by = c('region', 'years')) %>% 
      as.matrix()
    
  }, silent=T, error = function(e) {message('\nYearly Age-Period-Cohort model cannot be fit at the Admin1 level due to data sparsity.')})
  
  #### AP -----------------------------------------------
  
  tryCatch({
    ap.adm1.u5.fit <- 
      smoothAPC(data = births.2014.cv,
                Amat = admin1.mat,
                mod = 'ap',
                admin.level = 'Admin1',
                year.label = c(beg.year:end.proj.year.cv),
                stratified = TRUE,
                overdisp.mean = -7.5,
                overdisp.prec = 0.39) %>% 
      suppressWarnings() %>% 
      suppressMessages()
    
    ap.adm1.u5.res <- 
      getSmoothedAPC(inla.mod = ap.adm1.u5.fit,
                     n.sim = n.sim,
                     strata.weights = weight.strata.adm1.u5,
                     save.u5m.draws = TRUE) %>% 
      suppressWarnings() %>% 
      suppressMessages()
    
    ap.adm1.u5.cv[i,] <- 
      ap.adm1.u5.res$u5m.draws %>% 
      dplyr::filter(region == admin1.names$Internal[i], years == end.proj.year.cv) %>%
      dplyr::mutate(across(starts_with('theta:'), ~ .x * proportion)) %>%
      dplyr::summarise(across(starts_with('theta:'), sum),
                       .by = c('region', 'years')) %>% 
      as.matrix()
    
  }, silent=T, error = function(e) {message('\nYearly Age-Period model cannot be fit at the Admin1 level due to data sparsity.')})
  
  #### AC -----------------------------------------------
  
  tryCatch({
    ac.adm1.u5.fit <- 
      smoothAPC(data = births.2014.cv,
                Amat = admin1.mat,
                mod = 'ac',
                admin.level = 'Admin1',
                year.label = c(beg.year:end.proj.year.cv),
                stratified = TRUE,
                overdisp.mean = -7.5,
                overdisp.prec = 0.39) %>% 
      suppressWarnings() %>% 
      suppressMessages()
    
    ac.adm1.u5.res <- 
      getSmoothedAPC(inla.mod = ac.adm1.u5.fit,
                     n.sim = n.sim,
                     strata.weights = weight.strata.adm1.u5,
                     save.u5m.draws = TRUE) %>% 
      suppressWarnings() %>% 
      suppressMessages()
    
    ac.adm1.u5.cv[i,] <- 
      ac.adm1.u5.res$u5m.draws %>% 
      dplyr::filter(region == admin1.names$Internal[i], years == end.proj.year.cv) %>%
      dplyr::mutate(across(starts_with('theta:'), ~ .x * proportion)) %>%
      dplyr::summarise(across(starts_with('theta:'), sum),
                       .by = c('region', 'years')) %>% 
      as.matrix()
    
  }, silent=T, error = function(e) {message('\nYearly Age-Cohort model cannot be fit at the Admin1 level due to data sparsity.')})
  
  print(paste0(admin1.names$GADM[i], ' (', admin1.names$Internal[i], ') complete'))
  
}

colnames(apc.adm1.u5.cv) <-
  colnames(ap.adm1.u5.cv) <-
  colnames(ac.adm1.u5.cv) <-
  c('region', 'years', paste0('theta:', 1:n.sim))

# add sample distribution variation -----------------------------------------------

setwd(res.dir)

# direct estimates
load(file = paste0('Direct/U5MR/',country, '_direct_adm1_yearly_u5.rda'))

direct.admin1.yearly.u5.2013 <- 
  direct.admin1.yearly.u5 %>% 
  dplyr::filter(region != 'All', years == 2013) %>% 
  dplyr::select(region, years, logit.est, logit.prec)

apc.adm1.u5.cv <- 
  apc.adm1.u5.cv %>% 
  as.data.frame() %>% 
  dplyr::mutate(years = years %>% as.numeric(),
                dplyr::across(dplyr::starts_with('theta:'), ~ SUMMER::logit(as.numeric(.x)))) %>% 
  dplyr::left_join(., direct.admin1.yearly.u5.2013, by = c('region', 'years')) %>% 
  dplyr::mutate(dplyr::across(dplyr::starts_with('theta:'), ~ .x + rnorm(n = 1, mean = 0, sd = 1/sqrt(logit.prec)))) %>% 
  dplyr::mutate(dplyr::select(., dplyr::starts_with('theta:')) %>% 
                  apply(., 1, quantileSummary) %>% 
                  lapply(., data.frame) %>% 
                  do.call(rbind, .)) %>% 
  dplyr::relocate(dplyr::starts_with('theta:'), .after = dplyr::last_col())

ap.adm1.u5.cv <- 
  ap.adm1.u5.cv %>% 
  as.data.frame() %>% 
  dplyr::mutate(years = years %>% as.numeric(),
                dplyr::across(dplyr::starts_with('theta:'), ~ SUMMER::logit(as.numeric(.x)))) %>% 
  dplyr::left_join(., direct.admin1.yearly.u5.2013, by = c('region', 'years')) %>% 
  dplyr::mutate(dplyr::across(dplyr::starts_with('theta:'), ~ .x + rnorm(n = 1, mean = 0, sd = 1/sqrt(logit.prec)))) %>% 
  dplyr::mutate(dplyr::select(., dplyr::starts_with('theta:')) %>% 
                  apply(., 1, quantileSummary) %>% 
                  lapply(., data.frame) %>% 
                  do.call(rbind, .)) %>% 
  dplyr::relocate(dplyr::starts_with('theta:'), .after = dplyr::last_col())

ac.adm1.u5.cv <- 
  ac.adm1.u5.cv %>% 
  as.data.frame() %>% 
  dplyr::mutate(years = years %>% as.numeric(),
                dplyr::across(dplyr::starts_with('theta:'), ~ SUMMER::logit(as.numeric(.x)))) %>% 
  dplyr::left_join(., direct.admin1.yearly.u5.2013, by = c('region', 'years')) %>% 
  dplyr::mutate(dplyr::across(dplyr::starts_with('theta:'), ~ .x + rnorm(n = 1, mean = 0, sd = 1/sqrt(logit.prec)))) %>% 
  dplyr::mutate(dplyr::select(., dplyr::starts_with('theta:')) %>% 
                  apply(., 1, quantileSummary) %>% 
                  lapply(., data.frame) %>% 
                  do.call(rbind, .)) %>% 
  dplyr::relocate(dplyr::starts_with('theta:'), .after = dplyr::last_col())

save(apc.adm1.u5.cv, file = paste0('Age-Period-Cohort/U5MR/',country,'_cv_adm1_u5.rda'))
save(ap.adm1.u5.cv, file = paste0('Age-Period/U5MR/',country,'_cv_adm1_u5.rda'))
save(ac.adm1.u5.cv, file = paste0('Age-Cohort/U5MR/',country,'_cv_adm1_u5.rda'))
