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

## Load admin names -----------------------------------------------

setwd(data.dir)

load(paste0(poly.path,'/', country, '_Amat.rda'))
load(paste0(poly.path,'/', country, '_Amat_Names.rda'))

## Load data -----------------------------------------------

load(paste0(country,'_cluster_dat_APC.rda'), envir = .GlobalEnv)

survey_years <- unique(mod.dat$survey)
beg.year <- 2006
end.year <- 2013 
end.proj.year <- end.year + 5

births.2014 <- 
  mod.dat %>% 
  dplyr::filter(survey == survey_years[3], years %in% beg.year:end.year) %>% 
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

if(dir.exists(paths = paste0(data.dir, '/worldpop'))){
  setwd(paste0(data.dir, '/worldpop'))
  load(file = 'adm1_weights_u5.rda')
  
  weight.strata.adm1.u5.natl <- 
    weight.strata.adm1.u5 %>% 
    dplyr::left_join(., weight.adm1.u5, by = c('region', 'years')) %>% 
    dplyr::mutate(urban = urban*proportion,
                  rural = rural*proportion) %>% 
    dplyr::select(-proportion)
  
}

## Fit APC models  -----------------------------------------------

setwd(paste0(res.dir))

### National -----------------------------------------------

#### APC -----------------------------------------------

apc.natl.u5.fit <- smoothAPC(data = births.2014,
                             Amat = NULL,
                             mod = 'apc',
                             admin.level = 'National',
                             year.label = c(beg.year:end.proj.year),
                             slope.drop = 'c',
                             stratified = TRUE,
                             overdisp.mean = -7.5,
                             overdisp.prec = 0.39)

apc.natl.u5.res <- getSmoothedAPC(inla.mod = apc.natl.u5.fit,
                                  n.sim = 1000,
                                  save.u5m.draws = TRUE,
                                  save.linpred.draws = TRUE,
                                  strata.weights = weight.strata.natl.u5)

save(apc.natl.u5.fit, file = paste0('Age-Period-Cohort/U5MR/',country,'_fit_natl_u5.rda'))
save(apc.natl.u5.res, file = paste0('Age-Period-Cohort/U5MR/',country,'_res_natl_u5.rda'))

#### AP -----------------------------------------------

ap.natl.u5.fit <- smoothAPC(data = births.2014,
                            Amat = NULL,
                            mod = 'ap',
                            admin.level = 'National',
                            year.label = c(beg.year:end.proj.year),
                            stratified = TRUE,
                            overdisp.mean = -7.5,
                            overdisp.prec = 0.39)

ap.natl.u5.res <- getSmoothedAPC(inla.mod = ap.natl.u5.fit,
                                 n.sim = 1000,
                                 save.u5m.draws = TRUE,
                                 save.linpred.draws = TRUE,
                                 strata.weights = weight.strata.natl.u5)

save(ap.natl.u5.fit, file = paste0('Age-Period/U5MR/',country,'_fit_natl_u5.rda'))
save(ap.natl.u5.res, file = paste0('Age-Period/U5MR/',country,'_res_natl_u5.rda'))

#### AC -----------------------------------------------

ac.natl.u5.fit <- smoothAPC(data = births.2014,
                            Amat = NULL,
                            mod = 'ac',
                            admin.level = 'National',
                            year.label = c(beg.year:end.proj.year),
                            stratified = TRUE,
                            overdisp.mean = -7.5,
                            overdisp.prec = 0.39)

ac.natl.u5.res <- getSmoothedAPC(inla.mod = ac.natl.u5.fit,
                                 n.sim = 1000,
                                 save.u5m.draws = TRUE,
                                 save.linpred.draws = TRUE,
                                 strata.weights = weight.strata.natl.u5)

save(ac.natl.u5.fit, file = paste0('Age-Cohort/U5MR/',country,'_fit_natl_u5.rda'))
save(ac.natl.u5.res, file = paste0('Age-Cohort/U5MR/',country,'_res_natl_u5.rda'))

### Admin1 -----------------------------------------------

#### APC -----------------------------------------------

tryCatch({
  apc.adm1.u5.fit <- smoothAPC(data = births.2014,
                               Amat = admin1.mat,
                               mod = 'apc',
                               admin.level = 'Admin1',
                               year.label = c(beg.year:end.proj.year),
                               slope.drop = 'c',
                               stratified = TRUE,
                               overdisp.mean = -7.5,
                               overdisp.prec = 0.39)
  
  apc.adm1.u5.res <- getSmoothedAPC(inla.mod = apc.adm1.u5.fit,
                                    n.sim = 1000,
                                    save.u5m.draws = TRUE,
                                    save.linpred.draws = TRUE,
                                    strata.weights = weight.strata.adm1.u5)
  
  
  apc.adm1.u5.natl.res <- getSmoothedAPC(inla.mod = apc.adm1.u5.fit,
                                         n.sim = 1000,
                                         save.u5m.draws = TRUE,
                                         save.linpred.draws = TRUE,
                                         strata.weights = weight.strata.adm1.u5.natl)
  
  save(apc.adm1.u5.fit, file = paste0('Age-Period-Cohort/U5MR/',country,'_fit_adm1_u5.rda'))
  save(apc.adm1.u5.res, file = paste0('Age-Period-Cohort/U5MR/',country,'_res_adm1_u5.rda'))
  save(apc.adm1.u5.natl.res, file = paste0('Age-Period-Cohort/U5MR/',country,'_res_adm1_u5_natl.rda'))
  
}, silent=T, error = function(e) {message('\nYearly Age-Period-Cohort model cannot be fit at the Admin1 level due to data sparsity.')})

#### AP -----------------------------------------------

tryCatch({
  ap.adm1.u5.fit <- smoothAPC(data = births.2014,
                              Amat = admin1.mat,
                              mod = 'ap',
                              admin.level = 'Admin1',
                              year.label = c(beg.year:end.proj.year),
                              stratified = TRUE,
                              overdisp.mean = -7.5,
                              overdisp.prec = 0.39)
  
  ap.adm1.u5.res <- getSmoothedAPC(inla.mod = ap.adm1.u5.fit,
                                   n.sim = 1000,
                                   save.u5m.draws = TRUE,
                                   save.linpred.draws = TRUE,
                                   strata.weights = weight.strata.adm1.u5)
  
  ap.adm1.u5.natl.res <- getSmoothedAPC(inla.mod = ap.adm1.u5.fit,
                                        n.sim = 1000,
                                        save.u5m.draws = TRUE,
                                        save.linpred.draws = TRUE,
                                        strata.weights = weight.strata.adm1.u5.natl)
  
  save(ap.adm1.u5.fit, file = paste0('Age-Period/U5MR/',country,'_fit_adm1_u5.rda'))
  save(ap.adm1.u5.res, file = paste0('Age-Period/U5MR/',country,'_res_adm1_u5.rda'))
  save(ap.adm1.u5.natl.res, file = paste0('Age-Period/U5MR/',country,'_res_adm1_u5_natl.rda'))
  
}, silent=T, error = function(e) {message('\nYearly Age-Period model cannot be fit at the Admin1 level due to data sparsity.')})



#### AC -----------------------------------------------

tryCatch({
  ac.adm1.u5.fit <- smoothAPC(data = births.2014,
                              Amat = admin1.mat,
                              mod = 'ac',
                              admin.level = 'Admin1',
                              year.label = c(beg.year:end.proj.year),
                              stratified = TRUE,
                              overdisp.mean = -7.5,
                              overdisp.prec = 0.39)
  
  ac.adm1.u5.res <- getSmoothedAPC(inla.mod = ac.adm1.u5.fit,
                                   n.sim = 1000,
                                   save.u5m.draws = TRUE,
                                   save.linpred.draws = TRUE,
                                   strata.weights = weight.strata.adm1.u5)
  
  ac.adm1.u5.natl.res <- getSmoothedAPC(inla.mod = ac.adm1.u5.fit,
                                        n.sim = 1000,
                                        save.u5m.draws = TRUE,
                                        save.linpred.draws = TRUE,
                                        strata.weights = weight.strata.adm1.u5.natl)
  
  save(ac.adm1.u5.fit, file = paste0('Age-Cohort/U5MR/',country,'_fit_adm1_u5.rda'))
  save(ac.adm1.u5.res, file = paste0('Age-Cohort/U5MR/',country,'_res_adm1_u5.rda'))
  save(ac.adm1.u5.natl.res, file = paste0('Age-Cohort/U5MR/',country,'_res_adm1_u5_natl.rda'))
}, silent=T, error = function(e) {message('\nYearly Age-Cohort model cannot be fit at the Admin1 level due to data sparsity.')})