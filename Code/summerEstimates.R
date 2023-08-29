rm(list = ls())
# ENTER COUNTRY OF INTEREST -----------------------------------------------
# Please capitalize the first letter of the country name and replace ' ' in the country name to '_' if there is.
country <- 'Kenya'

# Setup
# Load libraries and info ----------------------------------------------------------

options(gsubfn.engine = 'R')
library(SUMMER)
library(tidyverse)

# extract file location of this script
code.path <- rstudioapi::getActiveDocumentContext()$path
code.path.splitted <- strsplit(code.path, '/')[[1]]

# retrieve directories
home.dir <- paste(code.path.splitted[1: (length(code.path.splitted)-2)], collapse = '/')
data.dir <- paste0(home.dir,'/Data/',country) # set the directory to store the data
res.dir <- paste0(home.dir,'/Results/',country) # set the directory to store the results (e.g. fitted R objects, figures, tables in .csv etc.)
info.name <- paste0(country, '_general_info.Rdata')
load(file = paste0(home.dir,'/Info/',info.name, sep='')) # load the country info

# Load polygon files  ------------------------------------------------------

setwd(data.dir)

load(paste0(poly.path,'/', country, '_Amat.rda'))
load(paste0(poly.path,'/', country, '_Amat_Names.rda'))
# load(paste0(country,'_cluster_dat_SUMMER.rda'), envir = .GlobalEnv)
load(paste0(country,'_cluster_dat_APC.rda'), envir = .GlobalEnv)
survey_years <- sort(unique(mod.dat$survey))

# define new beg year
beg.year <- 2006
end.year <- 2014

# Load data and separate each survey  ------------------------------------------------------

births.2014 <- 
  mod.dat %>% 
  dplyr::filter(survey == survey_years[3], years %in% beg.year:end.year) %>% 
  dplyr::mutate(years = years %>% as.character() %>% as.numeric(),
                v005 = v005/1e6,
                died = Y,
                total = total %>% as.numeric())

# Direct Estimates  ------------------------------------------------------

setwd(paste0(res.dir,'/Direct'))

## National ------------------------------------------------------

# yearly estimates
direct.natl.yearly.u5 <-  
  SUMMER::getDirect(births = births.2014, 
                    years = beg.year:end.year,
                    regionVar = 'admin1.char',
                    timeVar = 'years', 
                    clusterVar =  '~cluster',
                    ageVar = 'age', 
                    Ntrials = 'total',
                    weightsVar = 'v005',
                    national.only = T)

direct.natl.yearly.u5$survey <- 1
direct.natl.yearly.u5$surveyYears <- survey_years[3]
direct.natl.yearly.u5$region_num <- direct.natl.yearly.u5$region

# save national direct estimates
save(direct.natl.yearly.u5, file = paste0('U5MR/',country, '_direct_natl_yearly_u5.rda'))

## Admin1 ------------------------------------------------------

direct.admin1.yearly.u5 <-  
  SUMMER::getDirect(births = births.2014, 
                    years = beg.year:end.year,
                    regionVar = 'admin1.char',
                    timeVar = 'years', 
                    clusterVar =  '~cluster',
                    ageVar = 'age', 
                    Ntrials = 'total',
                    weightsVar = 'v005',
                    national.only = F)

direct.admin1.yearly.u5$survey <- 1
direct.admin1.yearly.u5$surveyYears <- survey_years[3]
direct.admin1.yearly.u5$region_num <- direct.admin1.yearly.u5$region

save(direct.admin1.yearly.u5, file = paste0('U5MR/',country, '_direct_adm1_yearly_u5.rda'))

# Smoothed direct estimates  ------------------------------------------------------

# time model
time.model <- c('rw2','ar1')[1]
# new end year for projection
end.year <- 2013
# new proj year
end.proj.year <- end.year + 5

## load in appropriate direct estimates  ------------------------------------------------------

setwd(paste0(res.dir,'/Direct'))

load(paste0('U5MR/',country, '_direct_natl_yearly_u5.rda'))
load(paste0('U5MR/',country, '_direct_adm1_yearly_u5.rda'))

## rename data ------------------------------------------------------

data.natl.yearly.u5 <- 
  direct.natl.yearly.u5 %>% 
  dplyr::filter(years %in% beg.year:end.year)
data.admin1.yearly.u5 <- 
  direct.admin1.yearly.u5 %>% 
  dplyr::filter(years %in% beg.year:end.year)

## National, yearly  ------------------------------------------------------

fit.natl.yearly.u5 <- 
  SUMMER::smoothDirect(data = data.natl.yearly.u5, 
                       geo = NULL, 
                       Amat = NULL,
                       year_label = as.character(beg.year:end.proj.year),
                       year_range = c(beg.year, end.proj.year), 
                       time.model = time.model,
                       control.inla = list(strategy = 'adaptive', int.strategy = 'auto'),  
                       is.yearly = F)
res.natl.yearly.u5 <- 
  SUMMER::getSmoothed(inla_mod = fit.natl.yearly.u5, 
                      year_range = c(2006, end.proj.year),
                      year_label = as.character(beg.year:end.proj.year))
res.natl.yearly.u5$years.num <- beg.year:end.proj.year
res.natl.yearly.u5$region.gadm <- country
save(res.natl.yearly.u5, file = paste0('U5MR/',country, '_res_natl_', time.model, '_yearly_u5_SmoothedDirect.rda')) # save the national yearly smoothed direct U5MR

## Admin 1, yearly  ------------------------------------------------------

tryCatch({
  data.admin1.yearly.u5 <- data.admin1.yearly.u5[data.admin1.yearly.u5$region!='All',]
  fit.admin1.yearly.u5 <- 
    SUMMER::smoothDirect(data = data.admin1.yearly.u5, 
                         Amat = admin1.mat, 
                         time.model = time.model,
                         year_label = as.character(beg.year:end.proj.year),
                         type.st = 4,
                         year_range = c(beg.year, end.proj.year), 
                         is.yearly = F)
  sd.admin1.yearly.u5 <- 
    SUMMER::getSmoothed(inla_mod = fit.admin1.yearly.u5,
                        Amat = admin1.mat,
                        year_label = as.character(beg.year:end.proj.year),
                        year_range = c(beg.year, end.proj.year),
                        save.draws = TRUE)
  sd.admin1.yearly.u5$region.gadm <- admin1.names$GADM[match(sd.admin1.yearly.u5$region, admin1.names$Internal)]
  save(sd.admin1.yearly.u5, file = paste0('U5MR/',country, '_res_adm1_', time.model, '_yearly_u5_SmoothedDirect.rda')) # save the admin1 yearly smoothed direct U5MR
  
 
}, silent=T, error = function(e) {message('\nYearly smoothed direct model cannot be fit at the Admin1 level due to data sparsity.')})

## smooth cluster  ------------------------------------------------------

### national  ------------------------------------------------------

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

if(dir.exists(paths = paste0(res.dir,'/UR/'))){
  setwd(paste0(res.dir,'/UR'))
  weight.strata.natl.u5 <- readRDS(paste0('U5_fraction/','natl_u5_urban_weights.rds'))
  weight.strata.natl.u5$rural <- 1-weight.strata.natl.u5$urban
  weight.strata.adm1.u5 <- readRDS(paste0('U5_fraction/','admin1_u5_urban_weights.rds'))
}

setwd(paste0(res.dir,'/Direct'))

# time model
time.model <- c('rw2','ar1')[1]
# new end year for projection
end.year <- 2013
# new proj year
end.proj.year <- end.year + 5

births.2014$strata <- births.2014$urban
births.2014$region <- 'All'

bb.fit.natl.yearly.u5 <- 
  SUMMER::smoothCluster(data = births.2014, 
                        family = 'betabinomial',
                        Amat = NULL, 
                        year_label = c(beg.year:end.proj.year),
                        strata.time.effect = TRUE,
                        overdisp.mean = -7.5, 
                        overdisp.prec = 0.39,
                        control.inla = list(strategy = "adaptive", int.strategy = "auto"))

bb.natl.yearly.u5 <- 
  SUMMER::getSmoothed(inla_mod = bb.fit.natl.yearly.u5, 
                      year_range = beg.year:end.proj.year, 
                      year_label = beg.year:end.proj.year,
                      nsim = 1000, 
                      weight.strata = weight.strata.natl.u5, 
                      save.draws = TRUE)

save(bb.natl.yearly.u5, file = paste0('U5MR/',country, '_res_natl_', time.model, '_yearly_u5_SmoothedCluster.rda'))
