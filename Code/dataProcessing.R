rm(list = ls())
# ENTER COUNTRY OF INTEREST -----------------------------------------------
# Please capitalize the first letter of the country name and replace ' ' in the country name to '_' if there is.
country <- 'Kenya'

# Load libraries and info ----------------------------------------------------------
library(SUMMER)
library(stringr)
library(tidyverse)
library(rdhs)
library(sf)
library(spdep)
library(haven)

# extract file location of this script
code.path <- rstudioapi::getActiveDocumentContext()$path
code.path.splitted <- strsplit(code.path, '/')[[1]]

# retrieve directories
home.dir <- paste(code.path.splitted[1: (length(code.path.splitted)-2)], collapse = '/')
data.dir <- paste0(home.dir,'/Data/',country) # set the directory to store the data
res.dir <- paste0(home.dir,'/Results/',country) # set the directory to store the results (e.g. fitted R objects, figures, tables in .csv etc.)
info.name <- paste0(country, '_general_info.Rdata')
load(file = paste0(home.dir,'/Info/',info.name, sep='')) # load the country info

# set API to get DHS data -- you will need to change this to your information!
rdhs::set_rdhs_config(email = 'c.gascoigne@imperial.ac.uk',
                project = 'Estimating subnational under-5 mortality rates in space and time using age-period-cohort models')

rdhs::update_rdhs_config(email = 'c.gascoigne@imperial.ac.uk', password = T,
                   project = 'Estimating subnational under-5 mortality rates in space and time using age-period-cohort models')


# Load polygon files ----------------------------------------------------------

setwd(data.dir)

poly.adm0 <- sf::st_read(dsn = poly.path, layer = as.character(poly.layer.adm0)) # load the national shape file
# use encoding to read special characters
poly.adm1 <- sf::st_read(dsn = poly.path, layer = as.character(poly.layer.adm1)) # load the shape file of admin-1 regions

# Create Adjacency Matrix ----------------------------------------------------------

# Adjacency matrix is a symmetric matrix with each entry of 1's or 0's indicating if the two administrative regions are adjacent. 
# Each row or column represents an administrative region.
# The codes below generates spatial adjacency matrix based on the spatial polygon file.

# create the adjacency matrix for admin1 regions.
admin1.mat <- spdep::poly2nb(as(poly.adm1, 'Spatial'))
admin1.mat <- spdep::nb2mat(admin1.mat, zero.policy = TRUE)
colnames(admin1.mat) <- rownames(admin1.mat) <- paste0('admin1_', 1:dim(admin1.mat)[1])
admin1.names <- data.frame(GADM = poly.adm1$NAME_1,
                           Internal = rownames(admin1.mat))

save(admin1.mat, file = paste0(poly.path,'/', country, '_Amat.rda'))
save(admin1.names, file = paste0(poly.path, '/', country, '_Amat_Names.rda'))

# Find DHS surveys ----------------------------------------------------------

#get country ID
countryId <- rdhs::dhs_countries()[rdhs::dhs_countries()$ISO3_CountryCode==toupper(gadm.abbrev),]

potential_surveys <- rdhs::dhs_datasets(countryIds = countryId$DHS_CountryCode, surveyYearStart = 2000) %>% 
  dplyr::filter((FileType == 'Births Recode' & FileFormat=='Stata dataset (.dta)') |(FileType == 'Geographic Data' & FileFormat =='Flat ASCII data (.dat)'))

#only keep surveys with both a births recode and a geographic dataset
dhs_survey_ids <- as.numeric(unique(potential_surveys$SurveyNum)[sapply(unique(potential_surveys$SurveyNum),
                                                                        function(num){
                                                                          if(sum(c('Births Recode','Geographic Data') %in% (potential_surveys %>% filter(SurveyNum==num))$FileType) ==2){return(T)
                                                                          }else(return(F))})])
surveys <- 
  potential_surveys %>% 
  dplyr::filter(SurveyNum %in% dhs_survey_ids) %>% 
  dplyr::arrange(SurveyYear, DatasetType)
dhs_survey_years <- as.numeric(unique(surveys$SurveyYear))

# CHECK THAT SURVEYS FOR CORRECT COUNTRY HAVE BEEN CHOSEN
unique(surveys$CountryName)

# Process data for each DHS survey year ----------------------------------------------------------

## for APC data ----------------------------------------------------------

# The codes below first loads the raw DHS data, then it assigns the GPS coordinates to each sampling cluster and admin regions where
# the sampling is conducted and assigns the admin regions where the clusters are located.

for(survey_year in dhs_survey_years){
  
  dhs.svy.ind <- which(dhs_survey_years==survey_year)
  message('Processing DHS data for ', country,' ', survey_year,'\n')
  
  data.paths.tmp <- rdhs::get_datasets(surveys[surveys$SurveyYear==survey_year,]$FileName, clear_cache = T)
  
  raw.dat.tmp <- readRDS(paste0(data.paths.tmp[2]))
  
  # convert some variables to factors
  alive <- attr(raw.dat.tmp$b5, which = 'labels')
  names(alive) <- tolower(names(alive))
  raw.dat.tmp$b5 <- ifelse(raw.dat.tmp$b5 == alive['yes'][[1]], 'yes', 'no')
  raw.dat.tmp$b5 <- factor(raw.dat.tmp$b5, levels = c('yes', 'no'))
  
  strat <- attr(raw.dat.tmp$v025,which='labels')
  names(strat) <- tolower(names(strat))
  raw.dat.tmp$v025 <- ifelse(raw.dat.tmp$v025 == strat['urban'][[1]],'urban', 'rural')
  raw.dat.tmp$v025 <- factor(raw.dat.tmp$v025, levels = c('urban', 'rural'))
  
  if(country=='Ethiopia'){
    cmc.adjust <- 92
  }else if(country=='Nepal'){
    cmc.adjust <- -681
  }else{cmc.adjust <- 0}
  
  # read DHS data
  dat.tmp <- SUMMER::getBirths(data = raw.dat.tmp,
                               surveyyear = survey_year,
                               year.cut = seq(beg.year, survey_year + 1, 1),
                               cmc.adjust = cmc.adjust)
  
  # organise temporal
  ## columns include age, time (period), dob (cohort)
  dat.tmp <- 
    dat.tmp[,c('v001', 'v005', 'v024', 'v025', 'strata', 'age', 'time', 'dob', 'died')] %>%
    dplyr:::mutate(dob = floor((dob-1)/12)+1900) %>%
    dplyr::summarise(died = sum(died),
                     total = n(),
                     .by = c('v001', 'v005', 'v024', 'v025', 'strata', 'age', 'time', 'dob'))
  
  # specify the name of DHS GPS file, which contains the GPS coordinates of the sampling cluster where the data is sampled
  points <- readRDS(paste0(data.paths.tmp[1]))
  
  # detect points in the DHS GPS file with mis-specified coordinates and remove them if any
  wrong.points <- which(points@data$LATNUM == 0.0 & points@data$LONGNUM == 0.0)
  if(!is.null(dim(wrong.points))){message('There are wrong GPS points: (Longitude, Latitude) = (0, 0)')}
  
  # remove wrong points in the data if any
  dat.tmp <- dat.tmp[!(dat.tmp$v001 %in% points@data$DHSCLUST[wrong.points]),]
  points@data$DHSCLUST[wrong.points] %in% unique(dat.tmp$v001)
  
  # add the column for GPS coordinate in the data
  dat.tmp$LONGNUM <- dat.tmp$LATNUM <- NA
  for(i in 1:dim(points)[1]){
    dat.tmp$LATNUM[dat.tmp$v001 == points@data$DHSCLUST[i]] <- points@data$LATNUM[i] # assign latitude to DHS cluster location
    dat.tmp$LONGNUM[dat.tmp$v001 == points@data$DHSCLUST[i]] <- points@data$LONGNUM[i] # assign longitude to DHS cluster location
  }
  
  # remove missing points in the data if any
  miss <- which(dat.tmp$LATNUM == 0 & dat.tmp$LONGNUM == 0)
  if(length(miss != 0)){
    dat.tmp <- dat.tmp[-miss,]
  }
  
  message('\n Assigned LAT & LONG')
  
  
  # assign admin regions based on coordinates and polygon files
  adm1.ind <- exists('poly.adm1')
  
  points.frame <- as.data.frame(dat.tmp[,c('LONGNUM', 'LATNUM')]) # retrieve GPS coordinates where data is sampled.
  points.frame <- sp::SpatialPoints(points.frame) # convert the GPS coordinates into 'sp' object.
  
  if(adm1.ind){
    poly.over.adm1 <- sp::SpatialPolygons(poly.adm1@polygons)
    sp::proj4string(points.frame) <- sp::proj4string(poly.over.adm1) <- sp::proj4string(poly.adm1) 
    admin1.key <- sp::over(points.frame, poly.over.adm1)
    miss.frame.adm1 <- unique(points.frame@coords[which(is.na(admin1.key)),])
    
    if(dim(miss.frame.adm1)[1] != 0){
      miss.poly.adm1 <- geosphere::dist2Line( miss.frame.adm1, poly.over.adm1)
      
      for(i in 1:dim(miss.poly.adm1)[1]){
        long.ids <- which(points.frame@coords[,c('LONGNUM')] %in% miss.frame.adm1[i,1])
        lat.ids <- which(points.frame@coords[,c('LATNUM')] %in% miss.frame.adm1[i,2])
        ids <- intersect(long.ids, lat.ids)
        admin1.key[ids] <- rep(miss.poly.adm1[i, 'ID'], length(ids))
      }
    }
    
    dat.tmp$admin1 <- admin1.key
    dat.tmp$admin1.char <- paste0('admin1_', admin1.key)
    dat.tmp$admin1.name <- as.character(eval(str2lang(poly.label.adm1)))[admin1.key]
  }else{
    dat.tmp$admin1 <- dat.tmp$admin1.name <- NA
    message('There is no Admin1 polygon to assign points to.')
  }  
  
  if(FALSE){
    check <- dat.tmp$strata
    check <- gsub(' - rural', '', check)
    check <- gsub(' - urban', '', check)
    inconsist <- which(check != tolower(dat.tmp$admin1.name))
    table(check[inconsist], dat.tmp$admin1.name[inconsist])
    unique(dat.tmp[which(check == 'neno' & dat.tmp$admin1.name == 'Balaka'), 'v001'])
  }
  
  # finish preparing data ###
  dat.tmp <- dat.tmp[,c('v001', 'age', 'time', 'dob', 'total', 'died', 'v005', 'v025', 'LONGNUM', 'LATNUM','strata',
                        'admin1', 'admin1.char', 'admin1.name')]
  colnames(dat.tmp) <- c('cluster', 'age', 'years', 'cohort', 'total', 'Y', 'v005', 'urban', 'LONGNUM', 'LATNUM','strata',
                         'admin1', 'admin1.char', 'admin1.name')
  
  dat.tmp$survey <- raw.dat.tmp$survey_year <- survey_year
  dat.tmp$survey.type <- 'DHS'
  
  if(survey_year==dhs_survey_years[1]){
    mod.dat <- dat.tmp
    raw.dat <- raw.dat.tmp[,c('caseid', 'b5', 'b7', 'survey_year')]
  }else{mod.dat <- rbind(mod.dat,dat.tmp)
  raw.dat <- rbind(raw.dat,raw.dat.tmp[,c('caseid', 'b5', 'b7', 'survey_year')])}
  
}

### Change cluster numbers to get rid of duplicates ----------------------------------------------------------

clusters <- unique(mod.dat[,c('cluster', 'survey')])
clusters$cluster.new <- 1:nrow(clusters)
mod.dat <- merge(mod.dat,clusters,by=c('cluster', 'survey'))
mod.dat$cluster <- mod.dat$cluster.new
mod.dat <- mod.dat[,!(names(mod.dat)=='cluster.new')]

survey_years <- sort(unique(mod.dat$survey))
mod.dat$survey.id<- unlist(sapply(1:nrow(mod.dat),function(x){which(mod.dat$survey[x] ==survey_years)}))


### Save processed data  ----------------------------------------------------------

save(mod.dat, file = paste0(country,'_cluster_dat_APC.rda'))

## for SUMMER data ----------------------------------------------------------

# The codes below first loads the raw DHS data, then it assigns the GPS coordinates to each sampling cluster and admin regions where
# the sampling is conducted and assigns the admin regions where the clusters are located.

for(survey_year in dhs_survey_years){
  
  dhs.svy.ind <- which(dhs_survey_years==survey_year)
  message('Processing DHS data for ', country,' ', survey_year,'\n')
  
  data.paths.tmp <- get_datasets(surveys[surveys$SurveyYear==survey_year,]$FileName, clear_cache = T)
  
  raw.dat.tmp <- readRDS(paste0(data.paths.tmp[2]))
  
  # convert some variables to factors
  alive <- attr(raw.dat.tmp$b5, which = 'labels')
  names(alive) <- tolower(names(alive))
  raw.dat.tmp$b5 <- ifelse(raw.dat.tmp$b5 == alive['yes'][[1]], 'yes', 'no')
  raw.dat.tmp$b5 <- factor(raw.dat.tmp$b5, levels = c('yes', 'no'))
  
  strat <- attr(raw.dat.tmp$v025,which='labels')
  names(strat) <- tolower(names(strat))
  raw.dat.tmp$v025 <- ifelse(raw.dat.tmp$v025 == strat['urban'][[1]],'urban','rural')
  raw.dat.tmp$v025 <- factor(raw.dat.tmp$v025, levels = c('urban','rural'))
  
  if(country=='Ethiopia'){
    cmc.adjust <- 92
  }else if(country=='Nepal'){
    cmc.adjust <- -681
  }else{cmc.adjust <- 0}
  
  # read DHS data
  dat.tmp <- getBirths(data=raw.dat.tmp,
                       surveyyear = survey_year,
                       year.cut = seq(beg.year, survey_year + 1, 1),
                       cmc.adjust = cmc.adjust,compact = T)
  
  
  # retrieve the some columns of the full data
  dat.tmp <- dat.tmp[ ,c('v001', 'v024', 'time', 'total',
                         'age', 'v005', 'v025', 'strata', 'died')]
  
  # specify the name of DHS GPS file, which contains the GPS coordinates of the sampling cluster where the data is sampled
  points <- readRDS(paste0(data.paths.tmp[1]))
  
  # detect points in the DHS GPS file with mis-specified coordinates and remove them if any
  wrong.points <- which(points@data$LATNUM == 0.0 & points@data$LONGNUM == 0.0)
  if(!is.null(dim(wrong.points))){message('There are wrong GPS points: (Longitude, Latitude) = (0, 0)')}
  
  # remove wrong points in the data if any
  dat.tmp <- dat.tmp[!(dat.tmp$v001 %in% points@data$DHSCLUST[wrong.points]),]
  points@data$DHSCLUST[wrong.points] %in% unique(dat.tmp$v001)
  
  # add the column for GPS coordinate in the data
  dat.tmp$LONGNUM <- dat.tmp$LATNUM <- NA
  for(i in 1:dim(points)[1]){
    dat.tmp$LATNUM[dat.tmp$v001 == points@data$DHSCLUST[i]] <- points@data$LATNUM[i] # assign latitude to DHS cluster location
    dat.tmp$LONGNUM[dat.tmp$v001 == points@data$DHSCLUST[i]] <- points@data$LONGNUM[i] # assign longitude to DHS cluster location
  }
  
  # remove missing points in the data if any
  miss <- which(dat.tmp$LATNUM == 0 & dat.tmp$LONGNUM == 0)
  if(length(miss != 0)){
    dat.tmp <- dat.tmp[-miss,]
  }
  
  message('\n Assigned LAT & LONG')
  
  # assign admin regions based on coordinates and polygon files
  adm1.ind <- exists('poly.adm1')
  
  points.frame <- as.data.frame(dat.tmp[,c('LONGNUM', 'LATNUM')]) # retrieve GPS coordinates where data is sampled.
  points.frame <- SpatialPoints(points.frame) # convert the GPS coordinates into 'sp' object.
  
  if(adm1.ind){
    poly.over.adm1 <- SpatialPolygons(poly.adm1@polygons)
    proj4string(points.frame) <- proj4string(poly.over.adm1) <- 
      proj4string(poly.adm1) 
    admin1.key <- over(points.frame, poly.over.adm1)
    miss.frame.adm1 <- unique(points.frame@coords[which(is.na(admin1.key)),])
    
    if(dim(miss.frame.adm1)[1] != 0){
      miss.poly.adm1 <- dist2Line( miss.frame.adm1, poly.over.adm1)
      
      for(i in 1:dim(miss.poly.adm1)[1]){
        long.ids <- which(points.frame@coords[,c('LONGNUM')] %in% miss.frame.adm1[i,1])
        lat.ids <- which(points.frame@coords[,c('LATNUM')] %in% miss.frame.adm1[i,2])
        ids <- intersect(long.ids, lat.ids)
        admin1.key[ids] <- rep(miss.poly.adm1[i, 'ID'], length(ids))
      }
    }
    
    dat.tmp$admin1 <- admin1.key
    dat.tmp$admin1.char <- paste0('admin1_', admin1.key)
    dat.tmp$admin1.name <- as.character(eval(str2lang(poly.label.adm1)))[admin1.key]
  }else{
    dat.tmp$admin1 <- dat.tmp$admin1.name <- NA
    message('There is no Admin1 polygon to assign points to.')
  }  
  
  if(FALSE){
    check <- dat.tmp$strata
    check <- gsub(' - rural', '', check)
    check <- gsub(' - urban', '', check)
    inconsist <- which(check != tolower(dat.tmp$admin1.name))
    table(check[inconsist], dat.tmp$admin1.name[inconsist])
    unique(dat.tmp[which(check == 'neno' & dat.tmp$admin1.name == 'Balaka'), 'v001'])
  }
  
  # finish preparing data ###
  dat.tmp <- dat.tmp[,c('v001', 'age', 'time', 'total', 'died', 'v005', 'v025', 'LONGNUM', 'LATNUM','strata',
                        'admin1', 'admin1.char', 'admin1.name')]
  colnames(dat.tmp) <- c('cluster', 'age', 'years', 'total', 'Y', 'v005', 'urban', 'LONGNUM', 'LATNUM','strata',
                         'admin1', 'admin1.char', 'admin1.name')
  
  dat.tmp$survey <- raw.dat.tmp$survey_year <-survey_year
  dat.tmp$survey.type <- 'DHS'
  
  if(survey_year==dhs_survey_years[1]){
    mod.dat <- dat.tmp
    raw.dat <- raw.dat.tmp[,c('caseid', 'b5','b7','survey_year')]
  }else{mod.dat <- rbind(mod.dat,dat.tmp)
  raw.dat <- rbind(raw.dat,raw.dat.tmp[,c('caseid', 'b5','b7','survey_year')])}
  
}

### Change cluster numbers to get rid of duplicates ----------------------------------------------------------

clusters <- unique(mod.dat[,c('cluster', 'survey')])
clusters$cluster.new <- 1:nrow(clusters)
mod.dat <- merge(mod.dat,clusters,by=c('cluster', 'survey'))
mod.dat$cluster <- mod.dat$cluster.new
mod.dat <- mod.dat[,!(names(mod.dat)=='cluster.new')]

survey_years <- sort(unique(mod.dat$survey))
mod.dat$survey.id<- unlist(sapply(1:nrow(mod.dat),function(x){which(mod.dat$survey[x] ==survey_years)}))


### Save processed data  ----------------------------------------------------------

save(mod.dat, file = paste0(country,'_cluster_dat_SUMMER.rda'))