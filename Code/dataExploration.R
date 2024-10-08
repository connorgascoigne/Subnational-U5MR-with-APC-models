rm(list=ls())

# ENTER COUNTRY OF INTEREST AND FINAL ESTIMATE INFO -----------------------------------------------
# Please capitalize the first letter of the country name and replace " " in the country name to "_" if there is.
country <- 'Kenya'

# Load libraries and info ----------------------------------------------------------

library(tidyverse)

# extract file location of this script
code.path <- rstudioapi::getActiveDocumentContext()$path
code.path.splitted <- strsplit(code.path, "/")[[1]]

# retrieve directories
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

# DHS data -----------------------------------------------

load(paste0(country,'_cluster_dat_APC.rda'), envir = .GlobalEnv)

# Load UR proportions -----------------------------------------------

if(dir.exists(paths = paste0(res.dir,'/UR/'))){
  setwd(paste0(res.dir,'/UR'))
  weight.strata.natl.u5 <- readRDS(paste0('U5_fraction/','natl_u5_urban_weights.rds'))
  weight.strata.natl.u5$rural <- 1-weight.strata.natl.u5$urban
  weight.strata.adm1.u5 <- readRDS(paste0('U5_fraction/','admin1_u5_urban_weights.rds'))
}

if(dir.exists(paths = paste0(data.dir, '/worldpop'))){
  setwd(paste0(data.dir, '/worldpop'))
  load(file = 'adm1_weights_u5.rda')
}

# Parameters ------------------------------------------------------

time.model <- c('rw2','ar1')[1]
beg.year <- 2006
end.year <- 2013 
end.proj.year <- end.year + 5
end.year.cv <- 2012 
end.proj.year.cv <- end.year.cv + 5
text.size <- 20 
height <- width <- 10

# Cluster locations -----------------------------------------------

setwd(res.dir)

adm1.cluster.location.data <-
  mod.dat %>% 
  dplyr::select(cluster, urban, LONGNUM, LATNUM, admin1.char) %>% 
  dplyr::distinct() %>% 
  dplyr::mutate(urban = urban %>% factor(.,
                                         levels = c('urban', 'rural'),
                                         labels = c('Urban', 'Rural')))
adm1.cluster.location.plot <-
  ggplot2::ggplot() +
  ggplot2::geom_sf(data = poly.adm1, fill = 'white') +
  ggplot2::geom_point(data = adm1.cluster.location.data,
                      aes(x = LONGNUM, y = LATNUM, colour = urban, shape = urban), size = 1, alpha = 0.5) +
  ggplot2::scale_shape_manual(values = c(8, 0)) +
  ggplot2::scale_colour_manual(values = c('red3', 'blue3')) + 
  mapTheme(text = element_text(size = text.size),
           legend.title = element_blank(),
           legend.position = 'top'); # adm1.cluster.location.plot

ggplot2::ggsave(plot = adm1.cluster.location.plot,
                filename = paste0('Figures/', country, '_adm1_cluster_location_plot.png'),
                height = height, width = width)

# adm1 proporions -----------------------------------------------

setwd(res.dir)

adm1.yearly.urban.proportions.data <-
  weight.strata.adm1.u5 %>% 
  tidyr::pivot_longer(cols = c('urban', 'rural'), names_to = 'strata', values_to = 'proportion') %>% 
  dplyr::mutate(strata = strata %>% factor(.,
                                         levels = c('urban', 'rural'),
                                         labels = c('Urban', 'Rural')),
                proportion = proportion * 100) %>% 
  dplyr::left_join(., admin1.names, by = c('region' = 'Internal')) %>% 
  dplyr::left_join(., poly.adm1, by = c('GADM' = 'NAME_1')) %>% 
  sf::st_as_sf()

adm1.yearly.urban.proportions.plot <-
  ggplot2::ggplot() +
  ggplot2::geom_sf(data = 
                     adm1.yearly.urban.proportions.data %>% 
                     dplyr::filter(years %in% seq(from = beg.year, to = end.proj.year, length.out = 4)), 
                   aes(fill = proportion, group = interaction(years, strata)), colour = 'black') +
  ggplot2::scale_fill_viridis_c(name = 'Strata proportions',
                                direction = -1, option = 'G', n.breaks = 5) + 
  ggplot2::facet_grid(strata ~ years) +
  mapTheme(text = element_text(size = text.size),
           legend.position = 'bottom') +
  ggplot2::guides(fill = guide_colourbar(title.position = 'top', title.hjust = 0.5)); # adm1.yearly.urban.proportions.plot

ggplot2::ggsave(plot = adm1.yearly.urban.proportions.plot,
                filename = paste0('Figures/', country, '_adm1_yearly_urban_proportions_plot.png'),
                height = height, width = width)

# adm0 proportions -----------------------------------------------

setwd(res.dir)

adm1.yearly.natl.proportions.data <-
  weight.adm1.u5 %>% 
  dplyr::mutate(proportion = proportion * 100,
                strata = 'County') %>% 
  dplyr::left_join(., admin1.names, by = c('region' = 'Internal')) %>% 
  dplyr::left_join(., poly.adm1, by = c('GADM' = 'NAME_1')) %>% 
  sf::st_as_sf()

adm1.yearly.natl.proportions.plot <-
  ggplot2::ggplot() +
  ggplot2::geom_sf(data = 
                     adm1.yearly.natl.proportions.data %>% 
                     dplyr::filter(years %in% seq(from = beg.year, to = end.proj.year, length.out = 4)), 
                   aes(fill = proportion, group = interaction(years, strata)), colour = 'black') +
  ggplot2::scale_fill_viridis_c(name = 'Strata proportions',
                                direction = -1, option = 'G', n.breaks = 5) + 
  ggplot2::facet_grid(strata ~ years) +
  mapTheme(text = element_text(size = text.size),
           legend.position = 'bottom') +
  ggplot2::guides(fill = guide_colourbar(title.position = 'top', title.hjust = 0.5)); # adm1.yearly.natl.proportions.plot

ggplot2::ggsave(plot = adm1.yearly.natl.proportions.plot,
                filename = paste0('Figures/', country, '_adm1_yearly_natl_proportions_plot.png'),
                height = height, width = width)