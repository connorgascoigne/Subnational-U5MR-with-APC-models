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

# Load polygon files ------------------------------------------------------

setwd(data.dir)

poly.adm0 <- sf::st_read(dsn = poly.path, layer = as.character(poly.layer.adm0)) # load the national shape file
poly.adm1 <- sf::st_read(dsn = poly.path, layer = as.character(poly.layer.adm1)) # load the shape file of admin-1 regions

load(paste0(poly.path,'/', country, '_Amat.rda'))
load(paste0(poly.path,'/', country, '_Amat_Names.rda'))

# Parameters ------------------------------------------------------

time.model <- c('rw2','ar1')[1]
beg.year <- 2006
end.year <- 2013 
end.proj.year <- end.year + 5
end.year.cv <- 2012 
text.size <- 20 
height <- width <- 10

# Load model results ------------------------------------------------------

setwd(res.dir)

# direct estimates
load(file = paste0('Direct/U5MR/',country, '_direct_natl_yearly_u5.rda'))
load(file = paste0('Direct/U5MR/',country, '_direct_adm1_yearly_u5.rda'))
# FH estimates
load(file = paste0('Direct/U5MR/',country, '_res_natl_', time.model, '_yearly_u5_SmoothedDirect.rda'))
load(file = paste0('Direct/U5MR/',country, '_res_adm1_', time.model, '_yearly_u5_SmoothedDirect.rda'))
# apc
load(file = paste0('Age-Period-Cohort/U5MR/',country,'_res_adm1_u5.rda'))
load(file = paste0('Age-Period-Cohort/U5MR/',country,'_res_adm1_u5_natl.rda'))
load(file = paste0('Age-Period-Cohort/U5MR/',country,'_fit_adm1_u5.rda'))
load(file = paste0('Age-Period-Cohort/U5MR/',country,'_cv_adm1_u5.rda'))
# ap
load(file = paste0('Age-Period/U5MR/',country,'_res_adm1_u5.rda'))
load(file = paste0('Age-Period/U5MR/',country,'_res_adm1_u5_natl.rda'))
load(file = paste0('Age-Period/U5MR/',country,'_fit_adm1_u5.rda'))
load(file = paste0('Age-Period/U5MR/',country,'_cv_adm1_u5.rda'))
# ac
load(file = paste0('Age-Cohort/U5MR/',country,'_res_adm1_u5.rda'))
load(file = paste0('Age-Cohort/U5MR/',country,'_res_adm1_u5_natl.rda'))
load(file = paste0('Age-Cohort/U5MR/',country,'_fit_adm1_u5.rda'))
load(file = paste0('Age-Cohort/U5MR/',country,'_cv_adm1_u5.rda'))

# Estimate plots ------------------------------------------------------

## national ------------------------------------------------------

### all line plots ------------------------------------------------------

setwd(res.dir)

natl.frame <-
  rbind(direct.natl.yearly.u5 %>%
          dplyr::select(region, years, lower, mean, upper) %>% 
          dplyr::filter(years %in% beg.year:end.year) %>% 
          dplyr::rename(median = mean) %>% # need to do this for naming as cannot get median from DE
          dplyr::mutate(type = 'Direct'),
        res.natl.yearly.u5 %>% 
          dplyr::select(region, years, lower, median, upper) %>% 
          dplyr::mutate(type = 'Fay-Herriot'),
        apc.adm1.u5.natl.res$national %>% 
          dplyr::select(region, years, lower, median, upper) %>% 
          dplyr::mutate(type = 'Age-Period-Cohort'),
        ap.adm1.u5.natl.res$national %>% 
          dplyr::select(region, years, lower, median, upper) %>% 
          dplyr::mutate(type = 'Age-Period'),
        ac.adm1.u5.natl.res$national %>% 
          dplyr::select(region, years, lower, median, upper) %>% 
          dplyr::mutate(type = 'Age-Cohort')) %>% 
  dplyr::mutate(years = years %>% as.numeric(),
                lower = 1000 * lower,
                median = 1000 * median,
                upper = 1000 * upper,
                type = type %>% factor(., levels = c('Direct', 'Fay-Herriot', 'Age-Period', 'Age-Cohort', 'Age-Period-Cohort'))) %>% 
  dplyr::filter(years %in% beg.year:end.proj.year)


natl.line.plot <- 
  ggplot2::ggplot(data = natl.frame, 
                  aes(x = years, y = median, colour = type, group = type, fill = type)) +
  ggplot2::geom_vline(aes(xintercept = end.year), linetype = 'dashed', colour = 'black') +
  ggplot2::geom_hline(aes(yintercept = 25), colour = 'black', linetype = 'dashed') +
  ggplot2::annotate('text', x = 2008, y = 26, label = 'SDG 3.2: 25 deaths',
                    size = 7.5) +
  ggplot2::geom_line() +
  ggplot2::geom_point() +
  ggplot2::scale_x_continuous(breaks = seq(from = beg.year, to = end.proj.year, by = 1)) +
  ggplot2::scale_colour_manual(values = c('red3', 'blue3', 'green4', 'orange2', 'purple3')) + 
  ggplot2::labs(x = 'Period', y = 'U5MR (deaths per 1000 live births)') +
  plotTheme(text = element_text(size = text.size),
            axis.text.x = element_text(angle = 90, hjust = 1),
            legend.title = element_blank(),
            legend.position = 'top'); natl.line.plot

natl.error.plot <-
  ggplot2::ggplot(data = natl.frame, 
                  aes(x = years, y = median, colour = type, group = type, fill = type)) +
  ggplot2::geom_vline(aes(xintercept = end.year + 0.5), linetype = 'dashed', colour = 'black') +
  ggplot2::geom_hline(aes(yintercept = 25), colour = 'black', linetype = 'dashed') +
  ggplot2::annotate('text', x = 2008, y = 28, label = 'SDG 3.2: 25 deaths',
                    size = 7.5) +
  ggplot2::geom_point(position = position_dodge(0.5)) +
  ggplot2::geom_errorbar(aes(ymin = lower, ymax = upper), width = .2, position = position_dodge(0.5)) +
  ggplot2::scale_x_continuous(breaks = seq(from = beg.year, to = end.proj.year, by = 1)) +
  ggplot2::scale_colour_manual(values = c('red3', 'blue3', 'green4', 'orange2', 'purple3')) + 
  ggplot2::labs(x = 'Period', y = 'U5MR (deaths per 1000 live births)') +
  plotTheme(text = element_text(size = text.size),
            axis.text.x = element_text(angle = 90, hjust = 1),
            legend.title = element_blank(),
            legend.position = 'top'); natl.error.plot

ggplot2::ggsave(plot = natl.line.plot,
                filename = paste0('Figures/', country, '_natl_line_plot.png'),
                height = height, width = width)

ggplot2::ggsave(plot = natl.error.plot,
                filename = paste0('Figures/', country, '_natl_error_plot.png'),
                height = height, width = width)

### apc vs other ------------------------------------------------------

setwd(res.dir)

natl.frame.wide <- 
  natl.frame %>% 
  dplyr::select(-lower, -upper) %>% 
  tidyr::pivot_wider(names_from = type, values_from = median) %>% 
  dplyr::filter(years %in% 2006:2013)

apc.direct.natl.u5.compairson.plot <- 
  ggplot2::ggplot(natl.frame.wide, aes(x = `Age-Period-Cohort`, y = Direct, label = years)) +
  ggplot2::geom_abline(intercept = 0) +
  ggplot2::coord_fixed(xlim = c(40, 60), ylim = c(40, 60)) +
  ggplot2::geom_point() +
  ggrepel::geom_text_repel() +
  ggplot2::labs(x = 'APC U5MR (per 1000 live births)', y = 'Direct Estimates U5MR (per 1000 live births)') +
  plotTheme(text = element_text(size = text.size)); apc.direct.natl.u5.compairson.plot

apc.ap.natl.u5.compairson.plot <- 
  ggplot2::ggplot(natl.frame.wide, aes(x = `Age-Period-Cohort`, y = `Age-Period`, label = years)) +
  ggplot2::geom_abline(intercept = 0) +
  ggplot2::coord_fixed(xlim = c(40, 60), ylim = c(40, 60)) +
  ggplot2::geom_point() +
  ggrepel::geom_text_repel() +
  ggplot2::labs(x = 'APC U5MR (per 1000 live births)', y = 'AP U5MR (per 1000 live births)') +
  plotTheme(text = element_text(size = text.size)); apc.ap.natl.u5.compairson.plot

apc.ac.natl.u5.compairson.plot <- 
  ggplot2::ggplot(natl.frame.wide, aes(x = `Age-Period-Cohort`, y = `Age-Cohort`, label = years)) +
  ggplot2::geom_abline(intercept = 0) +
  ggplot2::coord_fixed(xlim = c(40, 60), ylim = c(40, 60)) +
  ggplot2::geom_point() +
  ggrepel::geom_text_repel() +
  ggplot2::labs(x = 'APC U5MR (per 1000 live births)', y = 'AC U5MR (per 1000 live births)') +
  plotTheme(text = element_text(size = text.size)); apc.ac.natl.u5.compairson.plot

ggplot2::ggsave(plot = apc.direct.natl.u5.compairson.plot,
                filename = paste0('Figures/', country, '_apc_direct_natl_u5_compairson_plot.png'),
                height = height, width = width)

ggplot2::ggsave(plot = apc.ap.natl.u5.compairson.plot,
                filename = paste0('Figures/', country, '_apc_ap_natl_u5_compairson_plot.png'),
                height = height, width = width)

ggplot2::ggsave(plot = apc.ac.natl.u5.compairson.plot,
                filename = paste0('Figures/', country, '_apc_ac_natl_u5_compairson_plot.png'),
                height = height, width = width)

## admin1 ------------------------------------------------------

setwd(res.dir)

adm1.frame <-
  rbind(direct.admin1.yearly.u5 %>%
          dplyr::select(region, years, lower, mean, upper) %>% 
          dplyr::filter(region != 'All') %>% 
          dplyr::rename(median = mean) %>% # need to do this for naming as cannot get median from DE
          dplyr::mutate(type = 'Direct'),
        sd.admin1.yearly.u5 %>% 
          dplyr::select(region, years, lower, median, upper) %>% 
          dplyr::mutate(type = 'Fay-Herriot'),
        apc.adm1.u5.res$subnational %>% 
          dplyr::select(region, years, lower, median, upper) %>% 
          dplyr::mutate(type = 'Age-Period-Cohort'),
        ap.adm1.u5.res$subnational %>% 
          dplyr::select(region, years, lower, median, upper) %>% 
          dplyr::mutate(type = 'Age-Period'),
        ac.adm1.u5.res$subnational %>% 
          dplyr::select(region, years, lower, median, upper) %>% 
          dplyr::mutate(type = 'Age-Cohort')) %>% 
  dplyr::mutate(years = years %>% as.numeric(),
                lower = 1000 * lower,
                median = 1000 * median,
                upper = 1000 * upper,
                type = type %>% factor(., levels = c('Direct', 'Fay-Herriot', 'Age-Period', 'Age-Cohort', 'Age-Period-Cohort'))) %>% 
  dplyr::filter(years %in% beg.year:end.proj.year)

### model parameters ------------------------------------------------------

setwd(res.dir)

ap.adm1.u5.parameters <-
  dplyr::bind_rows(ap.adm1.u5.fit$fit$summary.fixed[,c(3:5)], 
                   ap.adm1.u5.fit$fit$summary.hyperpar[,c(3:5)])

ac.adm1.u5.parameters <-
  dplyr::bind_rows(ac.adm1.u5.fit$fit$summary.fixed[,c(3:5)], 
                   ac.adm1.u5.fit$fit$summary.hyperpar[,c(3:5)])

apc.adm1.u5.parameters <-
  dplyr::bind_rows(apc.adm1.u5.fit$fit$summary.fixed[,c(3:5)], 
                   apc.adm1.u5.fit$fit$summary.hyperpar[,c(3:5)])

all.adm1.u5.parameters <-
  data.frame(parameter = 
               c('age_id', 'period_id', 'cohort_id', 
                 'strata_idrural', 'strata_idurban',
                 'overdispersion for the betabinomial observations',
                 'Precision for age2_id', 'Precision for period2_id', 'Precision for cohort2_id', 
                 'Precision for space_id', 'Phi for space_id', 
                 'Precision for spaceTime_id'),
             parameterLabel = 
               c('Age slope', 'Period slope', 'Cohort slope', 
                 'Rural intercept', 'Urban intercept', 
                 'Betabinomial overdispersion',
                 'Age precision', 'Period precision', 'Cohort precision', 
                 'Region precision', 'Region mixing', 
                 'Space-period precision')) %>% 
  dplyr::left_join(., 
                   ap.adm1.u5.parameters %>% tibble::rownames_to_column(., var = 'parameter'),
                   by = 'parameter') %>% 
  dplyr::left_join(., 
                   ac.adm1.u5.parameters %>% tibble::rownames_to_column(., var = 'parameter'),
                   by = 'parameter') %>% 
  dplyr::left_join(., 
                   apc.adm1.u5.parameters %>% tibble::rownames_to_column(., var = 'parameter'),
                   by = 'parameter')

print(xtable::xtable(all.adm1.u5.parameters %>% dplyr::select(-parameter), digits = 2), include.rownames=FALSE)

options(scipen = 2)
print(xtable::xtable((all.adm1.u5.parameters %>%
                        dplyr::filter(parameterLabel == 'Betabinomial overdispersion') %>%
                        dplyr::select(-parameter, -parameterLabel))*100000 %>% 
                       round(., digits = 2)))


print(xtable::xtable(ap.adm1.u5.parameters, digits = 2))
print(xtable::xtable(ac.adm1.u5.parameters, digits = 2))
print(xtable::xtable(apc.adm1.u5.parameters, digits = 2))

### apc subnational maps ------------------------------------------------------

setwd(res.dir)

apc.adm1.u5.plot.data <-
  apc.adm1.u5.res$subnational %>% 
  dplyr::left_join(., admin1.names, by = c('region' = 'Internal')) %>% 
  dplyr::left_join(., poly.adm1, by = c('GADM' = 'NAME_1')) %>% 
  dplyr::mutate(lower = lower * 1000,
                median = median * 1000,
                upper = upper * 1000,
                width = upper - lower) %>% 
  sf::st_as_sf()

apc.adm1.u5.estimate.plot <-
  ggplot2::ggplot() +
  ggplot2::geom_sf(data = apc.adm1.u5.plot.data, aes(fill = median), colour = 'black') +
  ggplot2::scale_fill_viridis_c(name = 'U5MR (death per 1000 live births)',
                                direction = -1, option = 'G', n.breaks = 5) + 
  ggplot2::facet_wrap(~ years, nrow = 3) +
  mapTheme(text = element_text(size = text.size),
           legend.position = 'bottom') +
  ggplot2::guides(fill = guide_colourbar(title.position = 'top', title.hjust = 0.5)); #apc.adm1.u5.estimate.plot

apc.adm1.u5.width.plot <-
  ggplot2::ggplot() +
  ggplot2::geom_sf(data = apc.adm1.u5.plot.data, aes(fill = width), colour = 'black') +
  ggplot2::scale_fill_viridis_c(name = 'Width of Credible Interval',
                                direction = -1, option = 'A', n.breaks = 3) + 
  ggplot2::facet_wrap(~ years, nrow = 3) +
  mapTheme(text = element_text(size = text.size),
           legend.position = 'bottom') +
  ggplot2::guides(fill = guide_colourbar(title.position = 'top', title.hjust = 0.5)); #apc.adm1.u5.width.plot

ggplot2::ggsave(plot = apc.adm1.u5.estimate.plot,
                filename = paste0('Figures/', country, '_apc_adm1_u5_estimate_plot.png'),
                height = height, width = width)

ggplot2::ggsave(plot = apc.adm1.u5.width.plot,
                filename = paste0('Figures/', country, '_apc_adm1_u5_width_plot.png'),
                height = height, width = width)

### apc line plot ------------------------------------------------------

setwd(res.dir)

apc.adm1.u5.plot.data <-
  apc.adm1.u5.res$subnational %>% 
  dplyr::left_join(., admin1.names, by = c('region' = 'Internal')) %>% 
  dplyr::left_join(., poly.adm1, by = c('GADM' = 'NAME_1')) %>% 
  dplyr::mutate(lower = lower * 1000,
                median = median * 1000,
                upper = upper * 1000,
                width = upper - lower)

apc.adm1.u5.lineplot <- 
  ggplot2::ggplot(data = apc.adm1.u5.plot.data, aes(x = years, y = median, group = GADM, colour = GADM)) +
  ggplot2::geom_vline(aes(xintercept = end.year), linetype = 'dashed') +
  ggplot2::geom_line() +
  ggplot2::scale_x_continuous(breaks = seq(from = beg.year, to = end.proj.year, by = 1)) +
  ggplot2::labs(x = 'Period', y = 'U5MR (deaths per 1000 live births)') +
  plotTheme(text = element_text(size = text.size),
            axis.text.x = element_text(angle = 90, hjust = 1),
            legend.position = 'none',
            legend.title = element_blank()); apc.adm1.u5.lineplot

ggplot2::ggsave(plot = apc.adm1.u5.lineplot,
                filename = paste0('Figures/', country, '_apc_adm1_u5_lineplot.png'),
                height = height, width = width)

### apc age by period and strata ------------------------------------------------------

setwd(res.dir)

apc.adm1.u5.agePeriodStrata.data <-
  apc.adm1.u5.res$linpred.draws %>% 
  dplyr::summarise(dplyr::across(dplyr::starts_with('theta:'), mean),
                   .by = c('strata', 'age', 'years')) %>%
  dplyr::mutate(dplyr::select(., starts_with('theta:')) %>% 
                  # SUMMER::expit() %>% 
                  apply(., 1, quantileSummary) %>% 
                  lapply(., data.frame) %>%
                  do.call(rbind, .),
                Age = dplyr::case_when(age == 0 ~ '[0, 1)',
                                       age == 6 ~ '[1, 12)',
                                       age == 17.5 ~ '[12, 24)',
                                       age == 29.5 ~ '[24, 36)',
                                       age == 41.5 ~ '[36, 48)',
                                       TRUE ~ '[48, 60)')) %>%
  dplyr::select(-starts_with('theta:')) %>% 
  dplyr::mutate(strata = 
                  strata %>% 
                  factor(., 
                         levels = c('urban', 'rural'),
                         labels = c('Urban', 'Rural')))

apc.adm1.u5.agePeriodStrata.lineplot <- 
  ggplot2::ggplot(data = apc.adm1.u5.agePeriodStrata.data, aes(x = years, y = median, group = interaction(Age, strata), colour = Age, shape = Age, linetype = Age)) +
  ggplot2::geom_vline(xintercept = end.year + 0.5, linetype = 'dashed') + 
  ggplot2::geom_point(position = position_dodge(1)) +
  ggplot2::geom_errorbar(aes(ymin = lower, ymax = upper), width = .2, position = position_dodge(1)) +
  ggplot2::facet_grid(~ strata) +
  ggplot2::scale_x_continuous(breaks = seq(from = min(apc.adm1.u5.agePeriodStrata.data$years), to = max(apc.adm1.u5.agePeriodStrata.data$years), by = 1)) +
  ggplot2::scale_colour_manual(values = c('red3', 'blue3', 'green4', 'orange2', 'purple3', 'pink2')) + 
  ggplot2::labs(x = 'Period', y = 'Logit Mortality') +
  plotTheme(text = element_text(size = text.size),
            axis.text.x = element_text(angle = 90, hjust = 1),
            legend.position = 'top') +
  ggplot2::guides(colour = guide_legend(nrow = 1)); apc.adm1.u5.agePeriodStrata.lineplot

ggplot2::ggsave(plot = apc.adm1.u5.agePeriodStrata.lineplot,
                filename = paste0('Figures/', country, '_apc_adm1_u5_agePeriodStrata_lineplot.png'),
                height = height, width = width)


### apc age by cohort and strata ------------------------------------------------------

setwd(res.dir)

apc.adm1.u5.ageCohortStrata.data <-
  apc.adm1.u5.res$linpred.draws %>% 
  dplyr::summarise(dplyr::across(dplyr::starts_with('theta:'), mean),
                   .by = c('strata', 'age', 'cohort')) %>%
  dplyr::mutate(dplyr::select(., starts_with('theta:')) %>% 
                  # SUMMER::expit() %>% 
                  apply(., 1, quantileSummary) %>% 
                  lapply(., data.frame) %>%
                  do.call(rbind, .),
                Age = dplyr::case_when(age == 0 ~ '[0, 1)',
                                       age == 6 ~ '[1, 12)',
                                       age == 17.5 ~ '[12, 24)',
                                       age == 29.5 ~ '[24, 36)',
                                       age == 41.5 ~ '[36, 48)',
                                       TRUE ~ '[48, 60)')) %>%
  dplyr::group_by(Age) %>% 
  dplyr::filter(dplyr::case_when(Age != '[0, 1)' ~ cohort < max(cohort),
                                 TRUE ~ cohort %in% min(cohort):max(cohort))) %>% 
  dplyr::ungroup() %>% 
  dplyr::select(-starts_with('theta:')) %>% 
  dplyr::mutate(strata = 
                  strata %>% 
                  factor(., 
                         levels = c('urban', 'rural'),
                         labels = c('Urban', 'Rural')))

apc.adm1.u5.ageCohortStrata.lineplot <- 
  ggplot2::ggplot(data = apc.adm1.u5.ageCohortStrata.data, aes(x = cohort, y = median, group = interaction(Age, strata), colour = Age, shape = Age, linetype = Age)) +
  # ggplot2::geom_vline(xintercept = end.year + 0.5, linetype = 'dashed') + 
  ggplot2::geom_point(position = position_dodge(1)) +
  ggplot2::geom_errorbar(aes(ymin = lower, ymax = upper), width = .2, position = position_dodge(1)) +
  ggplot2::facet_grid(~ strata) +
  ggplot2::scale_x_continuous(breaks = seq(from = min(apc.adm1.u5.ageCohortStrata.data$cohort), to = max(apc.adm1.u5.ageCohortStrata.data$cohort), by = 1)) +
  ggplot2::scale_colour_manual(values = c('red3', 'blue3', 'green4', 'orange2', 'purple3', 'pink2')) + 
  ggplot2::labs(x = 'Cohort', y = 'Logit Mortality') +
  plotTheme(text = element_text(size = text.size),
            axis.text.x = element_text(angle = 90, hjust = 1),
            # legend.title = element_blank(),
            legend.position = 'top') +
  ggplot2::guides(colour = guide_legend(nrow = 1)); apc.adm1.u5.ageCohortStrata.lineplot

ggplot2::ggsave(plot = apc.adm1.u5.ageCohortStrata.lineplot,
                filename = paste0('Figures/', country, '_apc_adm1_u5_ageCohortStrata_lineplot.png'),
                height = height, width = width)

### apc vs other line plots ------------------------------------------------------

setwd(res.dir)

adm1.frame.wide <- 
  adm1.frame %>% 
  dplyr::select(-lower, -upper) %>% 
  tidyr::pivot_wider(names_from = type, values_from = median) %>% 
  dplyr::filter(years %in% 2006:2013)

apc.direct.adm1.u5.compairson.plot <- 
  ggplot2::ggplot(adm1.frame.wide, aes(x = `Age-Period-Cohort`, y = Direct)) +
  ggplot2::geom_abline(intercept = 0) +
  ggplot2::coord_fixed(xlim = c(0, 210), ylim = c(0, 210)) +
  ggplot2::geom_point() +
  ggplot2::labs(x = 'APC U5MR (per 1000 live births)', y = 'Direct Estimates U5MR (per 1000 live births)') +
  plotTheme(text = element_text(size = text.size)); apc.direct.adm1.u5.compairson.plot

apc.ap.adm1.u5.compairson.plot <- 
  ggplot2::ggplot(adm1.frame.wide, aes(x = `Age-Period-Cohort`, y = `Age-Period`)) +
  ggplot2::geom_abline(intercept = 0) +
  ggplot2::coord_fixed(xlim = c(0, 210), ylim = c(0, 210)) +
  ggplot2::geom_point() +
  ggplot2::labs(x = 'APC U5MR (per 1000 live births)', y = 'AP U5MR (per 1000 live births)') +
  plotTheme(text = element_text(size = text.size)); apc.ap.adm1.u5.compairson.plot

apc.ac.adm1.u5.compairson.plot <- 
  ggplot2::ggplot(adm1.frame.wide, aes(x = `Age-Period-Cohort`, y = `Age-Cohort`)) +
  ggplot2::geom_abline(intercept = 0) +
  ggplot2::coord_fixed(xlim = c(0, 210), ylim = c(0, 210)) +
  ggplot2::geom_point() +
  ggplot2::labs(x = 'APC U5MR (per 1000 live births)', y = 'AC U5MR (per 1000 live births)') +
  plotTheme(text = element_text(size = text.size)); apc.ac.adm1.u5.compairson.plot

ggplot2::ggsave(plot = apc.direct.adm1.u5.compairson.plot,
                filename = paste0('Figures/', country, '_apc_direct_adm1_u5_compairson_plot.png'),
                height = height, width = width)

ggplot2::ggsave(plot = apc.ap.adm1.u5.compairson.plot,
                filename = paste0('Figures/', country, '_apc_ap_adm1_u5_compairson_plot.png'),
                height = height, width = width)

ggplot2::ggsave(plot = apc.ac.adm1.u5.compairson.plot,
                filename = paste0('Figures/', country, '_apc_ac_adm1_u5_compairson_plot.png'),
                height = height, width = width)


# cross validation plots ------------------------------------------------------

# these already have the direct estimates attached !
adm1.cv <-
  dplyr::bind_rows(
    apc.adm1.u5.cv %>% dplyr::mutate(type = 'Age-Period-Cohort'),
    ap.adm1.u5.cv %>% dplyr::mutate(type = 'Age-Period'),
    ac.adm1.u5.cv %>% dplyr::mutate(type = 'Age-Cohort')) %>% 
  dplyr::mutate(type = type %>% factor(., levels = c('Age-Period', 'Age-Cohort', 'Age-Period-Cohort')))

## scores table -----------------------------------------------

setwd(res.dir)

adm1.cv.point.scores <-
  adm1.cv %>% 
  dplyr::mutate(mae.dist = 
                  dplyr::across(dplyr::starts_with('theta:'), ~ abs(.x - logit.est)) %>% 
                  rowMeans(),
                mse.dist = 
                  dplyr::across(dplyr::starts_with('theta:'), ~ (.x - logit.est)^2) %>% 
                  rowMeans()) %>%
  dplyr::select(-starts_with('theta:'), 
                -logit.est, -logit.prec, -mean, -variance, -median, -upper, -lower) %>% 
  tidyr::pivot_longer(cols = c('mae.dist', 'mse.dist'),
                      names_to = 'metric', values_to = 'value') %>% 
  dplyr::mutate(metric = 
                  metric %>% 
                  factor(., 
                         levels = c('mae.dist', 'mse.dist'),
                         labels = c('MAE', 'MSE'))) %>%
  dplyr::summarise(value = value %>% mean(., na.rm = TRUE),
                   .by = c('type', 'metric')) %>%
  tidyr::pivot_wider(names_from = metric, values_from = 'value') %>% 
  dplyr::mutate(MAE = MAE*100,
                MSE = MSE*100); adm1.cv.point.scores

adm1.cv.dist.scores <- 
  adm1.cv %>% 
  dplyr::mutate(dplyr::select(., dplyr::starts_with('theta:')) %>% 
                  apply(., 1, quantileSummary, CI = 0.50) %>% 
                  lapply(., data.frame) %>% 
                  do.call(rbind, .),
                lower50 = lower,
                upper50 = upper) %>% 
  dplyr::mutate(dplyr::select(., dplyr::starts_with('theta:')) %>% 
                  apply(., 1, quantileSummary, CI = 0.95) %>% 
                  lapply(., data.frame) %>% 
                  do.call(rbind, .),
                lower95 = lower,
                upper95 = upper) %>%
  dplyr::summarise(intervalScore50 = intervalScore(lower = lower50, upper = upper50, true = logit.est, alpha = 0.50)$score,
                   coverage50 = intervalScore(lower = lower50, upper = upper50, true = logit.est, alpha = 0.50)$coverage*100,
                   intervalScore95 = intervalScore(lower = lower95, upper = upper95, true = logit.est, alpha = 0.05)$score,
                   coverage95 = intervalScore(lower = lower95, upper = upper95, true = logit.est, alpha = 0.05)$coverage*100, 
                   .by = 'type'); adm1.cv.dist.scores

adm1.cv.all.scores <-
  adm1.cv.point.scores %>% 
  dplyr::left_join(., adm1.cv.dist.scores, by = 'type') %>% 
  dplyr::arrange(type); adm1.cv.all.scores

print(xtable::xtable(adm1.cv.all.scores, digits = 1), include.rownames = FALSE)

## scores plot -----------------------------------------------

setwd(res.dir)

adm1.cv.point.score.boxplot.data <-
  adm1.cv %>% 
  dplyr::mutate(mae.dist = 
                  dplyr::across(dplyr::starts_with('theta:'), ~ abs(.x - logit.est)) %>% 
                  rowMeans(),
                mse.dist = 
                  dplyr::across(dplyr::starts_with('theta:'), ~ (.x - logit.est)^2) %>% 
                  rowMeans()) %>%
  dplyr::select(-starts_with('theta:'), 
                -logit.est, -logit.prec, -mean, -variance, -median, -upper, -lower) %>% 
  tidyr::pivot_longer(cols = c('mae.dist', 'mse.dist'),
                      names_to = 'metric', values_to = 'value') %>% 
  dplyr::mutate(metric = 
                  metric %>% 
                  factor(., 
                         levels = c('mae.dist', 'mse.dist'),
                         labels = c('MAE', 'MSE')))

adm1.cv.point.score.boxplot <- 
  ggplot2::ggplot(data = adm1.cv.point.score.boxplot.data, aes(x = type, y = value, shape = type, colour = type)) +
  ggplot2::geom_boxplot(outlier.shape = NA) +
  ggplot2::geom_jitter() +
  ggplot2::facet_wrap(~ metric) +
  ggplot2::labs(y = '', x = '') +
  ggplot2::scale_colour_manual(values = c('red3', 'blue3', 'green4')) + 
  ggplot2::scale_y_continuous(breaks = seq(from = 0, to = 3, by = 0.5)) +
  ggplot2::geom_hline(yintercept = 0, linetype = 'dashed') +
  plotTheme(legend.title = element_blank(),
            axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
            text = element_text(size = text.size)); adm1.cv.point.score.boxplot

ggplot2::ggsave(plot = adm1.cv.point.score.boxplot,
                filename = paste0('Figures/', country, '_adm1_cv_point_score_boxplot.png'),
                width = width, height = height)

## distribution ridge plot -----------------------------------------------

setwd(res.dir)

adm1.u5.cv.ridgeplot.data <- 
  adm1.cv %>% 
  dplyr::mutate(dplyr::across(dplyr::starts_with('theta:'), ~ 1000*SUMMER::expit(.x))) %>% 
  dplyr::left_join(., admin1.names, by = c('region' = 'Internal')) %>% 
  tidyr::pivot_longer(cols = dplyr::starts_with('theta:'), names_to = 'theta', values_to = 'value') %>% 
  dplyr::filter(!is.na(logit.est)) %>% 
  dplyr::arrange(median) %>% 
  dplyr::mutate(GADM = GADM %>% haven::as_factor())

adm1.u5.cv.ridgeplot.direct.data <-
  adm1.u5.cv.ridgeplot.data %>% 
  dplyr::select(GADM, logit.est) %>% 
  dplyr::distinct() %>% 
  dplyr::mutate(GADM = GADM %>% as.numeric(),
                logit.est = 1000*SUMMER::expit(logit.est),
                number  = dplyr::row_number())

adm1.u5.cv.ridgeplot <- 
  ggplot2::ggplot(adm1.u5.cv.ridgeplot.data, 
                  aes(x = value, y = GADM, fill = after_stat(x))) +
  ggridges::geom_density_ridges_gradient(scale = 3, size = 0.3, rel_min_height = 0.01) +
  ggplot2::scale_fill_viridis_c(name = 'Predicted U5MR \n(per 1000 live births)', 
                                direction = -1, option = 'G', 
                                limits = c(0, 300)) +
  ggplot2::geom_segment(data = adm1.u5.cv.ridgeplot.direct.data,
                        aes(y = GADM, yend = GADM + 1, x = logit.est, xend = logit.est),
                        color = 'red', linewidth = 1) +
  ggplot2::xlim(0, 250) +
  ggplot2::scale_y_discrete(limits = rev) +
  ggplot2::labs(y = '', x = '') +
  ggplot2::guides(fill = guide_colourbar(title.position = 'top', title.hjust = 0.5)) +
  ggplot2::facet_grid(~ type) +
  plotTheme(legend.position = 'bottom',
            text = element_text(size = text.size)); adm1.u5.cv.ridgeplot

ggplot2::ggsave(plot = adm1.u5.cv.ridgeplot,
                filename = paste0('Figures/', country, '_adm1_u5_cv_ridgeplot.png'),
                width = 1.25*width, height = 1.25*height)

## distribution map plots -----------------------------------------------

setwd(res.dir)

adm1.u5.cv.map.data <-
  adm1.cv %>% 
  dplyr::mutate(dplyr::select(., dplyr::starts_with('theta:')) %>%
                  SUMMER::expit(.) %>% 
                  apply(., 1, quantileSummary) %>% 
                  lapply(., data.frame) %>% 
                  do.call(rbind, .)) %>%
  dplyr::select(-dplyr::starts_with('theta:'), -logit.est, -logit.prec) %>%  
  dplyr::mutate(years = years %>% as.numeric(),
                lower = lower * 1000,
                median = median * 1000,
                upper = upper * 1000,
                width = upper - lower) %>% 
  dplyr::left_join(., admin1.names, by = c('region' = 'Internal')) %>% 
  dplyr::left_join(., poly.adm1, by = c('GADM' = 'NAME_1')) %>%
  sf::st_as_sf()

adm1.u5.cv.map.estimate.plot <-
  ggplot2::ggplot() +
  ggplot2::geom_sf(data = adm1.u5.cv.map.data, aes(fill = median, group = type), colour = 'black') +
  ggplot2::scale_fill_viridis_c(name = 'Predicted U5MR (death per 1000 live births)',
                                direction = -1, option = 'G', n.breaks = 4) +
  ggplot2::facet_grid(~ type) +
  mapTheme(text = element_text(size = text.size),
           legend.position = 'bottom') +
  ggplot2::guides(fill = guide_colourbar(title.position = 'top', title.hjust = 0.5));  #adm1.u5.cv.map.estimate.plot

adm1.u5.cv.map.width.plot <-
  ggplot2::ggplot() +
  ggplot2::geom_sf(data = adm1.u5.cv.map.data, aes(fill = width), colour = 'black') +
  ggplot2::scale_fill_viridis_c(name = 'Width of Credible Interval',
                                direction = -1, option = 'A', n.breaks = 3) + 
  ggplot2::facet_grid(~ type) +
  mapTheme(text = element_text(size = text.size),
           legend.position = 'bottom') +
  ggplot2::guides(fill = guide_colourbar(title.position = 'top', title.hjust = 0.5)); #adm1.u5.cv.map.width.plot

ggplot2::ggsave(plot = adm1.u5.cv.map.estimate.plot,
                filename = paste0('Figures/', country, '_adm1_u5_cv_map_estimate_plot.png'),
                height = height, width = width)

ggplot2::ggsave(plot = adm1.u5.cv.map.width.plot,
                filename = paste0('Figures/', country, '_adm1_u5_cv_map_width_plot.png'),
                height = height, width = width)

## cluster level versus age-period -----------------------------------------------

setwd(res.dir)

# BB estimates
load(file = paste0('Direct/U5MR/',country, '_res_natl_', time.model, '_yearly_u5_SmoothedCluster.rda'))
# ap
load(file = paste0('Age-Period/U5MR/',country,'_res_natl_u5.rda'))

natl.frame.comp <-
  rbind(bb.natl.yearly.u5$overall %>% 
          dplyr::select(years, lower, median, upper) %>% 
          dplyr::mutate(type = 'Cluster-Level',
                        years = years %>% as.character %>% as.numeric()),
        ap.natl.u5.res$national %>% 
          dplyr::select(years, lower, median, upper) %>% 
          dplyr::mutate(type = 'Age-Period')) %>% 
  dplyr::mutate(years = years %>% as.numeric(),
                lower = 1000 * lower,
                median = 1000 * median,
                upper = 1000 * upper,
                type = type %>% factor(., levels = c('Cluster-Level', 'Age-Period'))) %>% 
  dplyr::filter(years %in% beg.year:end.year)

natl.line.comp.plot <- 
  ggplot2::ggplot(data = natl.frame.comp, 
                  aes(x = years, y = median, colour = type, group = type, fill = type)) +
  ggplot2::geom_line() +
  ggplot2::geom_point() +
  ggplot2::scale_x_continuous(breaks = seq(from = beg.year, to = end.year, by = 1)) +
  ggplot2::scale_colour_manual(values = c('red3', 'blue3')) + 
  ggplot2::labs(x = 'Period', y = 'U5MR (deaths per 1000 live births)') +
  plotTheme(text = element_text(size = text.size),
            axis.text.x = element_text(angle = 90, hjust = 1),
            legend.title = element_blank(),
            legend.position = 'top'); natl.line.comp.plot

natl.error.comp.plot <-
  ggplot2::ggplot(data = natl.frame.comp, 
                  aes(x = years, y = median, colour = type, group = type, fill = type)) +
  ggplot2::geom_point(position = position_dodge(0.5)) +
  ggplot2::geom_errorbar(aes(ymin = lower, ymax = upper), width = .2, position = position_dodge(0.5)) +
  ggplot2::scale_x_continuous(breaks = seq(from = beg.year, to = end.year, by = 1)) +
  ggplot2::scale_colour_manual(values = c('red3', 'blue3')) + 
  ggplot2::labs(x = 'Period', y = 'U5MR (deaths per 1000 live births)') +
  plotTheme(text = element_text(size = text.size),
            axis.text.x = element_text(angle = 90, hjust = 1),
            legend.title = element_blank(),
            legend.position = 'top'); natl.error.comp.plot

ggplot2::ggsave(plot = natl.line.comp.plot,
                filename = paste0('Figures/', country, '_natl_line_comp_plot.png'),
                height = height, width = width)

ggplot2::ggsave(plot = natl.error.comp.plot,
                filename = paste0('Figures/', country, '_natl_error_comp_plot.png'),
                height = height, width = width)
