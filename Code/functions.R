# precision matric for ICAR
precision.icar <- function(amat){
  
  # structured
  Q <- amat
  
  # organise amat into ICAR precision
  if (sum(Q > 0 & Q < 1) != 0) {
    #### if 0 < Q_{ij} < 1 then make Q_{ij} = 1
    for (i in 1:nrow(Q)) {
      idx <- which(Q[i, ] > 0 & Q[i, ] < 1)
      Q[i, idx] <- 1
    }
  }
  
  # turn amat into an ICAR precision
  ## no 0s on diagonal
  ## main diagonal to be the number of neighbours in total
  ## off diagonal to be 0 for not neighbours and -1 for neighbours
  diag(Q) <- 0
  diag <- apply(Q, 1, sum)
  Q[Q != 0] <- -1
  diag(Q) <- diag
  
  # scale and set constraints
  
  ## inla functions
  inla.scale.model.bym <- utils::getFromNamespace('inla.scale.model.bym', 'INLA')
  inla.bym.constr.internal <- utils::getFromNamespace('inla.bym.constr.internal', 'INLA')
  
  ## scale
  Q <- inla.scale.model.bym(Q, adjust.for.con.comp = TRUE)
  
  ## set constraints
  Q.constr <- inla.bym.constr.internal(Q, adjust.for.con.comp = TRUE)
  
  return(Q)
  
}

# APC model
smoothAPC_new <- function(data, admin.level, Amat = NULL, 
                      age.group = c("0", "1-11", "12-23", "24-35", "36-47", "48-59"), 
                      year.label, stratified = FALSE, 
                      mod = c('apc', 'ac', 'ap', 'pc', 'a', 'p', 'c')[1], slope.drop = NULL, 
                      family = c('betabinomial', 'binomial')[1],
                      pc.u = 1, pc.alpha = 0.01, pc.u.phi = 0.5, pc.alpha.phi = 2/3,
                      overdisp.mean = 0, overdisp.prec = 0.4, 
                      type.st = 4, pc.st.u = NA, pc.st.alpha = NA,
                      control.inla = list(strategy = 'adaptive', int.strategy = 'auto'),
                      inla.mode = c('classic', 'twostage', 'experimental')[3],
                      control.compute = list(config = TRUE), verbose = FALSE, ...){
  
  # inputs (for testing) ----
  
  # data = births.2014
  # admin.level = 'Admin1'
  # Amat = admin1.mat
  # age.groups = c("0", "1-11", "12-23", "24-35", "36-47", "48-59")
  # year.label = c(beg.year:end.proj.year)
  # mod = 'apc'
  # slope.drop = 'c'
  # stratified = TRUE
  # family = c('betabinomial', 'binomial')[1]
  # pc.u = 1
  # pc.alpha = 0.01
  # pc.u.phi = 0.5
  # pc.alpha.phi = 2/3
  # overdisp.mean = -7.5
  # overdisp.prec = 0.39
  # type.st = 4
  # pc.st.u = NA
  # pc.st.alpha = NA
  # control.inla = list(strategy = 'adaptive', int.strategy = 'auto')
  # inla.mode = c('classic', 'twostage', 'experimental')[3]
  # control.compute = list(config = TRUE)
  # verbose = FALSE
  
  # check inputs ----
  
  if(!(admin.level %in% c('National','Admin1','Admin2'))){
    stop('Enter a valid admin level (National, Admin1, or Admin2)')
  }
  if((admin.level!='National') == (is.null(Amat) == T)){
    stop('Admin level and adjacency matrix are not compatible')}
  if(!"config" %in% names(control.compute)){
    message("config = TRUE is added to control.compute so that posterior draws can be taken.")
    control.compute$config <- TRUE
  }
  if(mod == 'apc' && is.null(slope.drop)){
    stop('Warning: When fitting an APC model need to drop linear column of one of age, period or cohort.')
  }
  if(mod == 'apc' &&! (slope.drop %in% c('a', 'p', 'c'))){
    stop('slope.drop needs to be one of a, p or c')
  }
  if (!is.null(Amat)) {
    if (is.null(rownames(Amat))) {
      stop('Row names of Amat needs to be specified to region names.')
    }
    if (is.null(colnames(Amat))) {
      stop('Column names of Amat needs to be specified to region names.')
    }
    if (sum(rownames(Amat) != colnames(Amat)) > 0) {
      stop('Row and column names of Amat needs to be the same.')
    }
    is.spatial <- TRUE
  } else {
    Amat <- matrix(1, 1, 1)
    colnames(Amat) <- rownames(Amat) <- "All"
    is.spatial <- FALSE
  }
  
  # set admin-level parameters ----
  
  if(admin.level == 'National'){
    data$region <- "All"
  }else if(admin.level == 'Admin1'){
    data$region <- data$admin1.char
  }else if(admin.level == 'Admin2'){
    data$region <- data$admin2.char
  }
  
  # data oragnisation ----
  
  ## correct temporal specification ----
  
  # age + period + cohort
  data2a <- 
    expand.grid(ageMonth = 0:59,
                periodMonth = (12*(year.label[1] - 1900) + 1):(12*(year.label[length(year.label)] - 1900) + 12)) %>% 
    dplyr::mutate(cohortMonth = periodMonth - ageMonth,
                  age = dplyr::case_when(ageMonth == 0 ~ 0,
                                         ageMonth %in% 1:11 ~ 6,
                                         ageMonth %in% 12:23 ~ 17.5,
                                         ageMonth %in% 24:35 ~ 29.5,
                                         ageMonth %in% 36:47 ~ 41.5,
                                         TRUE ~ 53.5),
                  period = floor((periodMonth-1)/12)+1900,
                  cohort= floor((cohortMonth-1)/12)+1900,
                  age_id = age,
                  period_id = period %>% factor() %>% as.numeric(),
                  cohort_id = cohort %>% factor() %>% as.numeric()) %>% 
    dplyr::select(-tidyr::ends_with('Month')) %>% 
    dplyr::distinct()
  A <- data2a$age %>% unique() %>% length()
  P<- data2a$period %>% unique() %>% length()
  C <- data2a$cohort %>% unique() %>% length()
  
  ## correct spatial specification ----
  
  R <- nrow(Amat)
  data2b <- data.frame(space = rownames(Amat), space_id = 1:R)
  
  ## correct spatio-temporal specification ----
  
  ### space + age + period + cohort
  data2c <- 
    dplyr::cross_join(data2a, data2b) %>% 
    dplyr::mutate(spaceTime = interaction(space, period),
                  spaceTime_id = interaction(space_id, period_id) %>% as.numeric()) %>% 
    dplyr::relocate(c('space', 'spaceTime'), .before = age_id)
  
  ## finalised data ----
  
  data3 <- 
    data  %>%
    # labels for INLA modelling
    dplyr::mutate(
      # new terms
      age = dplyr::case_when(age == '0' ~ 0,
                             age == '1-11' ~ 6,
                             age == '12-23' ~ 17.5,
                             age == '24-35' ~ 29.5,
                             age == '36-47' ~ 41.5,
                             TRUE ~ 53.5),
      period = years,
      space = region,
      strata = urban) %>% 
    # join with the full range of space and time
    dplyr::left_join(., data2c, by = c('age', 'period', 'cohort', 'space')) %>% 
    dplyr::mutate(
      # final id terms for inla
      cluster_id = cluster %>% as.factor() %>% as.numeric(),
      strata_id = strata %>% factor(., levels = c('rural', 'urban')), # want this to be a factor in the model as we stratify by it
      age2_id = age_id, 
      period2_id = period_id, 
      cohort2_id = cohort_id)
  
  # priors ----
  
  # pc hyper priors
  ## temproal
  hyper_pc_time <- list(prec = list(prior = 'pc.prec', param = c(pc.u, pc.alpha)))
  ## spatial
  hyper_pc_space <- list(prec = list(prior = 'pc.prec', param = c(pc.u, pc.alpha)), phi = list(prior = 'pc', param = c(pc.u.phi, pc.alpha.phi)))
  ## spatio-temproal
  pc.st.u <- ifelse(is.na(pc.st.u), pc.u, pc.st.u)
  pc.st.alpha <- ifelse(is.na(pc.st.alpha), pc.alpha, pc.st.alpha)
  hyper_pc_spaceTime <- list(prec = list(prior = "pc.prec", param = c(pc.st.u, pc.st.alpha)))
  
  # formula ----
  
  ## temporal (APC) ----
  
  ### temporal constraints ----
  
  # age
  age.values <- data3$age_id %>% unique() %>% sort()
  age.rank.def <- 2
  age.cMat <- rbind(1, age.values) 
  age.e <- rep(0, times = age.rank.def)
  # period
  per.values <- data3$period_id %>% unique() %>% sort()
  per.rank.def <- 2
  per.cMat <- rbind(1, per.values) 
  per.e <- rep(0, times = per.rank.def)
  # cohort
  coh.values <- data3$cohort_id %>% unique() %>% sort()
  coh.rank.def <- 2
  coh.cMat <- rbind(1, coh.values) 
  coh.e <- rep(0, times = coh.rank.def)
  
  ### APC selection
  
  # full temporal formula
  formula <- 
    Y ~ age_id + period_id + cohort_id +
    f(age2_id, model = 'rw2', hyper = hyper_pc_time, scale.model = TRUE, constr = FALSE, rankdef = age.rank.def, extraconstr = list(A = age.cMat, e = age.e, values = age.values)) +
    f(period2_id, model = 'rw2', hyper = hyper_pc_time, scale.model = TRUE, constr = FALSE, rankdef = per.rank.def, extraconstr = list(A = per.cMat, e = per.e, values = per.values)) +
    f(cohort2_id, model = 'rw2', hyper = hyper_pc_time, scale.model = TRUE, constr = FALSE, rankdef = coh.rank.def, extraconstr = list(A = coh.cMat, e = coh.e, values = coh.values))
  
  message("----------------------------------",
          "\nModel Specification", 
          "\n Temporal effect(s):               ",
          "\n  APC type:                    ",
          mod,
          appendLF = FALSE)
  # update pre and post forumla for the temporal models
  # # do not remove any linear terms from pre fit as they are all needed later in 
  # # re-parameterisation
  if (mod == 'apc'){
    message("\n  Slope dropped:                 ", slope.drop, appendLF = FALSE)
    if(slope.drop == 'a'){
      formula <- update(formula, ~. - age_id)
    } else if (slope.drop == 'p'){
      formula <- update(formula, ~. - period_id)
    } else if (slope.drop == 'c'){
      formula <- update(formula, ~. - cohort_id)
    }
  } else if (mod == 'ac'){
    formula <- update(formula, ~.
                      - period_id 
                      - f(period2_id, model = 'rw2', hyper = hyper_pc_time, scale.model = TRUE, constr = FALSE, rankdef = per.rank.def, extraconstr = list(A = per.cMat, e = per.e, values = per.values)))
  } else if (mod == 'pc'){
    formula <- update(formula, ~.
                      - age_id 
                      - f(age2_id, model = 'rw2', hyper = hyper_pc_time, scale.model = TRUE, constr = FALSE, rankdef = age.rank.def, extraconstr = list(A = age.cMat, e = age.e, values = age.values)))
  } else if (mod == 'ap'){
    formula <- update(formula, ~.
                      - cohort_id 
                      -  f(cohort2_id, model = 'rw2', hyper = hyper_pc_time, scale.model = TRUE, constr = FALSE, rankdef = coh.rank.def, extraconstr = list(A = coh.cMat, e = coh.e, values = coh.values)))
  }
  
  ## spatial & spatio-temporal ----
  
  if (is.spatial) {
    message("\n Spatial effect(s):               ",
            "\n  Structured:                 bym2",
            "\n  Interaction type:              ", type.st,
            appendLF = FALSE)
    
    ### spatial ----
    
    formula <- 
      update(formula, ~. +
               f(space_id, graph = Amat, model = 'bym2', hyper= hyper_pc_space, scale.model = TRUE, adjust.for.con.comp = TRUE))
    
    
    ### spatio-temporal ----
    
    inla.rw <- utils::getFromNamespace('inla.rw', 'INLA')
    order.rw <- 2
    
    if (type.st == 1) {
      
      # kronecker product for the Type I Interaction
      ## IID x IID
      formula <- update(formula, ~. + f(spaceTime_id, model="iid", hyper = hyper_pc_spaceTime))
      
    } else {
      
      if (type.st == 2) {
        # unstructured space x structured time
        Q.space <- diag(R)
        Q.time <- inla.rw(n = P, order = order.rw, scale.model = TRUE, sparse = TRUE)
        # constraint
        A1 <- kronecker(matrix(data = 1, nrow = 1, ncol = P), Matrix::Diagonal(R))
        A.constr <- A1[-1, ] %>% as.matrix()
        # rank
        full.rank <- (R)*(P)
        rank <- (R)*(P - order.rw)
        rank.def <- full.rank - rank
      } else if (type.st == 3) {
        # structured space x unstructured time
        Q.space <- precision.icar(Amat)
        Q.time <- diag(P)
        # constraint
        A1 <- kronecker(Matrix::Diagonal(P), matrix(data = 1, nrow = 1, ncol = R))
        A.constr <- A1[-1, ] %>% as.matrix()
        # rank
        full.rank <- (R)*(P)
        rank <- (R - 1)*(P)
        rank.def <- full.rank - rank
      } else {
        # structured space x structured time
        Q.space <- precision.icar(Amat)
        Q.time <- inla.rw(n = P, order = order.rw, scale.model = TRUE, sparse = TRUE)
        # constraint
        A1 <- kronecker(matrix(data = 1, nrow = 1, ncol = P), Matrix::Diagonal(R))
        A2 <- kronecker(Matrix::Diagonal(P), matrix(data = 1, nrow = 1, ncol = R))
        A.constr <- rbind(A1[-1, ], A2[-1, ]) %>% as.matrix()
        # rank
        full.rank <- (R)*(P)
        rank <- (R - 1)*(P - order.rw)
        rank.def <- full.rank - rank
      }
      
      Q_spaceTime <- kronecker(Q.time, Q.space)
      
      formula <- update(formula, ~. +
                          f(spaceTime_id,
                            model = "generic0",
                            Cmatrix = Q_spaceTime,
                            rankdef = rank.def,
                            constr = TRUE,
                            extraconstr = list(A = A.constr, e = rep(0, nrow(A.constr))),
                            hyper = hyper_pc_spaceTime))
      
    }
    
  }
  
  ## cluster (if binomial) ----
  
  # cluster random effect for binomial model
  if (family == 'binomial') {
    message("\n  cluster:                     iid", appendLF = FALSE)
    formula <- update(formula, ~. +  f(cluster_id, model = 'iid', hyper = hyper_pc_time))
  }
  
  ## stratification ----
  
  # sorting out the fixed effects
  ## urban/rural strata
  if(stratified){
    message("\n  Stratified:                  yes", appendLF = FALSE)
    message("\n  Global Intercept:             no", appendLF = FALSE)
    # update the formula
    formula <- update(formula, ~. - 1 + strata_id)
  } else {
    message("\n  Stratified:                   no", appendLF = FALSE)
    message("\n  Global intercept:            yes", appendLF = FALSE)
    data3$strata <- as.factor('All')
    data3$strata_id <- data3$strata_id <- 1
  }
  
  # model fit ----
  
  message("\nFitting model...",
          appendLF = FALSE)
  startClock <- proc.time() # Start the clock!
  # the extra information depending on the model
  if(family == 'betabinomial'){
    control.family <- list(hyper = list(rho = list(param = c(overdisp.mean, overdisp.prec), initial = overdisp.mean)))
    fit <-
      INLA::inla(formula, 
                 family = family, 
                 data = data3, 
                 Ntrials = data3$total,
                 control.compute = control.compute,
                 control.family = control.family,
                 control.predictor = list(compute = FALSE, link = 1),
                 control.inla = control.inla,
                 lincomb = NULL,
                 inla.mode = inla.mode,
                 verbose = verbose,
                 ...)
  } else {
    fit <-
      INLA::inla(formula, 
                 family = family, 
                 data = data3, 
                 Ntrials = data3$total,
                 control.compute = control.compute,
                 control.predictor = list(compute = FALSE, link = 1),
                 control.inla = control.inla,
                 lincomb = NULL,
                 inla.mode = inla.mode,
                 verbose = verbose,
                 ...)
  }
  endClock <- proc.time() - startClock # Stop the clock
  message("\n Time take to fit model(s):     ", endClock[3] %>% as.numeric %>% round(),
          "\n----------------------------------",
          appendLF = FALSE)
  
  # output ----
  
  priors <- list(pc.u = pc.u, pc.alpha = pc.alpha, pc.u.phi = pc.u.phi, 
                 pc.alpha.phi = pc.alpha.phi, 
                 pc.st.u = pc.st.u, 
                 pc.st.alpha = pc.st.alpha, 
                 overdisp.mean = overdisp.mean, overdisp.prec = overdisp.prec)
  
  return(list(model = formula, fit = fit, family = family, 
              Amat = Amat, newdata = data3, type.st = type.st,
              age.group = age.group, year.label = year.label,
              priors = priors))
}

# posterior sampling
getSmoothedAPC <- function(inla.mod, n.sim = 1000, 
                           use.approximation = TRUE, 
                           verbose = FALSE, CI = 0.95, 
                           save.u5m.draws = FALSE, 
                           strata.weights = NULL,
                           save.linpred.draws = FALSE, ...){
  
  # inla.mod = apc.adm1.u5.fit
  # strata.weights = NULL
  # use.approximation = TRUE
  # n.sim = 10
  # verbose = FALSE
  # CI = 0.95
  # save.u5m.draws = TRUE
  # save.linpred.draws = TRUE
  # strata.weights = weight.strata.natl.u5
  
  ## sample from posterior ----
  
  message("-----------------------------------",
          "\nStarting posterior sampling...   ",
          appendLF = F)
  start1 <- proc.time()
  sampAll <- INLA::inla.posterior.sample(n = n.sim, result = inla.mod$fit, intern = TRUE, verbose = verbose, ...)
  end1 <- proc.time() - start1
  
  start2 <- proc.time()
  
  ## define prediciton data ----
  ### age + period + cohort
  data2a <- 
    expand.grid(ageMonth = 0:59,
                periodMonth = (12*(inla.mod$year.label[1] - 1900) + 1):(12*(inla.mod$year.label[length(inla.mod$year.label)] - 1900) + 12)) %>% 
    dplyr::mutate(cohortMonth = periodMonth - ageMonth,
                  age = dplyr::case_when(ageMonth == 0 ~ 0,
                                         ageMonth %in% 1:11 ~ 6,
                                         ageMonth %in% 12:23 ~ 17.5,
                                         ageMonth %in% 24:35 ~ 29.5,
                                         ageMonth %in% 36:47 ~ 41.5,
                                         TRUE ~ 53.5),
                  period = floor((periodMonth-1)/12)+1900,
                  cohort= floor((cohortMonth-1)/12)+1900,
                  age_id = age, 
                  period_id = period %>% factor() %>% as.numeric(),
                  cohort_id = cohort %>% factor() %>% as.numeric(),
                  age2_id = age_id %>% as.factor() %>% as.numeric(), # AT THIS STAGE WANT THIS AS 1:A
                  period2_id = period_id, 
                  cohort2_id = cohort_id) %>% 
    dplyr::select(-tidyr::ends_with('Month')) %>% 
    dplyr::distinct()
  A <- data2a$age %>% unique() %>% length()
  P<- data2a$period %>% unique() %>% length()
  C <- data2a$cohort %>% unique() %>% length()
  ### space 
  R <- nrow(inla.mod$Amat)
  data2b <- 
    data.frame(space = rownames(inla.mod$Amat), space_id = 1:R)
  ### strata
  data2c <- 
    inla.mod$newdata %>% 
    dplyr::select(strata, strata_id) %>% 
    dplyr::distinct() %>% 
    # need this to include intercepts correctly
    dplyr::mutate(strata_idurban = dplyr::case_when(strata == 'urban' ~ 1, TRUE ~ 0),
                  strata_idrural = dplyr::case_when(strata == 'rural' ~ 1, TRUE ~ 0))
  ### join all together
  data3 <- 
    # age + period + cohort + space
    dplyr::cross_join(data2a, data2b) %>%
    dplyr::mutate(spaceTime = interaction(space, period),
                  spaceTime_id = interaction(space_id, period_id) %>% as.numeric()) %>%
    # ~ + strata
    dplyr::cross_join(., data2c) %>% 
    # order 
    dplyr::relocate(c('space', 'spaceTime', 'strata'), .before = age_id)
  
  ### locations
  fields <- rownames(sampAll[[1]]$latent) # to get the names for labelling
  
  # locations
  data3.loc <- 
    data3 %>%
    # chnage the variable names to match INLA output names
    mutate(intercept = '(Intercept):1',
           strata_idurban = 'strata_idurban:1',
           strata_idrural = 'strata_idrural:1',
           age_id = 'age_id:1',
           period_id = 'period_id:1',
           cohort_id = 'cohort_id:1') %>% 
    # match the names to INLA output positions
    mutate(intercept = match(intercept, fields),
           strata_idurban = match(strata_idurban, fields),
           strata_idrural = match(strata_idrural, fields),
           age_id = match(age_id, fields),
           period_id = match(period_id, fields),
           cohort_id = match(cohort_id, fields),
           age2_id = match(paste0('age2_id:', age2_id), fields),
           period2_id = match(paste0('period2_id:', period2_id), fields),
           cohort2_id = match(paste0('cohort2_id:', cohort2_id), fields),
           space_id = match(paste0('space_id:', space_id), fields),
           spaceTime_id = match(paste0('spaceTime_id:', spaceTime_id), fields))
  
  ## linear predictor ----
  theta <- matrix(0, nrow = n.sim, ncol = dim(data3)[1])
  for (i in 1:n.sim) {
    
    # i <- 1
    
    # extract sample i 
    draw <- sampAll[[i]]$latent
    
    # linear predictor
    ## drop the columns for the fixed effects
    ## only want samples from the random effect columns
    theta[i, ] <- 
      apply(data3.loc %>% 
              dplyr::select(age2_id, period2_id, cohort2_id, space_id, spaceTime_id), 
            1, 
            function(x, ff) { sum(ff[x], na.rm = TRUE) }, 
            draw)
    
    # now include the fixed effects
    ## want to have the slope multiplied by the value
    add.intercept <- draw[data3.loc$intercept]
    add.urban.slope <- draw[data3.loc$strata_idurban] * data3$strata_idurban
    add.rural.slope <- draw[data3.loc$strata_idrural] * data3$strata_idrural
    add.age.slope <- draw[data3.loc$age_id] * data3$age_id
    add.period.slope <- draw[data3.loc$period_id] * data3$period_id
    add.cohort.slope <- draw[data3.loc$cohort_id] * data3$cohort_id
    
    ## if there is any NA return give 0
    ### NA will be return for the linear slope that is dropped
    add.intercept[is.na(add.intercept)] <- 0
    add.urban.slope[is.na(add.urban.slope)] <- 0
    add.rural.slope[is.na(add.rural.slope)] <- 0
    add.age.slope[is.na(add.age.slope)] <- 0
    add.period.slope[is.na(add.period.slope)] <- 0
    add.cohort.slope[is.na(add.cohort.slope)] <- 0
    
    theta[i, ] <- theta[i, ] + add.intercept + add.urban.slope + add.rural.slope + add.age.slope + add.period.slope + add.cohort.slope
    if (inla.mod$family == 'binomial') {
      tau[i] <- exp(sampAll[[i]]$hyperpar[["Log precision for cluster_id"]]) # tau = exp(log(precision))
    }
  }
  
  ## process results ----
  message("\nCleaning the results...          ",
          appendLF = F)
  
  
  # function for summary
  quantileSummary <- function(x, CI = 0.95) {
    
    # x = theta1[,1]; lowerCI = lowerCI; upperCI = upperCI
    
    # lower and upper CI calculators
    lowerCI <- (1 - CI)/2
    upperCI <- 1 - lowerCI
    
    qntl <- quantile(x, probs = c(lowerCI, 0.5, upperCI))
    data.frame(mean = mean(x), variance = var(x), lower = qntl[1], median = qntl[2], upper = qntl[3])
    
  }
  
  # marginalisation function
  ## for binomial distribution
  ## use approximation or numerical integration
  marginalise.logit.binomial = function(tauTheta, use.approximation = TRUE, ...) {
    
    # tauTheta = cbind(tau, theta)[1,]; use.approximation = FALSE
    
    # tau = 1/sigma^2
    sigma <- sqrt(1/tauTheta[1])
    theta <- as.matrix(tauTheta[-1], ncol = 1)
    
    # calculation marginalistion
    if(sigma == 0){
      theta
    } else {
      if(any(is.na(c(theta, sigma)))){
        NA
      } else if(!use.approximation) {
        # numerically calculate the mean
        ## integrand to calculate numerical marginal
        integrand <- function(x, theta, sigma){ exp(plogis(x, log.p=TRUE) + dnorm(x, mean = theta, sd = sigma, log=TRUE)) }
        ## function to find integral
        ### requires theta and sigma to be single elements
        my.integrate <- function(fun, theta, sigma, ...){
          integrate(fun, lower = theta-10*sigma, upper = theta+10*sigma, theta = theta, sigma = sigma, abs.tol = 0, ...)$value
        }
        ## apply on my.integrate
        ### will use 1 sigma but cycle through thetas defined above
        sapply(X = theta, FUN = my.integrate, fun = integrand, sigma = sigma)
      } else {
        # use logistic approximation
        k = 16 * sqrt(3) / (15 * pi)
        theta / sqrt(1 + k^2 * sigma^2)
      }
    }
  }
  
  # correcting post marginalisation if binomial
  ## do we want to use approximation or numerical integration for marginalisation
  if (inla.mod$family == 'binomial') {
    if(use.approximation){
      message("\n Binomial marginalisation found using approximation", appendLF = F)
      theta2 <- apply(X = cbind(tau, theta), 1, FUN = marginalise.logit.binomial, use.approximation = TRUE)
    } else{
      message("\n Binomial marginalisation found using numerical integration", appendLF = F)
      theta2 <- apply(X = cbind(tau, theta), 1, FUN = marginalise.logit.binomial, use.approximation = FALSE)
    }
    
  } else {
    theta2 <- t(theta)
  }
  colnames(theta2) <- paste0('theta:', 1:ncol(theta2))
  theta3 <- 
    dplyr::bind_cols(data3 %>% 
                       dplyr::select(-tidyr::ends_with('_id')) %>% 
                       dplyr::rename(region = space, years = period),
                     theta2)
  
  # suppress summarise info
  options(dplyr.summarise.inform = FALSE)
  
  ### include strata weights ----
  if(!is.null(strata.weights)){
    if(nrow(inla.mod$Amat) == 1){ 
      strata.weigths2 <- 
        strata.weights %>% 
        dplyr::mutate(region = 'All') %>% 
        tidyr::pivot_longer(., cols = c('urban', 'rural'), names_to = 'strata', values_to = 'proportion')
      
      message("\nNational weights included", appendLF = F)
      
    } else {
      strata.weigths2 <- 
        strata.weights %>% 
        tidyr::pivot_longer(., cols = c('urban', 'rural'), names_to = 'strata', values_to = 'proportion')
      
      message("\nSubnational weights included", appendLF = F)
      
    }
    
    theta4 <- 
      theta3 %>% 
      dplyr::left_join(., strata.weigths2, by = c('strata', 'years', 'region'))
    
  } else {
    
    theta4 <-
      theta3 %>% 
      dplyr::mutate(proportion = 0)
    
    message("\n No strata weights supplied. Set all weights to zero", appendLF = F)
    
  } 
  
  ### labels from number of months in each age group ----
  labels <- inla.mod$age.group
  ns <- rep(1, length(labels))
  for (i in 1:length(labels)) {
    if (labels[i] == "0") {
      ns[i] <- 1
      next
    }
    tmp <- as.numeric(strsplit(as.character(labels[i]), "-")[[1]])
    ns[i] <- tmp[2] - tmp[1] + 1
  }
  theta5 <- 
    theta4 %>% 
    dplyr::mutate(ns = ns[age %>% as.factor() %>% as.numeric()]) %>% 
    dplyr::relocate(dplyr::starts_with('theta:'), .after = dplyr::last_col())
  
  
  ### u5mr formula ----
  # generic part of results
  ## true for overall and stratified
  ## true for national and subnational
  out.temp <- 
    theta5 %>% 
    dplyr::summarise(across(starts_with('theta:'), ~ SUMMER::expit(mean(.x))),
                     .by = c('region', 'strata', 'age', 'years', 'proportion', 'ns')) %>%
    # inside the product for survival
    dplyr::mutate(across(starts_with('theta:'), ~ ((1 - .x)^ns))) %>%
    # take product over all age groups
    dplyr::summarise(across(starts_with('theta:'), ~ prod(.x)),
                     .by = c('region', 'strata', 'years', 'proportion')) %>%
    # subtract the survival probability from one to get mortality
    dplyr::mutate(across(starts_with('theta:'), ~ 1 - .x))
  
  # stratified results
  ## by region and strata without proportions included
  out1 <- 
    out.temp %>%
    dplyr::mutate(dplyr::select(., dplyr::starts_with('theta:')) %>% 
                    apply(., 1, quantileSummary, CI = CI) %>% 
                    lapply(., data.frame) %>% 
                    do.call(rbind, .)) %>% 
    dplyr::relocate(dplyr::starts_with('theta:'), .after = dplyr::last_col())
  
  # subnational results
  ## proportionally summed over strata
  out2 <-
    out.temp %>%
    dplyr::mutate(across(starts_with('theta:'), ~ .x * proportion)) %>%
    # sum over proportions
    dplyr::summarise(across(starts_with('theta:'), sum),
                     .by = c('region', 'years')) %>%
    dplyr::mutate(dplyr::select(., dplyr::starts_with('theta:')) %>%
                    apply(., 1, quantileSummary, CI = CI) %>%
                    lapply(., data.frame) %>%
                    do.call(rbind, .)) %>% 
    dplyr::relocate(dplyr::starts_with('theta:'), .after = dplyr::last_col())
  
  # national results
  ## proportionally summed over strata
  out3 <-
    out.temp %>%
    dplyr::mutate(across(starts_with('theta:'), ~ .x * proportion),
                  region = 'All') %>%
    # sum over proportions
    dplyr::summarise(across(starts_with('theta:'), sum),
                     .by = c('region', 'years')) %>%
    dplyr::mutate(dplyr::select(., dplyr::starts_with('theta:')) %>%
                    apply(., 1, quantileSummary, CI = CI) %>%
                    lapply(., data.frame) %>%
                    do.call(rbind, .)) %>% 
    dplyr::relocate(dplyr::starts_with('theta:'), .after = dplyr::last_col())
  
  out <- list(stratified = out1, subnational = out2, national = out3)
  
  if(save.linpred.draws){
    out$linpred.draws <- theta3
  }
  
  if(save.u5m.draws){
    out$u5m.draws <- out.temp
  }
  
  end2 <- proc.time() - start2
  
  message("\nTime taken for posterior sampling:", 
          "\n Sampling (s):",  end1[3] %>% as.numeric %>% round(),
          "\n Cleaning (s):",  end2[3] %>% as.numeric %>% round(),
          "\n Total    (s):",  (end1[3] + end2[3]) %>% as.numeric %>% round(),
          "\n----------------------------------",
          appendLF = FALSE)
  
  return(out)
  
}

# function for summary
quantileSummary <- function(x, CI = 0.95) {
  
  # x = theta1[,1]; lowerCI = lowerCI; upperCI = upperCI
  
  # lower and upper CI calculators
  lowerCI <- (1 - CI)/2
  upperCI <- 1 - lowerCI
  
  qntl <- quantile(x, probs = c(lowerCI, 0.5, upperCI))
  data.frame(mean = mean(x), variance = var(x), lower = qntl[1], median = qntl[2], upper = qntl[3])
  
}

# interval score metrics
intervalScore <- function(lower, upper, true, alpha){
  
  # lower = adm1.cv.scores %>% filter(type == 'Age-Period-Cohort') %>% pull(lower)
  # upper = adm1.cv.scores %>% filter(type == 'Age-Period-Cohort') %>% pull(upper)
  # true = adm1.cv.scores %>% filter(type == 'Age-Period-Cohort') %>% pull(logit.est)
  # alpha = 0.05
  
  # each part seperately
  dispersion <- (upper - lower)
  overprediction <- 2/alpha * (lower - true) * (true < lower)
  underprediction <- 2/alpha * (true - upper) * (true > upper)
  
  score <- sum(dispersion, underprediction, overprediction, na.rm = TRUE) %>% mean()
  width <- dispersion %>% mean()
  coverage <- (lower <= true & true <= upper) %>% mean(., na.rm = TRUE)
  
  return(list(score = score, width = width, coverage = coverage))
}

# theme for plots
plotTheme<-function(...){
  theme(panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),
        panel.background=element_blank(),
        axis.line=element_line(colour='black'),
        # legend.title=element_blank(),
        legend.text = ggplot2::element_text(hjust = 0),
        legend.key=element_rect(fill=NA),
        ...)
}

# theme for map plots
mapTheme <- function(...){
  ggplot2::theme(axis.title.x = ggplot2::element_blank(),
                 axis.text.x = ggplot2::element_blank(),
                 axis.ticks.x = ggplot2::element_blank(),
                 axis.title.y = ggplot2::element_blank(),
                 axis.text.y = ggplot2::element_blank(),
                 axis.ticks.y = ggplot2::element_blank(),
                 legend.text = ggplot2::element_text(hjust = 0),
                 legend.key=element_rect(fill=NA),
                 panel.background=element_blank(),
                 panel.grid.major = ggplot2::element_blank(),
                 panel.grid.minor = ggplot2::element_blank(),
                 ...)
}

# APC model
smoothAPC <- function(data, admin.level, Amat = NULL, 
                      age.group = c("0", "1-11", "12-23", "24-35", "36-47", "48-59"), 
                      year.label, stratified = FALSE, 
                      mod = c('apc', 'ac', 'ap', 'pc', 'a', 'p', 'c')[1], slope.drop = NULL, 
                      family = c('betabinomial', 'binomial')[1],
                      pc.u = 1, pc.alpha = 0.01, pc.u.phi = 0.5, pc.alpha.phi = 2/3,
                      overdisp.mean = 0, overdisp.prec = 0.4, 
                      type.st = 4, pc.st.u = NA, pc.st.alpha = NA,
                      control.inla = list(strategy = 'adaptive', int.strategy = 'auto'),
                      inla.mode = c('classic', 'twostage', 'experimental')[3],
                      control.compute = list(config = TRUE), verbose = FALSE, ...){
  
  # inputs (for testing) ----
  
  # data = births.2014
  # admin.level = 'Admin1'
  # Amat = admin1.mat
  # age.groups = c("0", "1-11", "12-23", "24-35", "36-47", "48-59")
  # year.label = c(beg.year:end.proj.year)
  # mod = 'apc'
  # slope.drop = 'c'
  # stratified = TRUE
  # family = c('betabinomial', 'binomial')[1]
  # pc.u = 1
  # pc.alpha = 0.01
  # pc.u.phi = 0.5
  # pc.alpha.phi = 2/3
  # overdisp.mean = -7.5
  # overdisp.prec = 0.39
  # type.st = 4
  # pc.st.u = NA
  # pc.st.alpha = NA
  # control.inla = list(strategy = 'adaptive', int.strategy = 'auto')
  # inla.mode = c('classic', 'twostage', 'experimental')[3]
  # control.compute = list(config = TRUE)
  # verbose = FALSE
  
  # check inputs ----
  
  if(!(admin.level %in% c('National','Admin1','Admin2'))){
    stop('Enter a valid admin level (National, Admin1, or Admin2)')
  }
  if((admin.level!='National') == (is.null(Amat) == T)){
    stop('Admin level and adjacency matrix are not compatible')}
  if(!"config" %in% names(control.compute)){
    message("config = TRUE is added to control.compute so that posterior draws can be taken.")
    control.compute$config <- TRUE
  }
  if(mod == 'apc' && is.null(slope.drop)){
    stop('Warning: When fitting an APC model need to drop linear column of one of age, period or cohort.')
  }
  if(mod == 'apc' &&! (slope.drop %in% c('a', 'p', 'c'))){
    stop('slope.drop needs to be one of a, p or c')
  }
  if (!is.null(Amat)) {
    if (is.null(rownames(Amat))) {
      stop('Row names of Amat needs to be specified to region names.')
    }
    if (is.null(colnames(Amat))) {
      stop('Column names of Amat needs to be specified to region names.')
    }
    if (sum(rownames(Amat) != colnames(Amat)) > 0) {
      stop('Row and column names of Amat needs to be the same.')
    }
    is.spatial <- TRUE
  } else {
    Amat <- matrix(1, 1, 1)
    colnames(Amat) <- rownames(Amat) <- "All"
    is.spatial <- FALSE
  }
  
  # set admin-level parameters ----
  
  if(admin.level == 'National'){
    data$region <- "All"
  }else if(admin.level == 'Admin1'){
    data$region <- data$admin1.char
  }else if(admin.level == 'Admin2'){
    data$region <- data$admin2.char
  }
  
  # data oragnisation ----
  
  ## correct temporal specification ----
  
  # age + period + cohort
  data2a <- 
    expand.grid(ageMonth = 0:59,
                periodMonth = (12*(year.label[1] - 1900) + 1):(12*(year.label[length(year.label)] - 1900) + 12)) %>% 
    dplyr::mutate(cohortMonth = periodMonth - ageMonth,
                  age = dplyr::case_when(ageMonth == 0 ~ 0,
                                         ageMonth %in% 1:11 ~ 6,
                                         ageMonth %in% 12:23 ~ 17.5,
                                         ageMonth %in% 24:35 ~ 29.5,
                                         ageMonth %in% 36:47 ~ 41.5,
                                         TRUE ~ 53.5),
                  period = floor((periodMonth-1)/12)+1900,
                  cohort= floor((cohortMonth-1)/12)+1900,
                  age_id = age,
                  period_id = period %>% factor() %>% as.numeric(),
                  cohort_id = cohort %>% factor() %>% as.numeric()) %>% 
    dplyr::select(-tidyr::ends_with('Month')) %>% 
    dplyr::distinct()
  A <- data2a$age %>% unique() %>% length()
  P<- data2a$period %>% unique() %>% length()
  C <- data2a$cohort %>% unique() %>% length()
  
  ## correct spatial specification ----
  
  R <- nrow(Amat)
  data2b <- data.frame(space = rownames(Amat), space_id = 1:R)
  
  ## correct spatio-temporal specification ----
  
  ### space + age + period + cohort
  data2c <- 
    dplyr::cross_join(data2a, data2b) %>% 
    dplyr::mutate(spaceTime = interaction(space, period),
                  spaceTime_id = interaction(space_id, period_id) %>% as.numeric()) %>% 
    dplyr::relocate(c('space', 'spaceTime'), .before = age_id)
  
  ## finalised data ----
  
  data3 <- 
    data  %>%
    # labels for INLA modelling
    dplyr::mutate(
      # new terms
      age = dplyr::case_when(age == '0' ~ 0,
                             age == '1-11' ~ 6,
                             age == '12-23' ~ 17.5,
                             age == '24-35' ~ 29.5,
                             age == '36-47' ~ 41.5,
                             TRUE ~ 53.5),
      period = years,
      space = region,
      strata = urban) %>% 
    # join with the full range of space and time
    dplyr::left_join(., data2c, by = c('age', 'period', 'cohort', 'space')) %>% 
    dplyr::mutate(
      # final id terms for inla
      cluster_id = cluster %>% as.factor() %>% as.numeric(),
      strata_id = strata %>% factor(., levels = c('rural', 'urban')), # want this to be a factor in the model as we stratify by it
      age2_id = age_id, 
      period2_id = period_id, 
      cohort2_id = cohort_id)
  
  # priors ----
  
  # pc hyper priors
  ## temproal
  hyper_pc_time <- list(prec = list(prior = 'pc.prec', param = c(pc.u, pc.alpha)))
  ## spatial
  hyper_pc_space <- list(prec = list(prior = 'pc.prec', param = c(pc.u, pc.alpha)), phi = list(prior = 'pc', param = c(pc.u.phi, pc.alpha.phi)))
  ## spatio-temproal
  pc.st.u <- ifelse(is.na(pc.st.u), pc.u, pc.st.u)
  pc.st.alpha <- ifelse(is.na(pc.st.alpha), pc.alpha, pc.st.alpha)
  hyper_pc_spaceTime <- list(prec = list(prior = "pc.prec", param = c(pc.st.u, pc.st.alpha)))
  
  # formula ----
  
  ## temporal (APC) ----
  
  ### temporal constraints ----
  
  # age
  age.values <- data3$age_id %>% unique() %>% sort()
  age.rank.def <- 2
  age.cMat <- rbind(1, age.values) 
  age.e <- rep(0, times = age.rank.def)
  # period
  per.values <- data3$period_id %>% unique() %>% sort()
  per.rank.def <- 2
  per.cMat <- rbind(1, per.values) 
  per.e <- rep(0, times = per.rank.def)
  # cohort
  coh.values <- data3$cohort_id %>% unique() %>% sort()
  coh.rank.def <- 2
  coh.cMat <- rbind(1, coh.values) 
  coh.e <- rep(0, times = coh.rank.def)
  
  ### APC selection
  
  # full temporal formula
  formula <- 
    Y ~ age_id + period_id + cohort_id +
    f(age2_id, model = 'rw2', hyper = hyper_pc_time, scale.model = TRUE, constr = FALSE, rankdef = age.rank.def, extraconstr = list(A = age.cMat, e = age.e, values = age.values)) +
    f(period2_id, model = 'rw2', hyper = hyper_pc_time, scale.model = TRUE, constr = FALSE, rankdef = per.rank.def, extraconstr = list(A = per.cMat, e = per.e, values = per.values)) +
    f(cohort2_id, model = 'rw2', hyper = hyper_pc_time, scale.model = TRUE, constr = FALSE, rankdef = coh.rank.def, extraconstr = list(A = coh.cMat, e = coh.e, values = coh.values))
  
  message("----------------------------------",
          "\nModel Specification", 
          "\n Temporal effect(s):               ",
          "\n  APC type:                    ",
          mod,
          appendLF = FALSE)
  # update pre and post forumla for the temporal models
  # # do not remove any linear terms from pre fit as they are all needed later in 
  # # re-parameterisation
  if (mod == 'apc'){
    message("\n  Slope dropped:                 ", slope.drop, appendLF = FALSE)
    if(slope.drop == 'a'){
      formula <- update(formula, ~. - age_id)
    } else if (slope.drop == 'p'){
      formula <- update(formula, ~. - period_id)
    } else if (slope.drop == 'c'){
      formula <- update(formula, ~. - cohort_id)
    }
  } else if (mod == 'ac'){
    formula <- update(formula, ~.
                      - period_id 
                      - f(period2_id, model = 'rw2', hyper = hyper_pc_time, scale.model = TRUE, constr = FALSE, rankdef = per.rank.def, extraconstr = list(A = per.cMat, e = per.e, values = per.values)))
  } else if (mod == 'pc'){
    formula <- update(formula, ~.
                      - age_id 
                      - f(age2_id, model = 'rw2', hyper = hyper_pc_time, scale.model = TRUE, constr = FALSE, rankdef = age.rank.def, extraconstr = list(A = age.cMat, e = age.e, values = age.values)))
  } else if (mod == 'ap'){
    formula <- update(formula, ~.
                      - cohort_id 
                      -  f(cohort2_id, model = 'rw2', hyper = hyper_pc_time, scale.model = TRUE, constr = FALSE, rankdef = coh.rank.def, extraconstr = list(A = coh.cMat, e = coh.e, values = coh.values)))
  }
  
  ## spatial & spatio-temporal ----
  
  if (is.spatial) {
    message("\n Spatial effect(s):               ",
            "\n  Structured:                 bym2",
            "\n  Interaction type:              ", type.st,
            appendLF = FALSE)
    
    ### spatial ----
    
    formula <- 
      update(formula, ~. +
               f(space_id, graph = Amat, model = 'bym2', hyper= hyper_pc_space, scale.model = TRUE, adjust.for.con.comp = TRUE))
    
    
    ### spatio-temporal ----
    
    inla.rw <- utils::getFromNamespace('inla.rw', 'INLA')
    inla.scale.model <- utils::getFromNamespace('inla.scale.model', 'INLA')
    
    # precision matrices
    ## time
    ### unstructured
    Q_timeUnstruc <- diag(P)
    ### structured
    Q_timeStruc <- inla.rw(n = P, order = 2, scale.model = TRUE, sparse = TRUE)
    ## space
    ### unstructured
    Q_spaceUnstruc <- diag(R)
    ### structured
    Q_spaceStruc <- Amat
    #### organise Amat into ICAR precision
    if (sum(Q_spaceStruc > 0 & Q_spaceStruc < 1) != 0) {
      #### if 0 < Q_{ij} < 1 then make Q_{ij} = 1
      for (i in 1:nrow(Q_spaceStruc)) {
        idx <- which(Q_spaceStruc[i, ] > 0 & Q_spaceStruc[i, ] < 1)
        Q_spaceStruc[i, idx] <- 1
      }
    }
    #### turn Amat into an ICAR precision
    ##### no 0s on diagonal; main diag to be the number of neighbours in total; off diagonal to be 0 for not neighbours and -1 for neighbours
    diag(Q_spaceStruc) <- 0
    diag <- apply(Q_spaceStruc, 1, sum)
    Q_spaceStruc[Q_spaceStruc != 0] <- -1
    diag(Q_spaceStruc) <- diag
    Q_spaceStruc <- inla.scale.model(Q_spaceStruc, list(A = matrix(1, 1, dim(Q_spaceStruc)[1]), e = 0))
    
    if (type.st == 1) {
      
      # kronecker product for the Type I Interaction
      ## IID x IID
      formula <- update(formula, ~. + f(spaceTime_id, model="iid", hyper = hyper_pc_spaceTime))
      
    } else {
      
      if (type.st == 2) {
        # kronecker product for the Type II Interaction
        ## RW2 x IID
        Q_spaceTime <- Q_timeStruc %x% Q_spaceUnstruc
      } else if (type.st == 3) {
        # kronecker product for the Type III Interaction
        ## IID x ICAR
        Q_spaceTime <- Q_timeUnstruc %x% Q_spaceStruc
      } else {
        # kronecker product for the Type IV Interaction
        ## RW2 x ICAR
        Q_spaceTime <- Q_timeStruc %x% Q_spaceStruc
      }
      
      # constraints for space-time term
      ## sum to zero constraints for identifiability
      ## rank def is due to 0 eigenvalues
      ## define constraint matrix from the eigenvectors relating to zero eigenvalues
      ## rank def is the number of 0 eigenvalues
      eigenQ <- eigen(Q_spaceTime)
      ids <- which(eigenQ$values < 1e-10)
      cMat <- t(eigenQ$vectors[,ids])
      
      # # check the dimensions line up correctly
      # R*(P-2) + nrow(cMat) == nrow(Q_spaceTime) # Type II
      # (R-1)*P + nrow(cMat) == nrow(Q_spaceTime) # Type III
      # (R-1)*(P-2) + nrow(cMat) == nrow(Q_spaceTime) # Type IV
      
      # add spatio-temporal term to formula
      formula <- update(formula, ~. +
                          f(spaceTime_id,
                            model = "generic0",
                            Cmatrix = Q_spaceTime,
                            extraconstr = list(A = cMat, e = rep(0, nrow(cMat))),
                            rankdef = nrow(cMat),
                            hyper = hyper_pc_spaceTime))
      
    }
    
  }
  
  ## cluster (if binomial) ----
  
  # cluster random effect for binomial model
  if (family == 'binomial') {
    message("\n  cluster:                     iid", appendLF = FALSE)
    formula <- update(formula, ~. +  f(cluster_id, model = 'iid', hyper = hyper_pc_time))
  }
  
  ## stratification ----
  
  # sorting out the fixed effects
  ## urban/rural strata
  if(stratified){
    message("\n  Stratified:                  yes", appendLF = FALSE)
    message("\n  Global Intercept:             no", appendLF = FALSE)
    # update the formula
    formula <- update(formula, ~. - 1 + strata_id)
  } else {
    message("\n  Stratified:                   no", appendLF = FALSE)
    message("\n  Global intercept:            yes", appendLF = FALSE)
    data3$strata <- as.factor('All')
    data3$strata_id <- data3$strata_id <- 1
  }
  
  # model fit ----
  
  message("\nFitting model...",
          appendLF = FALSE)
  startClock <- proc.time() # Start the clock!
  # the extra information depending on the model
  if(family == 'betabinomial'){
    control.family <- list(hyper = list(rho = list(param = c(overdisp.mean, overdisp.prec), initial = overdisp.mean)))
    fit <-
      INLA::inla(formula, 
                 family = family, 
                 data = data3, 
                 Ntrials = data3$total,
                 control.compute = control.compute,
                 control.family = control.family,
                 control.predictor = list(compute = FALSE, link = 1),
                 control.inla = control.inla,
                 lincomb = NULL,
                 inla.mode = inla.mode,
                 verbose = verbose,
                 ...)
  } else {
    fit <-
      INLA::inla(formula, 
                 family = family, 
                 data = data3, 
                 Ntrials = data3$total,
                 control.compute = control.compute,
                 control.predictor = list(compute = FALSE, link = 1),
                 control.inla = control.inla,
                 lincomb = NULL,
                 inla.mode = inla.mode,
                 verbose = verbose,
                 ...)
  }
  endClock <- proc.time() - startClock # Stop the clock
  message("\n Time take to fit model(s):     ", endClock[3] %>% as.numeric %>% round(),
          "\n----------------------------------",
          appendLF = FALSE)
  
  # output ----
  
  priors <- list(pc.u = pc.u, pc.alpha = pc.alpha, pc.u.phi = pc.u.phi, 
                 pc.alpha.phi = pc.alpha.phi, 
                 pc.st.u = pc.st.u, 
                 pc.st.alpha = pc.st.alpha, 
                 overdisp.mean = overdisp.mean, overdisp.prec = overdisp.prec)
  
  return(list(model = formula, fit = fit, family = family, 
              Amat = Amat, newdata = data3, type.st = type.st,
              age.group = age.group, year.label = year.label,
              priors = priors))
}

smoothAPC_noInteraction <- function(data, admin.level, Amat = NULL, 
                      age.group = c("0", "1-11", "12-23", "24-35", "36-47", "48-59"), 
                      year.label, stratified = FALSE, 
                      mod = c('apc', 'ac', 'ap', 'pc', 'a', 'p', 'c')[1], slope.drop = NULL, 
                      family = c('betabinomial', 'binomial')[1],
                      pc.u = 1, pc.alpha = 0.01, pc.u.phi = 0.5, pc.alpha.phi = 2/3,
                      overdisp.mean = 0, overdisp.prec = 0.4, 
                      type.st = 4, pc.st.u = NA, pc.st.alpha = NA,
                      control.inla = list(strategy = 'adaptive', int.strategy = 'auto'),
                      inla.mode = c('classic', 'twostage', 'experimental')[3],
                      control.compute = list(config = TRUE), verbose = FALSE, ...){
  
  # inputs (for testing) ----
  
  # data = births.2014
  # admin.level = 'Admin1'
  # Amat = admin1.mat
  # age.groups = c("0", "1-11", "12-23", "24-35", "36-47", "48-59")
  # year.label = c(beg.year:end.proj.year)
  # mod = 'apc'
  # slope.drop = 'c'
  # stratified = TRUE
  # family = c('betabinomial', 'binomial')[1]
  # pc.u = 1
  # pc.alpha = 0.01
  # pc.u.phi = 0.5
  # pc.alpha.phi = 2/3
  # overdisp.mean = -7.5
  # overdisp.prec = 0.39
  # type.st = 4
  # pc.st.u = NA
  # pc.st.alpha = NA
  # control.inla = list(strategy = 'adaptive', int.strategy = 'auto')
  # inla.mode = c('classic', 'twostage', 'experimental')[3]
  # control.compute = list(config = TRUE)
  # verbose = FALSE
  
  # check inputs ----
  
  if(!(admin.level %in% c('National','Admin1','Admin2'))){
    stop('Enter a valid admin level (National, Admin1, or Admin2)')
  }
  if((admin.level!='National') == (is.null(Amat) == T)){
    stop('Admin level and adjacency matrix are not compatible')}
  if(!"config" %in% names(control.compute)){
    message("config = TRUE is added to control.compute so that posterior draws can be taken.")
    control.compute$config <- TRUE
  }
  if(mod == 'apc' && is.null(slope.drop)){
    stop('Warning: When fitting an APC model need to drop linear column of one of age, period or cohort.')
  }
  if(mod == 'apc' &&! (slope.drop %in% c('a', 'p', 'c'))){
    stop('slope.drop needs to be one of a, p or c')
  }
  if (!is.null(Amat)) {
    if (is.null(rownames(Amat))) {
      stop('Row names of Amat needs to be specified to region names.')
    }
    if (is.null(colnames(Amat))) {
      stop('Column names of Amat needs to be specified to region names.')
    }
    if (sum(rownames(Amat) != colnames(Amat)) > 0) {
      stop('Row and column names of Amat needs to be the same.')
    }
    is.spatial <- TRUE
  } else {
    Amat <- matrix(1, 1, 1)
    colnames(Amat) <- rownames(Amat) <- "All"
    is.spatial <- FALSE
  }
  
  # set admin-level parameters ----
  
  if(admin.level == 'National'){
    data$region <- "All"
  }else if(admin.level == 'Admin1'){
    data$region <- data$admin1.char
  }else if(admin.level == 'Admin2'){
    data$region <- data$admin2.char
  }
  
  # data oragnisation ----
  
  ## correct temporal specification ----
  
  # age + period + cohort
  data2a <- 
    expand.grid(ageMonth = 0:59,
                periodMonth = (12*(year.label[1] - 1900) + 1):(12*(year.label[length(year.label)] - 1900) + 12)) %>% 
    dplyr::mutate(cohortMonth = periodMonth - ageMonth,
                  age = dplyr::case_when(ageMonth == 0 ~ 0,
                                         ageMonth %in% 1:11 ~ 6,
                                         ageMonth %in% 12:23 ~ 17.5,
                                         ageMonth %in% 24:35 ~ 29.5,
                                         ageMonth %in% 36:47 ~ 41.5,
                                         TRUE ~ 53.5),
                  period = floor((periodMonth-1)/12)+1900,
                  cohort= floor((cohortMonth-1)/12)+1900,
                  age_id = age,
                  period_id = period %>% factor() %>% as.numeric(),
                  cohort_id = cohort %>% factor() %>% as.numeric()) %>% 
    dplyr::select(-tidyr::ends_with('Month')) %>% 
    dplyr::distinct()
  A <- data2a$age %>% unique() %>% length()
  P<- data2a$period %>% unique() %>% length()
  C <- data2a$cohort %>% unique() %>% length()
  
  ## correct spatial specification ----
  
  R <- nrow(Amat)
  data2b <- data.frame(space = rownames(Amat), space_id = 1:R)
  
  ## correct spatio-temporal specification ----
  
  ### space + age + period + cohort
  data2c <- 
    dplyr::cross_join(data2a, data2b) %>% 
    dplyr::mutate(spaceTime = interaction(space, period),
                  spaceTime_id = interaction(space_id, period_id) %>% as.numeric()) %>% 
    dplyr::relocate(c('space', 'spaceTime'), .before = age_id)
  
  ## finalised data ----
  
  data3 <- 
    data  %>%
    # labels for INLA modelling
    dplyr::mutate(
      # new terms
      age = dplyr::case_when(age == '0' ~ 0,
                             age == '1-11' ~ 6,
                             age == '12-23' ~ 17.5,
                             age == '24-35' ~ 29.5,
                             age == '36-47' ~ 41.5,
                             TRUE ~ 53.5),
      period = years,
      space = region,
      strata = urban) %>% 
    # join with the full range of space and time
    dplyr::left_join(., data2c, by = c('age', 'period', 'cohort', 'space')) %>% 
    dplyr::mutate(
      # final id terms for inla
      cluster_id = cluster %>% as.factor() %>% as.numeric(),
      strata_id = strata %>% factor(., levels = c('rural', 'urban')), # want this to be a factor in the model as we stratify by it
      age2_id = age_id, 
      period2_id = period_id, 
      cohort2_id = cohort_id)
  
  # priors ----
  
  # pc hyper priors
  ## temproal
  hyper_pc_time <- list(prec = list(prior = 'pc.prec', param = c(pc.u, pc.alpha)))
  ## spatial
  hyper_pc_space <- list(prec = list(prior = 'pc.prec', param = c(pc.u, pc.alpha)), phi = list(prior = 'pc', param = c(pc.u.phi, pc.alpha.phi)))
  ## spatio-temproal
  pc.st.u <- ifelse(is.na(pc.st.u), pc.u, pc.st.u)
  pc.st.alpha <- ifelse(is.na(pc.st.alpha), pc.alpha, pc.st.alpha)
  hyper_pc_spaceTime <- list(prec = list(prior = "pc.prec", param = c(pc.st.u, pc.st.alpha)))
  
  # formula ----
  
  ## temporal (APC) ----
  
  ### temporal constraints ----
  
  # age
  age.values <- data3$age_id %>% unique() %>% sort()
  age.rank.def <- 2
  age.cMat <- rbind(1, age.values) 
  age.e <- rep(0, times = age.rank.def)
  # period
  per.values <- data3$period_id %>% unique() %>% sort()
  per.rank.def <- 2
  per.cMat <- rbind(1, per.values) 
  per.e <- rep(0, times = per.rank.def)
  # cohort
  coh.values <- data3$cohort_id %>% unique() %>% sort()
  coh.rank.def <- 2
  coh.cMat <- rbind(1, coh.values) 
  coh.e <- rep(0, times = coh.rank.def)
  
  ### APC selection
  
  # full temporal formula
  formula <- 
    Y ~ age_id + period_id + cohort_id +
    f(age2_id, model = 'rw2', hyper = hyper_pc_time, scale.model = TRUE, constr = FALSE, rankdef = age.rank.def, extraconstr = list(A = age.cMat, e = age.e, values = age.values)) +
    f(period2_id, model = 'rw2', hyper = hyper_pc_time, scale.model = TRUE, constr = FALSE, rankdef = per.rank.def, extraconstr = list(A = per.cMat, e = per.e, values = per.values)) +
    f(cohort2_id, model = 'rw2', hyper = hyper_pc_time, scale.model = TRUE, constr = FALSE, rankdef = coh.rank.def, extraconstr = list(A = coh.cMat, e = coh.e, values = coh.values))
  
  message("----------------------------------",
          "\nModel Specification", 
          "\n Temporal effect(s):               ",
          "\n  APC type:                    ",
          mod,
          appendLF = FALSE)
  # update pre and post forumla for the temporal models
  # # do not remove any linear terms from pre fit as they are all needed later in 
  # # re-parameterisation
  if (mod == 'apc'){
    message("\n  Slope dropped:                 ", slope.drop, appendLF = FALSE)
    if(slope.drop == 'a'){
      formula <- update(formula, ~. - age_id)
    } else if (slope.drop == 'p'){
      formula <- update(formula, ~. - period_id)
    } else if (slope.drop == 'c'){
      formula <- update(formula, ~. - cohort_id)
    }
  } else if (mod == 'ac'){
    formula <- update(formula, ~.
                      - period_id 
                      - f(period2_id, model = 'rw2', hyper = hyper_pc_time, scale.model = TRUE, constr = FALSE, rankdef = per.rank.def, extraconstr = list(A = per.cMat, e = per.e, values = per.values)))
  } else if (mod == 'pc'){
    formula <- update(formula, ~.
                      - age_id 
                      - f(age2_id, model = 'rw2', hyper = hyper_pc_time, scale.model = TRUE, constr = FALSE, rankdef = age.rank.def, extraconstr = list(A = age.cMat, e = age.e, values = age.values)))
  } else if (mod == 'ap'){
    formula <- update(formula, ~.
                      - cohort_id 
                      -  f(cohort2_id, model = 'rw2', hyper = hyper_pc_time, scale.model = TRUE, constr = FALSE, rankdef = coh.rank.def, extraconstr = list(A = coh.cMat, e = coh.e, values = coh.values)))
  }
  
  ## spatial & spatio-temporal ----
  
  if (is.spatial) {
    message("\n Spatial effect(s):               ",
            "\n  Structured:                 bym2",
            # "\n  Interaction type:              ", type.st,
            appendLF = FALSE)
    
    ### spatial ----
    
    formula <- 
      update(formula, ~. +
               f(space_id, graph = Amat, model = 'bym2', hyper= hyper_pc_space, scale.model = TRUE, adjust.for.con.comp = TRUE))
    
    
    # ### spatio-temporal ----
    # 
    # inla.rw <- utils::getFromNamespace('inla.rw', 'INLA')
    # inla.scale.model <- utils::getFromNamespace('inla.scale.model', 'INLA')
    # 
    # # precision matrices
    # ## time
    # ### unstructured
    # Q_timeUnstruc <- diag(P)
    # ### structured
    # Q_timeStruc <- inla.rw(n = P, order = 2, scale.model = TRUE, sparse = TRUE)
    # ## space
    # ### unstructured
    # Q_spaceUnstruc <- diag(R)
    # ### structured
    # Q_spaceStruc <- Amat
    # #### organise Amat into ICAR precision
    # if (sum(Q_spaceStruc > 0 & Q_spaceStruc < 1) != 0) {
    #   #### if 0 < Q_{ij} < 1 then make Q_{ij} = 1
    #   for (i in 1:nrow(Q_spaceStruc)) {
    #     idx <- which(Q_spaceStruc[i, ] > 0 & Q_spaceStruc[i, ] < 1)
    #     Q_spaceStruc[i, idx] <- 1
    #   }
    # }
    # #### turn Amat into an ICAR precision
    # ##### no 0s on diagonal; main diag to be the number of neighbours in total; off diagonal to be 0 for not neighbours and -1 for neighbours
    # diag(Q_spaceStruc) <- 0
    # diag <- apply(Q_spaceStruc, 1, sum)
    # Q_spaceStruc[Q_spaceStruc != 0] <- -1
    # diag(Q_spaceStruc) <- diag
    # Q_spaceStruc <- inla.scale.model(Q_spaceStruc, list(A = matrix(1, 1, dim(Q_spaceStruc)[1]), e = 0))
    # 
    # if (type.st == 1) {
    #   
    #   # kronecker product for the Type I Interaction
    #   ## IID x IID
    #   formula <- update(formula, ~. + f(spaceTime_id, model="iid", hyper = hyper_pc_spaceTime))
    #   
    # } else {
    #   
    #   if (type.st == 2) {
    #     # kronecker product for the Type II Interaction
    #     ## RW2 x IID
    #     Q_spaceTime <- Q_timeStruc %x% Q_spaceUnstruc
    #   } else if (type.st == 3) {
    #     # kronecker product for the Type III Interaction
    #     ## IID x ICAR
    #     Q_spaceTime <- Q_timeUnstruc %x% Q_spaceStruc
    #   } else {
    #     # kronecker product for the Type IV Interaction
    #     ## RW2 x ICAR
    #     Q_spaceTime <- Q_timeStruc %x% Q_spaceStruc
    #   }
    #   
    #   # constraints for space-time term
    #   ## sum to zero constraints for identifiability
    #   ## rank def is due to 0 eigenvalues
    #   ## define constraint matrix from the eigenvectors relating to zero eigenvalues
    #   ## rank def is the number of 0 eigenvalues
    #   eigenQ <- eigen(Q_spaceTime)
    #   ids <- which(eigenQ$values < 1e-10)
    #   cMat <- t(eigenQ$vectors[,ids])
    #   
    #   # # check the dimensions line up correctly
    #   # R*(P-2) + nrow(cMat) == nrow(Q_spaceTime) # Type II
    #   # (R-1)*P + nrow(cMat) == nrow(Q_spaceTime) # Type III
    #   # (R-1)*(P-2) + nrow(cMat) == nrow(Q_spaceTime) # Type IV
    #   
    #   # add spatio-temporal term to formula
    #   formula <- update(formula, ~. +
    #                       f(spaceTime_id,
    #                         model = "generic0",
    #                         Cmatrix = Q_spaceTime,
    #                         extraconstr = list(A = cMat, e = rep(0, nrow(cMat))),
    #                         rankdef = nrow(cMat),
    #                         hyper = hyper_pc_spaceTime))
    #   
    # }
    
  }
  
  ## cluster (if binomial) ----
  
  # cluster random effect for binomial model
  if (family == 'binomial') {
    message("\n  cluster:                     iid", appendLF = FALSE)
    formula <- update(formula, ~. +  f(cluster_id, model = 'iid', hyper = hyper_pc_time))
  }
  
  ## stratification ----
  
  # sorting out the fixed effects
  ## urban/rural strata
  if(stratified){
    message("\n  Stratified:                  yes", appendLF = FALSE)
    message("\n  Global Intercept:             no", appendLF = FALSE)
    # update the formula
    formula <- update(formula, ~. - 1 + strata_id)
  } else {
    message("\n  Stratified:                   no", appendLF = FALSE)
    message("\n  Global intercept:            yes", appendLF = FALSE)
    data3$strata <- as.factor('All')
    data3$strata_id <- data3$strata_id <- 1
  }
  
  # model fit ----
  
  message("\nFitting model...",
          appendLF = FALSE)
  startClock <- proc.time() # Start the clock!
  # the extra information depending on the model
  if(family == 'betabinomial'){
    control.family <- list(hyper = list(rho = list(param = c(overdisp.mean, overdisp.prec), initial = overdisp.mean)))
    fit <-
      INLA::inla(formula, 
                 family = family, 
                 data = data3, 
                 Ntrials = data3$total,
                 control.compute = control.compute,
                 control.family = control.family,
                 control.predictor = list(compute = FALSE, link = 1),
                 control.inla = control.inla,
                 lincomb = NULL,
                 inla.mode = inla.mode,
                 verbose = verbose,
                 ...)
  } else {
    fit <-
      INLA::inla(formula, 
                 family = family, 
                 data = data3, 
                 Ntrials = data3$total,
                 control.compute = control.compute,
                 control.predictor = list(compute = FALSE, link = 1),
                 control.inla = control.inla,
                 lincomb = NULL,
                 inla.mode = inla.mode,
                 verbose = verbose,
                 ...)
  }
  endClock <- proc.time() - startClock # Stop the clock
  message("\n Time take to fit model(s):     ", endClock[3] %>% as.numeric %>% round(),
          "\n----------------------------------",
          appendLF = FALSE)
  
  # output ----
  
  priors <- list(pc.u = pc.u, pc.alpha = pc.alpha, pc.u.phi = pc.u.phi, 
                 pc.alpha.phi = pc.alpha.phi, 
                 pc.st.u = pc.st.u, 
                 pc.st.alpha = pc.st.alpha, 
                 overdisp.mean = overdisp.mean, overdisp.prec = overdisp.prec)
  
  return(list(model = formula, fit = fit, family = family, 
              Amat = Amat, newdata = data3, type.st = type.st,
              age.group = age.group, year.label = year.label,
              priors = priors))
}

