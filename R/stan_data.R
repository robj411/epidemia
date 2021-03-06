# Parses arguments into complete list of objects ready to pass as data
# arguments into stanmodels::epidemia_base
#
# @inheritParams epim
# @param group, x Objects returned by parse_mm
# @param link Not yet used.
standata_all <- function(rt,
                         obs,
                         data,
                         pops,
                         si,
                         seed_days,
                         group_subset,
                         prior_tau,
                         prior_PD) {

  # standata for general model params
  out <- standata_data(data)
  out <- c(
    out,
    list(
      si = pad(si, out$NS, 0, TRUE),
      N0 = seed_days,
      prior_PD = prior_PD,
      pop = as.array(pops$pop)
    ),
    standata_model_priors(prior_tau),
    standata_obs(
      obs = obs,
      groups = out$groups,
      nsim = out$NS,
      begin = out$begin
    ),
    standata_rt(rt)
  )
  
  ## crop phylo. NB: vectors are going backwards in time.
  if(out$N_phy > 0){
    tocrop <- c("ltt_terms_phy","nco_phy","intervalLength_phy","ne_bin_phy")
    # the tree data start at the earliest on the same date as the rest of the data and tree_delay=0. Otherwise, there is a delay > 1
    tree_delay <- out$tree_delay
    maxbin <- max(out$ne_bin_phy)
    N2 <- out$N2
    # ne_bin_phy is the position of the I vector. Position 1 is the first date for which there is data. 
    out$ne_bin_phy <- maxbin + tree_delay - out$ne_bin_phy + 1
    for(tc in tocrop) out[[tc]] <- out[[tc]][out$ne_bin_phy<=N2]
    out$N_phy <- length(out$ne_bin_phy)
  }
  
  return(out)
}

# Generate standata for autocorrelation terms
#
# @inheritParams epim
# @param formula, data Same as in epim
standata_autocor <- function(object) {
  out <- list()
  
  if (!inherits(object, "epirt_"))
    stop("Bug found. 'object' must have class 'epirt_'.")
  
  formula <- formula(object)
  
  if (is_autocor(formula)) {
    autocor <- object$autocor
    out$ac_nterms <- length(autocor$nproc)
    out$ac_ntime <- as.array(autocor$ntime)
    out$ac_q <- sum(autocor$ntime)
    out$ac_nproc <- sum(autocor$nproc)
    # todo: implement this as an option
    out$ac_prior_scales <- as.array(autocor$prior_scale)
    
    # add sparse matrix representation
    parts <- rstan::extract_sparse_parts(autocor$Z)
    out$ac_v <- parts$v - 1L
    out$ac_nnz <- length(out$ac_v)
  } else {
    out$ac_nterms <- out$ac_q <- out$ac_nproc <- out$ac_nnz <- 0
    out$ac_prior_scales <- out$ac_v <- out$ac_ntime <- numeric()
  }
  return(out)
}

# Generate standata for phylo terms
#
# @inheritParams epim
# @param formula, data Same as in epim
standata_phylo <- function(object) {
  out <- list()
  
  if (!inherits(object, "epirt_"))
    stop("Bug found. 'object' must have class 'epirt_'.")
  #formula <- formula(object)
  if ('phylo'%in%names(object)) {
    phylo <- object$phylo
    for(nm in names(phylo))
      out[[nm]] <- phylo[[nm]]
  } else {
    out$N_phy <- 0
    out$intervalLength_phy <- out$ne_bin_phy <- out$nco_phy <- out$ltt_terms_phy <- numeric()
  }
  return(out)
}

# Generate relevant standata from data. Used internally in epim.
#
# @param data The result of checkData
standata_data <- function(data) {
  groups <- sort(levels(data$group))
  M <- length(groups)
  NC <- as.numeric(table(data$group))
  max_sim <- max(NC)
  # compute first date
  starts <- aggregate(date ~ group, data = data, FUN = min)$date
  begin <- min(starts)
  # integer index of start (1 being 'begin')
  starts <- as.numeric(starts - begin + 1)
  return(list(
    groups = groups,
    M = M,
    NC = as.array(NC),
    NS = max_sim,
    N2 = max(starts + NC - 1),
    starts = as.array(starts),
    begin = begin
  ))
}

# Parses rt argument into data ready for stan
#
# @param rt An epirt_ object
# @param data data argument to epim
standata_rt <- function(rt) {
  out <- list()
  out$r0 <- rt$r0
  out <- c(
    out,
    standata_reg(rt)
  )
  out$rt_prior_info = out$prior_info
  return(out)
}

# Parses obs argument into data ready for stan
#
# @param obs A list of epiobs_ objects
# @param groups An ordered list of all modeled groups
# @param nsim The total simulation period
# @param begin The first simulation date
standata_obs <- function(obs, groups, nsim, begin) {
  maxtypes <- 10
  types <- length(obs)
  out <- list()

  if (types > maxtypes) {
    stop("Currently up to 10 observation types are
     supported. Please check 'obs'", .call = FALSE)
  }

  if (types > 0) {
    oN <- sapply(obs, function(x) nobs(x))
    oN <- array(pad(oN, maxtypes, 0))

    pvecs <- as.array(lapply(obs,
      function(x) pad_i2o(x, len = nsim)))

    dat <- array(unlist(
      lapply(
        obs,
        function(x) get_obs(x)
      )
    ))
    obs_group <- array(unlist(
      lapply(
        obs,
        function(x) match(get_gr(x), groups)
      )
    ))
    obs_date <- unlist(
      lapply(
        obs,
        function(x) as.character(get_time(x))
      )
    )
    obs_date <- array(as.Date(obs_date) - begin + 1)
    obs_type <- array(unlist(Map(rep, 1:maxtypes, oN)))

    # compute regression quantities
    reg <- lapply(obs, standata_reg)
    oK <- sapply(reg, function(x) x$K)
    oK <- array(pad(oK, maxtypes, 0))
    K_all <- sum(oK)

    # intercepts
    has_ointercept <- sapply(reg, function(x) x$has_intercept)
    num_ointercepts <- sum(has_ointercept)
    has_ointercept <- array(has_ointercept * cumsum(has_ointercept))

    # covariates
    oxbar <- array(unlist(lapply(reg, function(x) x$xbar)))
    nms <- paste("oX", 1:maxtypes, sep = "")
    for (i in seq_along(nms)) {
      if (i <= types) {
        out[[nms[i]]] <- reg[[i]]$X
      }
      else {
        out[[nms[i]]] <- array(0, dim = c(0, 0))
      }
    }

    # regression hyperparameters
    prior_omean <- array(unlist(
      lapply(
        reg,
        function(x) x$prior_mean
      )
    ))
    prior_oscale <- array(unlist(lapply(
      reg,
      function(x) x$prior_scale
    )))
    prior_mean_for_ointercept <- array(sapply(
      reg,
      function(x) x$prior_mean_for_intercept
    ))
    prior_scale_for_ointercept <- array(sapply(
      reg,
      function(x) x$prior_scale_for_intercept
    ))

    # auxiliary params
    ofamily <- array(sapply(reg, function(x) x$family))
    olink <- array(sapply(reg, function(x) x$link))
    has_oaux <- ofamily != 1
    num_oaux <- sum(has_oaux)
    has_oaux <- array(has_oaux * cumsum(has_oaux))

    nms_aux <- c(
      "prior_dist_for_oaux",
      "prior_mean_for_oaux",
      "prior_scale_for_oaux",
      "prior_df_for_oaux"
    )
    for (i in nms_aux){
      temp <- unlist(lapply(reg, function(x) x[[i]]))
      assign(i, array(as.numeric(temp) %ORifNULL% rep(0,0)))
    }

    has_offset <- array(sapply(obs, function(x) x$has_offset * 1))
    offset_ <- array(unlist(lapply(obs, function(x) x$offset)))

    obs_prior_info <- lapply(reg, function(x) x$prior_info)
  }
  else { # set to zero values
    N_obs <- K_all <- num_ointercepts <- num_oaux <- 0
    dat <- prior_omean <- prior_oscale <-
    prior_mean_for_ointercept <- prior_scale_for_ointercept <-
    prior_mean_for_oaux <- prior_scale_for_oaux <- 
    prior_df_for_oaux <- offset_ <- rep(0,0)
    obs_group <- obs_date <- obs_type  <- oxbar <-
    has_ointercept <- prior_dist_for_oaux <- has_offset <- 
    has_oaux <- olink <- ofamily <- integer(0)
    oN <- oK <- rep(0, maxtypes)
    pvecs <- array(0, dim = c(0, nsim))
    obs_prior_info <- NULL

     # covariates
    nms <- paste("oX", 1:maxtypes, sep = "")
    for (i in seq_along(nms)) {
        out[[nms[i]]] <- array(0, dim = c(0, 0))
    }
  }

  out <- c(out, loo::nlist(
    obs_prior_info,
    N_obs = sum(oN),
    R = types,
    oN,
    obs = dat,
    obs_group,
    obs_date,
    obs_type,
    oK,
    K_all,
    oxbar,
    has_ointercept,
    num_ointercepts,
    prior_omean,
    prior_oscale,
    prior_mean_for_ointercept,
    prior_scale_for_ointercept,
    ofamily,
    olink,
    has_oaux,
    num_oaux,
    prior_dist_for_oaux,
    prior_mean_for_oaux,
    prior_scale_for_oaux,
    prior_df_for_oaux,
    pvecs,
    has_offset,
    offset_
  ))
  return(out)
}

# Parse model priors (currently just tau) to standata representation.
# Used internally in epim.
#
# @param prior_tau See \code{\link{epim}}
standata_model_priors <- function(prior_tau) {

  # for passing R CMD Check
  prior_scale_for_tau <- NULL

  prior_tau_stuff <- handle_glm_prior(
    prior = prior_tau,
    nvars = 1,
    default_scale = 1 / 0.03,
    link = "dummy",
    ok_dists = loo::nlist("exponential")
  )

  names(prior_tau_stuff) <- paste0(names(prior_tau_stuff), "_for_tau")

  for (i in names(prior_tau_stuff)) {
    assign(i, prior_tau_stuff[[i]])
  }

  return(list(
    prior_scale_for_tau = as.numeric(prior_scale_for_tau)
  ))
}


# Takes a vector and extends to required length
#
# @param vec The vector to pad
# @param len The exact required length
# @param tol The value to impute for extended entries
# @param sv If TRUE, normalise for sum of 1
pad <- function(x, len, a, sv = FALSE) {
  out <- c(x, rep(a, max(len - length(x), 0)))
  out <- out[1:len]
  if (sv) {
    out <- out / sum(out)
  }
  return(out)
}