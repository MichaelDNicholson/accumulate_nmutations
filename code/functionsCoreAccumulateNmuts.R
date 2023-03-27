#functions for k type fluctuation analysis

#### Distribution functions ####

#Omega recur #
get_scale_tail_vec  = function(params){
  #input vec is parameter data.frame/tibble
  alpha_vec = params$brate
  beta_vec = params$drate
  nu_vec = params$mutrate
  lambda_vec = alpha_vec - beta_vec
  delta_vec = sapply(1:length(lambda_vec), function(k) max(lambda_vec[1:k]))
  rn_vec = sapply(1:length(lambda_vec), function(k) sum(lambda_vec[1:k] == delta_vec[k]))
  
  
  omega_vec = numeric(length(alpha_vec))
  omega_vec[1] = (1-beta_vec[1]/alpha_vec[1])^(-1)
  
  for (n in 1:(length(omega_vec)-1)){
    
    
    if (lambda_vec[n+1]>delta_vec[n]){ #fitness increase
      
      an = (lambda_vec[n+1]/alpha_vec[n+1])^(1-delta_vec[n]/lambda_vec[n+1]) / (lambda_vec[n+1]*delta_vec[n]^(rn_vec[n]-1)) *
        pi*sin(pi*delta_vec[n]/lambda_vec[n+1])
      
      omega_vec[n+1] = (omega_vec[n]*nu_vec[n]*log(nu_vec[n]^(-1))^(rn_vec[n]-1)*an)^(lambda_vec[n+1]/delta_vec[n])                                                          
      
    } else if (lambda_vec[n+1] == delta_vec[n]){ #neutral
      
      omega_vec[n+1] = nu_vec[n]/(rn_vec[n])*omega_vec[n]
      
    } else { #deleterious
      omega_vec[n+1] = nu_vec[n]/(delta_vec[n]-lambda_vec[n+1])*omega_vec[n]
    }
    
  }
  
  return(list(omega_vec, delta_vec))
  
}

get_mutantcounts_scaledby_tfn <- function(params,sizes){
  alpha_vec = params$brate
  beta_vec = params$drate
  lambda_vec = alpha_vec - beta_vec
  delta_vec = sapply(1:length(lambda_vec), function(k) max(lambda_vec[1:k]))
  rn_vec = sapply(1:length(lambda_vec), function(k) sum(lambda_vec[1:k] == delta_vec[k]))
  
  det_timefn <- exp(delta_vec*params$fintime)*(params$fintime)^(rn_vec-1)
  scale_sizes <- lapply(1:ncol(sizes), function(k){
    return(sizes[,k]/ det_timefn[k])
  }) %>% do.call(rbind,.)
  return(t(scale_sizes))
}



get_probs_nomutant <- function(params){
  alpha_vec = params$brate
  beta_vec = params$drate
  nu_vec <- params$mutrate
  lambda_vec = alpha_vec - beta_vec
  delta_vec = sapply(1:length(lambda_vec), function(k) max(lambda_vec[1:k]))
  rn_vec = sapply(1:length(lambda_vec), function(k) sum(lambda_vec[1:k] == delta_vec[k]))
  
  scale_tail_list <- get_scale_tail_vec(params)
  scale_vec <- scale_tail_list[[1]]
  thalf_vec = sapply(1:(length(scale_vec )-1), function(k){
    rel_thalf_kp1 =  1/delta_vec[k]*log (delta_vec[k]/(scale_vec[k]*nu_vec[k]*log( nu_vec[k] ^(-1/delta_vec[k])  )^(rn_vec[k]-1)   ))
    return(rel_thalf_kp1 )
  })
  
  
  prob_types_dont_exist <- sapply(1:length(thalf_vec), function(x){
    fintime <- params$fintime[x]
    return( ( 1+ exp(lambda_vec[1]*(fintime-thalf_vec[x])) )^(-1) )
  } )
  
  return(prob_types_dont_exist)
}


#### Likelihood functions ####

#note that mittag-leffler has 0 probability for  0 counts
get_loglikelihood_mutantcounts <- function(params,sizes){
  
  # k=131
  # params <- params
  # params$mutrate[1:(nrow(params)-1)] = rep(commonmu_grid[k],
  #                                                (nrow(params)-1))
  
  scale_tail_list <- get_scale_tail_vec(params)
  scale_vec <- scale_tail_list[[1]]
  tail_vec <- scale_tail_list[[2]]
  
  scaled_mutantcnts <- get_mutantcounts_scaledby_tfn(params,
                                                     sizes) 
  # scaled_mutantcnts_byscale_vect <- sapply(1:length(scale_vec), function(k){
  #   print((scaled_mutantcnts[,k]/scale_vec[k]) %>% mean)
  #   return(scaled_mutantcnts[,k]/scale_vec[k])
  # })
  # sum(dexp(scaled_mutantcnts_byscale_vect[,2],
  #          rate = 1,
  #          log =T))
  # 
  # sum(dexp(rexp(1000,rate = 1/0.0122584),
  #          rate = 1,
  #          log =T))
  # 
  loglike_per_type <- sapply(1:ncol(sizes), function(x){
    
    if (tail_vec[x] ==1){
      
      toreturn <- sum(dexp(scaled_mutantcnts[,x],
                           rate = 1/scale_vec[x],
                           log =T))
      return(toreturn)
    } else {
      
      toreturn <- sum(dml(scaled_mutantcnts[,x],
                          tail =tail_vec[x], 
                          scale = scale_vec[x],log = T))
      return(toreturn)
    }
    
  })
  
  return( loglike_per_type)
}

get_p0_loglikelihood <- function(params,sizes){
  
  probs_no_mutant <- get_probs_nomutant(params)
  #don't do p0 method of type 1 (always exists)
  loglikelihoods <- sapply(2:ncol(sizes), function(k){
    return(dbinom(sum(sizes[,k]==0),
                  size = nrow(sizes),
                  prob = probs_no_mutant[k-1],
                  log = T))
  })
  return(loglikelihoods)
}

#### Do likelihood inference ####
get_loglikelihood_vals_overmugrid <- function(mu_grid,params,sizes){
  loglikelihood_vals_over_mugrid <- lapply(1:length(mu_grid), function(k){
    
    relparams <- params
    relparams$mutrate[1:(nrow(relparams)-1)] = rep(mu_grid[k],
                                                   (nrow(relparams)-1))
    
    logLike_ML <- get_loglikelihood_mutantcounts(relparams,
                                                 sizes)
    logLike_p0 <- get_p0_loglikelihood(relparams,sizes)
    return(list(logLike_ML,logLike_p0))
  })
  return(loglikelihood_vals_over_mugrid)
}

get_mle_ci_fromlike <- function(loglike,mugrid){
  # loglike <-logLike_ML[1,]
  # mugrid <- commonmu_grid
  maxl <- max(loglike)
  index_mle <- which(loglike == maxl)
  normed_loglike <- loglike - maxl
  index_lowCI <- which((normed_loglike< -1.92)[1:(index_mle-1)]) %>% tail(.,n=1)
  index_upCI <- which((normed_loglike< -1.92)[(index_mle+1):length(loglike)])%>% 
    head(.,n=1)+index_mle
  mle =  mugrid[index_mle]
  lowCI = mugrid[index_lowCI]
  upCI = mugrid[index_upCI]
  
  if (length(mle) != 1){
    mle = NA
  } 
  if (length(lowCI) != 1){
    lowCI = NA
  } 
  if (length(upCI) != 1){
    upCI = NA
  } 
  
  
  toreturn <- c(mle,
                lowCI,
                upCI)
  
  return(toreturn)
}