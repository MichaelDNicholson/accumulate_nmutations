#functions to take in hitting times of multitype branching
#process simulations and compare against analytic results
#takes in feather input from python sims


projectRoot <-  "/home/mnichol3/ownCloud/accumulate_nmutations/"
source(paste0(projectRoot,"code/functionsCoreAccumulateNmuts.R"))

get_legend<-function(myggplot){
  tmp <- ggplot_gtable(ggplot_build(myggplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)
}

#analytic functions


get_analytic_thalf <- function(params){
  
  alpha_vec = params$brate
  beta_vec = params$drate
  nu_vec = params$mutrate
  lambda_vec = alpha_vec - beta_vec
  delta_vec = sapply(1:length(lambda_vec), function(k) max(lambda_vec[1:k]))
  rn_vec = sapply(1:length(lambda_vec), function(k) sum(lambda_vec[1:k] == delta_vec[k]))
  
  scale_tail_list <- get_scale_tail_vec(params)
  scale_vec <- scale_tail_list[[1]]
  thalf_vec = sapply(1:(length(scale_vec )-1), function(k){
    rel_thalf_kp1 =  1/delta_vec[k]*log (delta_vec[k]/(scale_vec[k]*nu_vec[k]*log( nu_vec[k] ^(-1/delta_vec[k])  )^(rn_vec[k]-1)   ))
    return(rel_thalf_kp1 )
  })
  return(thalf_vec)
}

#### density functions ####

time_density = function(thalf,delta1){
  function(x){
    
    return(delta1*exp(delta1*(x-thalf))/(1+exp(delta1*(x-thalf)))^2)
    
  }
}

get_density_plot <- function(params_times){
  
  params = params_times[,1:3]
  colnames(params) = c('brate','drate','mutrate')
  hitting_times = params_times[,4:ncol(params_times)]
  hitting_time2 = hitting_times[2,] %>% as.numeric()
  hitting_time3 = hitting_times[3,] %>% as.numeric()
  hitting_time4 = hitting_times[4,] %>% as.numeric()
  
  delta1 <- params$brate[1] - params$drate[1]
  
  thalf_vec <- get_analytic_thalf(params)
  
  dens2_fn = time_density(thalf_vec[1], delta1)
  dens3_fn = time_density(thalf_vec[2], delta1)
  dens4_fn = time_density(thalf_vec[3], delta1)
  
  df_htimes234 = data.frame("times2" =hitting_time2 ,
                            "times3" =hitting_time3 ,
                            "times4" =hitting_time4 )
  
  
  colors = c("red","gold1","#009E73")
  plot_dens234 = ggplot(df_htimes234)+
    geom_histogram(color="black", fill=colors[1],
                   alpha = .5,
                   binwidth = 1,
                   aes(x=times2,y=..density..))+
    stat_function(fun = dens2_fn,
                  lwd = 2, 
                  col = colors[1],
                  aes(x= times2))+
    geom_histogram(color="black", fill=colors[2],
                   alpha = .5,
                   binwidth = 1,
                   aes(x=times3,y=..density..))+
    stat_function(fun = dens3_fn,
                  lwd = 2, 
                  col = colors[2],
                  aes(x= times3))+
    theme_classic()+
    theme(axis.text = element_text(size = 14),
          axis.title = element_text(size=18))+
    xlab('Time (t)')+
    ylab(paste('Density of arrival time','\n',
               'for type n cell',sep=""))+
    geom_histogram(color="black", fill=colors[3],
                   binwidth = 1, alpha = .5,
                   aes(x=times4,y=after_stat(density)))+
    stat_function(fun = dens4_fn,
                  lwd = 2, 
                  col = colors[3],
                  aes(x= times4))
  
  return(plot_dens234)
}


#### cdf functions ####


time_cdf = function(thalf,delta1){
  function(x){
    
    return(1-1/(1+exp(delta1*(x-thalf))))
    
  }
}


get_cdf_plot <- function(params_times){
  
  params = params_times[,1:3]
  colnames(params) = c('brate','drate','mutrate')
  hitting_times = params_times[,4:ncol(params_times)]
  hitting_time2 = hitting_times[2,] %>% as.numeric()
  hitting_time3 = hitting_times[3,] %>% as.numeric()
  hitting_time4 = hitting_times[4,] %>% as.numeric()
  nsims <- ncol(hitting_times)
  
  #get cdfs from sims
  xmin = min(hitting_time2)-1
  xmax = max(hitting_time4)
  xgrid = seq(xmin,xmax,by=.1)
  
  cdf_time2 = sapply(xgrid, function(val) sum(hitting_time2<=val )/nsims)
  cdf_time3 = sapply(xgrid, function(val) sum(hitting_time3<=val )/nsims)
  cdf_time4 = sapply(xgrid, function(val) sum(hitting_time4<=val )/nsims)
  
  thalf_vec <- get_analytic_thalf(params)
  delta1 <- params$brate[1] - params$drate[1]
  
  analytic2_cdf = time_cdf(thalf_vec[1],delta1)(xgrid)
  analytic3_cdf = time_cdf(thalf_vec[2],delta1)(xgrid)
  analytic4_cdf = time_cdf(thalf_vec[3],delta1)(xgrid)
  
  
  df_cdfs = data.frame("xgrid" = xgrid,
                       "cdf_time2" = cdf_time2,
                       "cdf_time3" = cdf_time3,
                       "cdf_time4" = cdf_time4,
                       "analytic2_cdf" = analytic2_cdf,
                       "analytic3_cdf" = analytic3_cdf,
                       "analytic4_cdf"= analytic4_cdf)
  
  df_cdfs_long = melt(df_cdfs, id = "xgrid")
  df_cdfs_long$type = "Theory"
  logic_sims = sapply(df_cdfs_long$variable, function(str){
    str %>% as.character %>% grepl( "time",., fixed = TRUE)
  })
  df_cdfs_long$type[logic_sims] = "Simulation"
  
  logic_t3 = sapply(df_cdfs_long$variable, function(str){
    str %>% as.character %>% grepl( "3",., fixed = TRUE)
  })
  logic_t4 = sapply(df_cdfs_long$variable, function(str){
    str %>% as.character %>% grepl( "4",., fixed = TRUE)
  })
  df_cdfs_long$group = "2"
  df_cdfs_long$group[logic_t3] = "3"
  df_cdfs_long$group[logic_t4] = "4"
  
  
  
  plot_cdfs = ggplot( df_cdfs_long,
                      aes(x=xgrid,color =group,y=value,linetype=type))+
    geom_line(linewidth= .75) +
    theme_classic()+
    xlab('Time (t)')+
    ylab(paste("Probability type n cell",'\n',
               'arrived by  t',sep = ""))+
    scale_color_manual("kk",  
                       values = colors,
                       labels=c("n = 2","n = 3","n = 4"))   + 
    scale_linetype_manual("",values=c(1,2),
                          labels=c("Theory","Simulation")) +
    theme(legend.justification = c(0, 1), 
          legend.position = c(.75, .9),
          legend.title=element_blank(),
          legend.text=element_text(size=10),
          legend.key.width=unit(.8,"cm"))
  
  return(plot_cdfs)
  
}