#script to assess whether mutation rate can be inferred on 
#simulated data
library(feather)
library(magrittr)
library(ggplot2)
library(reshape2)
library(cowplot)
library(MittagLeffleR)



projectRoot <-  "/home/mnichol3/ownCloud/accumulate_nmutations/"
plotdir = paste0(projectRoot,"images/")
fluctSimsDir <- paste0(projectRoot,"results/simout_fluctuation_assay/")
saveplots = F #set to F, for running interactively

source(paste0(projectRoot,"code/functionsCoreAccumulateNmuts.R"))

setwd(fluctSimsDir)
input_files <- list.files()
mugrid <- rev(10^(-seq(0,4.5,.01))) 



stats_over_files <- lapply(input_files, function(x){
  params_sizes = arrow::read_feather(x)
  nsims = ncol(params_sizes) - 4
  params = params_sizes[,1:4]
  sizes = params_sizes[,5:ncol(params_sizes)] %>% t()
  colnames(params) = c('brate','drate','mutrate','fintime')
  ntype <- ncol(sizes)
  
  list_loglikes <- get_loglikelihood_vals_overmugrid(mugrid,params,sizes)
  logLikes_ML_fintype <- sapply(list_loglikes, function(l) l[[1]][ntype])
  logLikes_p0_fintype <- sapply(list_loglikes, function(l) l[[2]][ntype-1])
  
  stats_p0 <- get_mle_ci_fromlike(logLikes_p0_fintype,mugrid)
  stats_ML <- get_mle_ci_fromlike(logLikes_ML_fintype,mugrid)
  
  stats_p0_mu <- c(stats_p0,params$mutrate[1])
  stats_ML_mu <- c(stats_ML,params$mutrate[1])
  

  return(list(stats_p0_mu,stats_ML_mu))
  
}) 

stats_over_files_p0 = lapply(stats_over_files, function(x)x[[1]]) %>% do.call(rbind,.) %>% as.data.frame()
stats_over_files_ML = lapply(stats_over_files, function(x)x[[2]]) %>% do.call(rbind,.) %>% as.data.frame()

colnames(stats_over_files_p0) = c("mle","ci_low", "ci_up","nu")
colnames(stats_over_files_ML) = c("mle","ci_low", "ci_up","nu")



stats_over_files_p0_log <- stats_over_files_p0 %>% log10
plt_stats_p0method <- ggplot(stats_over_files_p0_log)+
  geom_point(aes(x=nu,y=mle))+
  geom_line(aes(x=nu,y=nu),linetype = "dashed")+
  geom_errorbar(aes(x=nu,y=mle,ymin=ci_low, ymax=ci_up), width=.1,
                position=position_dodge(0.05))+
  theme_bw()+
  xlab(bquote(True~log[10](nu)))+
  ylab(bquote(Inferred~log[10](nu)))+
  geom_segment(aes(x = -3, y = -1, xend = -2.8, yend = -1),
               linetype = "dashed")+
  geom_text(x=-2.68, y=-1, label="y=x")
  # ggtitle("P0 method")


stats_over_files_ML_log <- stats_over_files_ML %>% log10
plt_stats_MLmethod <- ggplot(stats_over_files_ML_log)+
  geom_point(aes(x=nu,y=mle))+
  geom_line(aes(x=nu,y=nu),linetype = "dashed")+
  geom_errorbar(aes(x=nu,y=mle,ymin=ci_low, ymax=ci_up), width=.1,
                position=position_dodge(0.05))+
  theme_bw()+
  xlab(bquote(True~log[10](nu)))+
  ylab(bquote(Inferred~log[10](nu)))+
  geom_segment(aes(x = -3, y = -1.75, xend = -2.8, yend = -1.75),
               linetype = "dashed")+
  geom_text(x=-2.68, y=-1.75, label="y=x")
  # ggtitle("Maximum likelihood on\nmutant counts")

plot_fluct_inference_both <- plot_grid(plt_stats_p0method,plt_stats_MLmethod,
                                       labels = c("A","B"))



if (saveplots == T){
  if (dir.exists(plotdir)==F){
    dir.create(plotdir)
  }
  
  plt_stats_p0method_filename <- paste0(plotdir,"plt_stats_p0method.pdf")
  save_plot(plt_stats_p0method_filename, plt_stats_p0method,
            base_asp = 1.6,base_height = 2.3)
  
  plt_stats_MLmethod_filename <- paste0(plotdir,"plt_stats_MLmethod.pdf")
  save_plot(plt_stats_MLmethod_filename, plt_stats_MLmethod,
            base_asp = 1.6,base_height = 2.3)
  
  plot_fluct_inference_both_filename <- paste0(plotdir,"plot_fluct_inference_both.pdf")
  save_plot(plot_fluct_inference_both_filename, plot_fluct_inference_both ,
           base_height = 2.5,base_width = 3*2 )
  
}
 
