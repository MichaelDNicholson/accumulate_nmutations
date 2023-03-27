##take in output of python simulations make hitting time plots

library(feather)
library(magrittr)
library(ggplot2)
library(cowplot)
library(reshape2)

projectRoot <-  "/home/mnichol3/ownCloud/accumulate_nmutations/"
plotdir = "/home/mnichol3/ownCloud/accumulate_nmutations/images/"
dirSimData <- paste0(projectRoot,"results/simout_arrival_times/")

#source plotting functions
source(paste0(projectRoot,"code/functionsPlotHittingTimes.R"))

listDensPlots <- lapply(list.files(dirSimData), function(f){
  params_times = read_feather(paste0(dirSimData,f))
  return(get_density_plot( params_times))
})

listCdfPlots <- lapply(list.files(dirSimData), function(f){
  params_times = read_feather(paste0(dirSimData,f))
  return(get_cdf_plot( params_times))
})

listCdfPlotsNoLeg <- lapply(listCdfPlots, function(p){
  return(p+theme(legend.position = "none"))
})
cdfPlotLegendTop <- get_legend(listCdfPlots[[1]]+
                              theme(legend.position = "top",
                                    legend.text=element_text(size=12)))


#save in this order so plots with same birth/death parameters but different
#mutation rates are in B D
gridPlotCdfsAllNoLeg <- plot_grid( plotlist= listCdfPlotsNoLeg[c(1,3,2,4)],
                                 nrow = 2,
                                 labels = c('A',
                                            'B',
                                            'C',
                                            'D'),
                                 vjust = -1)

gridPlotCdfsAllLeg <- plot_grid(cdfPlotLegendTop,gridPlotCdfsAllNoLeg,
                                  rel_heights = c(.2,1),
                                  nrow = 2)


if (saveplots == T){
  
  if (dir.exists(plotdir)==F){
    dir.create(plotdir)
  }
  
  #save density plots
  #just save third density plot
  indexSaveDens <- 3
  densNameId = (list.files(dirSimData)[indexSaveDens] %>% strsplit(.,"[.]"))[[1]][1]
  plot_dens_filename = paste(plotdir,
                                "dens234_", densNameId,
                                ".pdf",sep ="")
  
  save_plot(filename =  plot_dens_filename,
            listDensPlots[[indexSaveDens]],
            base_height = 4,
            base_asp = 2.5)


#save cdf plots

  plot_cdf_filename = paste(plotdir,
                             "cdf234_all",
                             ".pdf",sep ="")
  
  save_plot(filename =    plot_cdf_filename,
            gridPlotCdfsAllLeg,
            base_height = 4.2,
            base_asp = 1.6)
}