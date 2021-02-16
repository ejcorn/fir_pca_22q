rm(list=ls())

args <- commandArgs(TRUE)
name_root <- args[1]
basedir <- args[2]
component_design <- args[3]
fin <- 6
st <- 1

# name_root <- 'CPCA_IDSchaefer200Z1xcp_6p_noFilter'
# basedir <- '~/Dropbox/Cornblath_Bassett_Projects/BrainStates22q/fir_pca_22q/'
# component_design <- 'ThreatNonthreatAllStimuliStratified'

setwd(basedir)

source(paste0(basedir,'code/miscfxns/packages.R'))
library(nlme)
masterdir <- paste(basedir,'results/',name_root,'/',sep='')

source(paste(basedir,'code/plottingfxns/GeomSplitViolin.R',sep=''))
source(paste0(basedir,'code/plottingfxns/plottingfxns.R'))
source(paste0(basedir,'code/miscfxns/miscfxns.R'))
source(paste0(basedir,'code/statfxns/statfxns.R'))

demo <- read.csv(paste(basedir,'data/Demographics',name_root,'.csv',sep=''),stringsAsFactors=F)
rownames(demo) <- as.character(demo$scanid)
grp.colors <- getGroupColors()

# load FIR betas

TR <- 3
covariates <- c('scanage_months','sex','BrainSegVol','idemo_meanrelrms','handedness')
savedir <- paste0(masterdir,'analyses/fir/cpc_timecourse_fin',fin,'st',st,'/',component_design,'/')
dir.create(paste0(savedir,'lme'),recursive=T)
stim.types <- list(threat=1,nonthreat=0)
response.types <- list(correct=1,incorrect=0)

# load mixed effects modeling results
load(file = paste0(savedir,'lme_all_trials/',component_design,'_PCTimeCoursesLMEEffectsCorrectThreat.RData'))
source('code/statfxns/lme_msfxns.R') # functions for mixed effects model selection

if(grepl('xcp_36p_despike',name_root)){
  ncomps <- 6
  w.multiplier <- 1
  results <- results[1:ncomps]
} else if(grepl('xcp_6p_noFilter',name_root)){
  ncomps <- 6
  w.multiplier <- 6/5
  results <- results[1:ncomps] 
}

t.axis <- seq(from=TR,to=TR*fin,by=TR) # labels for time on x-axis
all.p.ints.22q <- lapply(results, function(X) X$p.ints.22q)
all.p.ints.time <- lapply(results, function(X) X$mdl.best.stats[rownames(X$p.ints.correct.threat),])
p.labels <- lapply(all.p.ints.22q, function(X) paste0(paste0(rownames(X),': p = ',signif(X,2)),collapse='\n'))

# lists contain elements for each PC; looping through PCs and doing the same thing to each
# concatenate actual data with predicted values for each observation (i.e. 1 time point from 1 subject from 1 stimulus) 
# from fitted values that belong to model object
df.list.all <- lapply(results, function(X) cbind(X$mdl.best$data[complete.cases(X$mdl.best$data),],Score.fixed=X$mdl.best$fitted[,'fixed']))
for(response.type.name in names(response.types)){   # loop through response types -- chose to split on response type instead of stimulus type because more interactions involving response type
    response.type <- response.types[[response.type.name]]
    df.list <- lapply(df.list.all, function(X) X[X$Correct == response.type,]) # loop through PCs and subset data from the desired response type (1 level of my interactions)
    # average all of the variables, including both data and fitted values, over all levels of time, group, and stimulus type (after having subsetted by response type)
    df.grp.time.mean <- lapply(df.list, function(df) aggregate(df,by=list(df$Time,df$Is22q,df$Threat),mean)) 
    # this contains the group average data and group average fitted values (should be same as computing fitted values with random effects set to 0, but easier to compute this way)
    for(j in 1:length(results)){ # overwrite these new variables in the original results list
      results[[j]]$grp.average.pred <- df.grp.time.mean[[j]]
      results[[j]]$df.data <- df.list[[j]]
    }
    # this function attempts to compute the group average fitted line by ignoring random effects... may be deprecated
    #for(j in 1:length(results)){results[[j]]$grp.average.pred <- lme.fixed.predict.ave(results[[j]]$mdl.best,c('Is22q','Time'))}
    
    # plot data -- this never goes in paper. plots actual data, not fitted values
    p.data.all <- lapply(results, function(X) ggplot() + geom_line(data= X$df.data[!outlier.mask(X$df.data$Score),],aes(x = Time,y = Score, group = paste0(scanid,Threat),color=as.character(Is22q)),alpha=0.2,size=0.25) + 
                           geom_line(data=X$grp.average.pred,aes(x=Time,y=Score,color=as.character(Is22q),linetype=as.character(Threat)),size=1)+
                           #ggtitle(X$PC)+ xlab('Time (s)') + ylab('Estimated HDR')+
                           ggtitle(NULL) + xlab(NULL) + ylab(NULL) + # easier to just add titles in illustrator
                           scale_x_continuous(limits=c(TR,fin*TR),breaks=t.axis,labels = t.axis) +
                           scale_color_manual(values=unname(grp.colors),limits=c('0','1'),guide='none')+
                           scale_linetype_manual(guide='none',values=c('dashed','solid'),limits=c('0','1'))+ # make threat a solid line
                           theme_classic() + theme(plot.title = element_text(hjust=0.5),text=element_text(size=8),plot.margin = unit(c(0,0,0,0),'cm')))
    
    # plot predicted values and add p-value annotation
    # first geom_line command plots the subject specific fitted values (e.g. results[[2]]$df.data, same as df.list[[2]] above)
    # second geom_line command plots the group average fitted values (e.g. results[[2]]$grp.average.pred, same as df.grp.time.mean[[2]] above)
    p.pred.all<- lapply(results, function(X) ggplot() + geom_line(data= X$df.data,aes(x = Time,y = pred.best, group = paste0(scanid,Threat),color=as.character(Is22q)),alpha=0.2,size=0.25) + 
                          geom_line(data=X$grp.average.pred,aes(x=Time,y=Score.fixed,color=as.character(Is22q),linetype=as.character(Threat)),size=1)+
                          #ggtitle(paste0(X$PC))+xlab('Time (s)')+ylab('Fitted Value')+
                          ggtitle(NULL) + xlab(NULL) + ylab(NULL) + # easier to just add titles in illustrator
                          scale_x_continuous(limits=c(TR,fin*TR),breaks=t.axis,labels = t.axis) +
                          scale_color_manual(values=grp.colors,limits=c('0','1'),guide='none')+
                          scale_linetype_manual(guide='none',values=c('dashed','solid'),limits=c('0','1'))+ # make threat a solid line
                          theme_classic() + theme(plot.title = element_text(hjust=0.5),text=element_text(size=8),plot.margin = unit(c(0,0,0,0),'cm')))
    #p.pred.all <- lapply(1:length(p.pred.all), function(PC) p.pred.all[[PC]] + annotate(geom='text',x=Inf,y=-Inf,hjust=1,vjust=-0.5,label=p.labels[[PC]],size=1))
    
    p.all <- plot_grid(plotlist = c(p.data.all,p.pred.all),align='hv',ncol=ncomps)
    #p.data <- plot_grid(plotlist = p.data.all,align='hv',nrow=1)
    #p.all <- plot_grid(plotlist = list(p.pred,p.data),nrow = 2)
    ggsave(plot = p.all, filename = paste0(savedir,'lme_all_trials/',response.type.name,'_PCTimeCoursesObsPred.pdf'),
           width = 17*w.multiplier, height = 6.5, units = "cm",useDingbats=F) # if you want to use all components can do this to keep scale
    p.pred.all.col <- plot_grid(plotlist = p.pred.all,ncol=1,align='hv')
    p.data.all.col <- plot_grid(plotlist = p.data.all,ncol=1,align='hv')
    p.all <- plot_grid(plotlist = list(p.data.all.col,p.pred.all.col),ncol=2)
    ggsave(plot = p.all, filename = paste0(savedir,'lme_all_trials/',response.type.name,'_PCTimeCoursesObsPredVertical.pdf'),
           width =6.5, height = 17*w.multiplier, units = "cm",useDingbats=F) # if you want to use all components can do this to keep scale
    
  }

