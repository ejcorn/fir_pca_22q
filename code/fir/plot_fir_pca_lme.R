rm(list=ls())

args <- commandArgs(TRUE)
name_root <- args[1]
basedir <- args[2]

name_root <- 'ScanIDSchaefer200Z1_22q_PreQC_XCP36pdespike_us_100reps'
basedir <- '~/Dropbox/Cornblath_Bassett_Projects/BrainStates22q/brain_states_22q/'

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
component_design <- 'ThreatNonthreatAllStimuliStratified'
ncomps <- 6
fin <- 6
TR <- 3
covariates <- c('scanage_months','sex','BrainSegVol','idemo_meanrelrms','handedness')
t.names <- as.character(TR*1:fin)
savedir <- paste0(masterdir,'analyses/fir/subject_fir_correct_incorrect_pca/cpc_timecourse/',component_design,'/')
dir.create(paste0(savedir,'lme'),recursive=T)
stim.types <- list(all='all',threat='threat',nonthreat='nonthreat')
response.types <- list(correct='correct',incorrect='incorrect')

results <- list()
load(file = paste0(savedir,'lme/',component_design,'_PCTimeCoursesLMEEffects.RData'))

for(stim.type in stim.types){
  for(response.type in response.types){
    # correct group comparisons for multiple comparisons
    results.sr <- results[[stim.type]][[response.type]]
    all.p.ints <- list.posthoc.correct(lapply(results.sr, function(X) X$p.ints),'fdr')
    all.p.ints <- lapply(results.sr, function(X) X$p.int)
    p.labels <- lapply(all.p.ints, function(X) paste0(paste0(rownames(X),': p = ',signif(X,2)),collapse='\n'))
    
    df.list <- lapply(results.sr, function(X) cbind(X$mdl.best$data[complete.cases(X$mdl.best$data),],Score.fixed=X$mdl.best$fitted[,'fixed']))
    df.grp.time.mean <- lapply(df.list, function(df) aggregate(df,by=list(df$Time,df$Is22q),mean))
    for(j in 1:length(results.sr)){results.sr[[j]]$grp.average.pred <- df.grp.time.mean[[j]]}
    #for(j in 1:length(results.sr)){results.sr[[j]]$grp.average.pred <- lme.fixed.predict.ave(results.sr[[j]]$mdl.best,c('Is22q','Time'))}
    
    # plot data
    p.data.all <- lapply(results.sr, function(X) ggplot() + geom_line(data= X$mdl.best$data,aes(x = Time,y = Score, group = scanid,color=as.character(Is22q)),alpha=0.2,size=0.25) + 
      geom_line(data=X$grp.average.pred,aes(x=Time,y=Score,color=as.character(Is22q)),size=1)+
      #ggtitle(X$PC)+ xlab('Time (s)') + ylab('Estimated HDR')+
      ggtitle(NULL) + xlab(NULL) + ylab(NULL) + # easier to just add titles in illustrator
      scale_color_manual(values=grp.colors,limits=c('0','1'),guide='none')+
      theme_classic() + theme(plot.title = element_text(hjust=0.5),text=element_text(size=8),plot.margin = unit(c(0,0,0,0),'cm')))
    
    # plot predicted values and add p-value annotation
    p.pred.all<- lapply(results.sr, function(X) ggplot() + geom_line(data= X$mdl.best$data,aes(x = Time,y = pred.best, group = scanid,color=as.character(Is22q)),alpha=0.2,size=0.25) + 
                          geom_line(data=X$grp.average.pred,aes(x=Time,y=Score.fixed,color=as.character(Is22q)),size=1)+
                          #ggtitle(paste0(X$PC))+xlab('Time (s)')+ylab('Fitted Value')+
                          ggtitle(NULL) + xlab(NULL) + ylab(NULL) + # easier to just add titles in illustrator
                          scale_color_manual(values=grp.colors,limits=c('0','1'),guide='none')+
                          theme_classic() + theme(plot.title = element_text(hjust=0.5),text=element_text(size=8),plot.margin = unit(c(0,0,0,0),'cm')))
    p.pred.all <- lapply(1:length(p.pred.all), function(PC) p.pred.all[[PC]] + annotate(geom='text',x=Inf,y=-Inf,hjust=1,vjust=-0.5,label=p.labels[[PC]],size=1))
    
    p.all <- plot_grid(plotlist = c(p.data.all,p.pred.all),align='hv',nrow=2)
    #p.data <- plot_grid(plotlist = p.data.all,align='hv',nrow=1)
    #p.all <- plot_grid(plotlist = list(p.pred,p.data),nrow = 2)
    ggsave(plot = p.all, filename = paste0(savedir,'lme/',stim.type,response.type,'_PCTimeCoursesObsPred.pdf'),
           width = 17.5,height = 6.5, units = "cm",useDingbats=F)
    
  }
}

