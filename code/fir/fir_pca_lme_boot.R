rm(list=ls())

args <- commandArgs(TRUE)
name_root <- args[1]
basedir <- args[2]
fin <- 6
st <- 1

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

TR <- 3
covariates <- c('scanage_months','sex','BrainSegVol','idemo_meanrelrms','handedness')
t.names <- as.character(TR*1:fin)
savedir <- paste0(masterdir,'analyses/fir/cpc_timecourse_fin',fin,'st',st,'/',component_design,'/')
dir.create(paste0(savedir,'lme_all_trials/'),recursive=T)
stim.types <- list(threat='threat',nonthreat='nonthreat')
response.types <- list(correct='correct',incorrect='incorrect')

matData <- lapply(stim.types, function(S)
  lapply(response.types, function(R) readMat(paste0(savedir,S,R,'FIRBetas_CPCScores.mat'))))

source('code/statfxns/lme_msfxns.R') # functions for mixed effects model selection
results <- data.list <- list()

# load data from threat nonthreat, correct incorrect, then test for an interaction with those factors
for(stim.type in stim.types){
  data.list[[stim.type]] <- list()
  for(response.type in response.types){
    PC.scores <- matData[[stim.type]][[response.type]]$betaFIR.reshape
    scanids.mat <- as.numeric(unname(unlist(matData[[stim.type]][[response.type]]$subjInd.scanID)))
    data.list[[stim.type]][[response.type]] <- list()
    for(PC in 1:ncomps){
      PC.scores.j <- setNames(as.data.frame(PC.scores[,,PC]),t.names)
      PC.scores.j <- cbind(PC.scores.j,scanid=scanids.mat,Is22q=as.numeric(demo[as.character(scanids.mat),'study']=='22q'),stringsAsFactors=F)
      PC.scores.j <- cbind(PC.scores.j,demo[as.character(scanids.mat),covariates])
      PC.df <- collapse.columns(PC.scores.j,cnames = t.names,groupby = c('scanid','Is22q',covariates))
      colnames(PC.df)[1:2] <- c('Score','Time')
      PC.df$Time <- as.numeric(PC.df$Time)
      
      for(td in 1:5){PC.df[,paste0('Time',td)] <- PC.df$Time^td}
      data.list[[stim.type]][[response.type]][[PC]] <- PC.df
    }
  }
}

# load results of full sample model selection
load(file = paste0(savedir,'lme_all_trials/',component_design,'_PCTimeCoursesLMEEffectsCorrectThreat.RData'))
source('code/statfxns/lme_msfxns.R') # functions for mixed effects model selection

iterator <- as.data.frame(t(expand.grid(stim=stim.types,resp=response.types)))
PC <- 1
# try to refit the model identified by the model selection procedure in null data
nboots <- 1000
results.boot <- list()
for(PC in 1:ncomps){
  # extract data on one PC at a time from all stimuli
  df.all <- do.call(rbind,lapply(iterator, function(X) data.frame(data.list[[X$stim]][[X$resp]][[PC]],
                                                                  Threat=as.numeric(X$stim=='threat'),Correct=as.numeric(X$resp == 'correct')))) # add binary indicators for correct and threat
  # resample subjects with replacement
  ids <- unique(df.all$scanid)
  ids <- setNames(ids,ids)
  df.all.by.id <- lapply(ids, function(id) df.all[df.all$scanid==id,]) # split df by id
  # now randomly choose ids with replacement and reconstruct a bunch of new data frames
  df.all.sample <- lapply(1:nboots, function(x) do.call(rbind,lapply(sample(ids,replace=T),function(X) df.all.by.id[[as.character(X)]])))
  #sapply(df.all.sample,function(X) length(unique(X$scanid))) # confirm there are slightly different numbers of unique scan ids in each resampled data frame
  # refit model on bootstrap resampled data
  results.boot[[PC]] <- list()
  for(P in 1:nboots){
    print(paste('boot',P))
    mdl.boot <- tryCatch({mdl.boot <- lme.ms(results[[PC]]$cfg.best,df.all.sample[[P]])},
                         error=function(err) {print('unable to fit');return('asdf')})
    if(is.object(mdl.boot)){results.boot[[PC]][[P]] <- mdl.boot
    } else {results.boot[[PC]][[P]] <- 'unable to fit'}
  }
}


boot.pvals <- boot.plots <- boot.stats.list<- list()
for(PC in 1:ncomps){
  boot.stats <- results.boot[[PC]][sapply(results.boot[[PC]],is.list)] # delete bootstraps unable to fit
  boot.stats.list[[PC]] <- lapply(boot.stats, function(X) summary(X)$tTable[,'Value',drop=F])
  boot.dist <- do.call(cbind,boot.stats.list[[PC]])
  boot.pvals[[PC]] <- sapply(rownames(boot.dist),function(j) pval.2tail.np(0,boot.dist[j,]))
}
# tried to put the below snippet in a loop and it wasn't working
boot.plots <- lapply(1:ncomps, function(PC) ggplot() + geom_text(aes(x=0,y=1:nrow(results[[PC]]$mdl.best.stats),label=rownames(results[[PC]]$mdl.best.stats),
                                                                     color=factor(ifelse(boot.pvals[[PC]]<(0.05/ncomps),yes='sig',no='nonsig'),levels = c('sig','nonsig'))),size=2)+
                       scale_color_manual(values=c('red','grey50'),guide='none')+ggtitle(paste0('PC',PC))+
                       theme_void() + xlab('') + ylab(expression(beta)) + theme(text=element_text(size=8),
                                                                                   plot.title = element_text(hjust=0.5,size=8,face='bold')))
p.all <- plot_grid(plotlist=boot.plots,nrow=2,ncol =3, align='hv')
ggsave(plot = p.all, filename = paste0(savedir,'lme_all_trials/',component_design,'_LMEBetaBootstrap.pdf'),
       width = 18,height = 12, units = "cm",useDingbats=F)
save(boot.stats.list,file = paste0(savedir,'lme_all_trials/',component_design,'_LMEBetaBootstrap.RData'))
