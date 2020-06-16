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
#demo$study <- sample(demo$study,replace = F) # shuffle group labels and refit as a sanity check permutation null
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

matData <- lapply(stim.types, function(S)
  lapply(response.types, function(R) readMat(paste0(savedir,S,R,'FIRBetas_CPCScores.mat'))))

source('code/statfxns/lme_msfxns.R') # functions for mixed effects model selection
results <- list()

for(stim.type in stim.types){
  results[[stim.type]] <- list()
 for(response.type in response.types){
    results[[stim.type]][[response.type]] <- list()
    PC.scores <- matData[[stim.type]][[response.type]]$betaFIR.reshape
    scanids.mat <- as.numeric(unname(unlist(matData[[stim.type]][[response.type]]$subjInd.scanID)))
    p.data.all <- p.pred.all <- list()
    for(PC in 1:ncomps){
      results.PC <- list()
      results.PC$PC <- paste0('PC',PC)
      PC.scores.j <- setNames(as.data.frame(PC.scores[,,PC]),t.names)
      PC.scores.j <- cbind(PC.scores.j,scanid=scanids.mat,Is22q=as.numeric(demo[as.character(scanids.mat),'study']=='22q'),stringsAsFactors=F)
      PC.scores.j <- cbind(PC.scores.j,demo[as.character(scanids.mat),covariates])
      PC.df <- collapse.columns(PC.scores.j,cnames = t.names,groupby = c('scanid','Is22q',covariates))
      colnames(PC.df)[1:2] <- c('Score','Time')
      PC.df$Time <- as.numeric(PC.df$Time)
      
      # plot raw data of PC time courses
      p.data.all[[PC]] <- ggplot() + geom_line(data= PC.df,aes(x = Time,y = Score, group = scanid,color=as.character(Is22q)),alpha=0.5,size=0.25) + 
        ggtitle(paste0('PC',PC,':'))+ xlab('Time (s)') + ylab('Estimated HDR')+
        scale_color_manual(values=grp.colors,limits=c('0','1'),guide='none')+
        theme_bw() + theme(plot.title = element_text(hjust=0.5),text=element_text(size=8))
      
      
      for(td in 1:5){PC.df[,paste0('Time',td)] <- PC.df$Time^td}
      # Procedure for model selection: add time, add random effects, then add covariates. each time, deciding whether to add using LR test p <0.05
      # uses lme.ms, lme.compare, lme.stepup in lme_msfxns.R
      
      #1. Base model, start with no time, just covariates and group
      cfg.base <- list(fixed=c(covariates,'Is22q'),random='1',response='Score',id='scanid')
      #cfg.base <- list(fixed='1',random='1',response='Score',id='scanid')
      mdl.gold <- lme.ms(cfg.base,PC.df) # in every step you are comparing a gold standard model to a test model
      # 2. Add increasing order polynomials of time, fixed effects
      t.poly.max <- 5
      test.time.cfg <- lapply(1:t.poly.max, function(o.max) list(fixed=c(cfg.base$fixed,paste0('Time',1:o.max)), random=c(cfg.base$random),response=cfg.base$response,id=cfg.base$id)) # define sequential list of models to iterate through
      # start at first and advance as long as : (c1) p < 0.05 for anova(m1,m2) AND (c2) AIC drops AND (c3) added coefficient is significant
      list[mdl.gold,cfg.gold] <- lme.selectbest(mdl.gold,cfg.base,test.time.cfg)
      
      t.poly.max.sig <- 0
      if(any(grepl('Time',cfg.gold$fixed))){
        # 3. now add random effects of increasing time order up to max fixed effect
        t.poly.max.sig <- max(substr(cfg.gold$fixed[grepl('Time',cfg.gold$fixed)],5,5)) # only test interactions up to the maximum time you've selected
        test.time.ran.cfg <- lapply(1:t.poly.max.sig, function(o.max) list(fixed=cfg.gold$fixed, random=c(cfg.gold$random,paste0('Time',1:o.max)),response=cfg.gold$response,id=cfg.gold$id)) # define sequential list of models to iterate through
        list[mdl.gold,cfg.gold] <- lme.stepup(mdl.gold,test.time.ran.cfg)

        # 4. now test additions of fixed effect interactions with Is22q and Time, only if sig random effects
        test.time.by.22q.cfg <- lapply(1:t.poly.max.sig, function(o.max) list(fixed=c(cfg.gold$fixed,paste0('Time',1:o.max,'*Is22q')), random=c(cfg.gold$random),response=cfg.gold$response,id=cfg.gold$id)) # define sequential list of models to iterate through
        list[mdl.gold,cfg.gold] <- lme.selectbest(mdl.gold,cfg.gold,test.time.by.22q.cfg)
        #lme.ms(test.time.by.22q.cfg[[1]],PC.df)
      }
      
      results.PC$mdl.best.max.tsig <- t.poly.max.sig
      print(paste0(stim.type,'-',response.type,', PC',PC,', best model maximum time order: ',results.PC$mdl.best.max.tsig))
      results.PC$mdl.best <- mdl.gold
      results.PC$cfg.best <- cfg.gold
      
      # get predicted values from best model
      results.PC$mdl.best$data$pred.best <- predict(results.PC$mdl.best)
      
      # get stats for interactions, i.e. group differences, from best model
      results.PC$mdl.best.stats <- summary(results.PC$mdl.best)$tTable
      find.ints <- c(grep('Time*:Is22q',rownames(results.PC$mdl.best.stats)),grep('Is22q:Time*',rownames(results.PC$mdl.best.stats)))
      results.PC$p.ints <- results.PC$mdl.best.stats[find.ints,'p-value',drop=FALSE] # p-values for interactions
      
      results[[stim.type]][[response.type]][[PC]] <- results.PC
    }
    }
}

save(results, file = paste0(savedir,'lme/',component_design,'_PCTimeCoursesLMEEffects.RData'))
