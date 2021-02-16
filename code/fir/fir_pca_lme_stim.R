rm(list=ls())

args <- commandArgs(TRUE)
name_root <- args[1]
basedir <- args[2]
component_design <- args[3]

# basedir <- '~/Dropbox/Cornblath_Bassett_Projects/BrainStates22q/fir_pca_22q/'
name_root <- 'CPCA_IDSchaefer200Z1xcp_6p_noFilter'
basedir <- '/cbica/home/cornblae/ecornblath/fir_pca_22q/'
component_design <- 'ThreatNonthreatAllStimuliStratified'

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
ncomps <- 7
fin <- 6
TR <- 3
covariates <- c('scanage_months','sex','BrainSegVol','idemo_meanrelrms','handedness')
t.names <- as.character(TR*1:fin)
savedir <- paste0(masterdir,'analyses/fir/subject_fir_correct_incorrect_pca/cpc_timecourse/',component_design,'/')
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

iterator <- as.data.frame(t(expand.grid(stim=stim.types,resp=response.types)))
PC <- 4
for(PC in 1:ncomps){
  results.PC <- list()
  results.PC$PC <- paste0('PC',PC)
  # extract data on one PC at a time from all stimuli
  df.all <- do.call(rbind,lapply(iterator, function(X) data.frame(data.list[[X$stim]][[X$resp]][[PC]],
                                                                  Threat=as.numeric(X$stim=='threat'),Correct=as.numeric(X$resp == 'correct')))) # add binary indicators for correct and threat
  #1. Base model, start with no time
  cfg.base <- list(fixed=c(covariates,'Is22q'),random='1',response='Score',id='scanid')
  #cfg.base <- list(fixed='1',random='1',response='Score',id='scanid')
  mdl.gold <- lme.ms(cfg.base,df.all) # in every step you are comparing a gold standard model to a test model
  # 2. Add increasing order polynomials of time
  t.poly.max <- 5
  test.time.cfg <- lapply(1:t.poly.max, function(o.max) list(fixed=c(cfg.base$fixed,paste0('Time',1:o.max)), random=c(cfg.base$random),response=cfg.base$response,id=cfg.base$id)) # define sequential list of models to iterate through
  # start at first and advance as long as : (c1) p < 0.05 for anova(m1,m2) AND (c2) AIC drops AND (c3) added coefficient is significant
  list[mdl.gold,cfg.gold] <- lme.selectbest(mdl.gold,cfg.base,test.time.cfg)
  
  # 3. do a step-down model selection to add 2-way and 3-way correct-by-threat-by-time interactions
  # start w most complex and drop down if necessary
  t.poly.max.sig <- max(substr(cfg.gold$fixed[grepl('Time',cfg.gold$fixed)],5,5)) # only test interactions up to the maximum time you've selected
  if(!is.na(t.poly.max.sig)){
    test.time.by.c.by.t.cfg <- lapply(t.poly.max.sig:1, function(o.max) list(fixed=c(cfg.gold$fixed,paste0('Time',1:o.max,'*Correct*Threat')), random=c(cfg.gold$random),response=cfg.gold$response,id=cfg.gold$id)) # define sequential list of models to iterate through
    list[mdl.gold,cfg.gold] <- lme.stepdown(mdl.gold,test.time.by.c.by.t.cfg)
    # test 2-way interactions for correct-by-time and threat-by-time for time terms with out 3-ways
    t.c.t.max.sig <- max(substr(cfg.gold$fixed[grepl('Time[0-9]\\*Correct\\*Threat',cfg.gold$fixed,perl = T)],5,5)) # only add 2-way interactions for times with no 3-way interactions
    if(is.na(t.c.t.max.sig)){t.c.t.max.sig <- 1} # if no 3-ways, then test all possible two ways
    if(t.c.t.max.sig < t.poly.max.sig){ # if there are any time terms without 3-way interactions with correct or threat (NA means no 3-way interactions)
      test.time.by.c.cfg <- lapply(t.poly.max.sig:t.c.t.max.sig, function(o.max) list(fixed=c(cfg.gold$fixed,paste0('Time',1:o.max,'*Correct')), random=c(cfg.gold$random),response=cfg.gold$response,id=cfg.gold$id)) # define sequential list of models to iterate through
      test.time.by.t.cfg <- lapply(t.poly.max.sig:t.c.t.max.sig, function(o.max) list(fixed=c(cfg.gold$fixed,paste0('Time',1:o.max,'*Threat')), random=c(cfg.gold$random),response=cfg.gold$response,id=cfg.gold$id)) # define sequential list of models to iterate through
      list[mdl.gold,cfg.gold] <- lme.stepdown(mdl.gold,test.time.by.c.cfg)
      list[mdl.gold,cfg.gold] <- lme.stepdown(mdl.gold,test.time.by.t.cfg)
    } 
    
    #t.poly.max.sig <- 0
    if(any(grepl('Time',cfg.gold$fixed))){
      # 4. now add random effects of increasing time order
      t.poly.max.sig <- max(substr(cfg.gold$fixed[grepl('Time',cfg.gold$fixed)],5,5)) # only add random effects up to the maximum time you've selected
      test.time.ran.cfg <- lapply(1:t.poly.max.sig, function(o.max) list(fixed=cfg.gold$fixed, random=c(cfg.gold$random,paste0('Time',1:o.max)),response=cfg.gold$response,id=cfg.gold$id)) # define sequential list of models to iterate through
      list[mdl.gold,cfg.gold] <- lme.stepup(mdl.gold,test.time.ran.cfg)
      
      # 5. now test additions of fixed effect interactions with Is22q and Time
      test.time.by.22q.cfg <- lapply(1:t.poly.max.sig, function(o.max) list(fixed=c(cfg.gold$fixed,paste0('Is22q*Time',1:o.max)), random=c(cfg.gold$random),response=cfg.gold$response,id=cfg.gold$id)) # define sequential list of models to iterate through
      list[mdl.gold,cfg.gold] <- lme.selectbest(mdl.gold,cfg.gold,test.time.by.22q.cfg)
      # test.time.by.22q.cfg <- lapply(t.poly.max.sig:1, function(o.max) list(fixed=c(cfg.gold$fixed,paste0('Is22q*Time',1:o.max)), random=c(cfg.gold$random),response=cfg.gold$response,id=cfg.gold$id)) # define sequential list of models to iterate through
      # list[mdl.gold,cfg.gold] <- lme.stepdown(mdl.gold,test.time.by.22q.cfg)
      
      # 6. now add four way interactions between 22q*time*correct*threat, using step-down approach
      # this monstrous term asks whether effects of 22q on time differ by correct responses/threat responses
      test.time.by.22q.by.c.by.t.cfg <- lapply(t.poly.max.sig:1, function(o.max) list(fixed=c(cfg.gold$fixed,paste0('Is22q*Time',1:o.max,'*Correct*Threat')), random=c(cfg.gold$random),response=cfg.gold$response,id=cfg.gold$id)) # define sequential list of models to iterate through
      list[mdl.gold,cfg.gold] <- lme.stepdown(mdl.gold,test.time.by.22q.by.c.by.t.cfg)
      # if no 4-ways, test 3-ways just with 22q*Time*threat or 22q*Time*correct, and stop there
      t.c.t.22q.max.sig <- max(substr(cfg.gold$fixed[grepl('Is22q\\*Time[0-9]\\*Correct\\*Threat',cfg.gold$fixed,perl = T)],5,5)) # only add 3-way interactions for times with no 2-way interactions
      if(is.na(t.c.t.22q.max.sig)){t.c.t.22q.max.sig <- 1} # if no 4-ways, then test all possible 3 ways
      if(t.c.t.22q.max.sig < t.poly.max.sig){ # if there are any time terms without 4-way interactions with correct or threat (NA means no 4-way interactions)
        test.time.by.c.by.22q.cfg <- lapply(t.poly.max.sig:t.c.t.22q.max.sig, function(o.max) list(fixed=c(cfg.gold$fixed,paste0('Is22q*Time',1:o.max,'*Correct')), random=c(cfg.gold$random),response=cfg.gold$response,id=cfg.gold$id)) # define sequential list of models to iterate through
        test.time.by.t.by.22q.cfg <- lapply(t.poly.max.sig:t.c.t.22q.max.sig, function(o.max) list(fixed=c(cfg.gold$fixed,paste0('Is22q*Time',1:o.max,'*Threat')), random=c(cfg.gold$random),response=cfg.gold$response,id=cfg.gold$id)) # define sequential list of models to iterate through
        list[mdl.gold,cfg.gold] <- lme.stepdown(mdl.gold,test.time.by.c.by.22q.cfg)
        list[mdl.gold,cfg.gold] <- lme.stepdown(mdl.gold,test.time.by.t.by.22q.cfg)
      } 
    }
  }
    # if(t.poly.max.sig == 1){ # if time is still linear then try to add in time again -- happens for PC 6
    #   test.time.by.c.by.t.cfg <- lapply(t.poly.max:(as.numeric(t.poly.max.sig)+1), function(o.max) list(fixed=c(cfg.gold$fixed,paste0('Correct*Threat*Time',1:o.max)), random=c(cfg.gold$random),response=cfg.gold$response,id=cfg.gold$id)) # define sequential list of models to iterate through
    #   list[mdl.gold,cfg.gold] <- lme.stepdown(mdl.gold,test.time.by.c.by.t.cfg)
    #   t.c.t.max.sig <- max(substr(cfg.gold$fixed[grepl('Correct\\*Threat\\*Time[0-9]',cfg.gold$fixed,perl = T)],20,20)) # only add 3-way interactions for times with no 2-way interactions
    #   # add random effects back in if possible
    #   test.time.ran.cfg <- lapply(1:t.c.t.max.sig, function(o.max) list(fixed=cfg.gold$fixed, random=c(cfg.gold$random,paste0('Time',1:o.max)),response=cfg.gold$response,id=cfg.gold$id)) # define sequential list of models to iterate through
    #   list[mdl.gold,cfg.gold] <- lme.stepup(mdl.gold,test.time.ran.cfg)
    #   # now try to add in interaction with 22q
    #   test.time.by.22q.by.c.by.t.cfg <- lapply(t.c.t.max.sig:1, function(o.max) list(fixed=c(cfg.gold$fixed,paste0('Correct*Threat*Is22q*Time',1:o.max)), random=c(cfg.gold$random),response=cfg.gold$response,id=cfg.gold$id)) # define sequential list of models to iterate through
    #   list[mdl.gold,cfg.gold] <- lme.stepdown(mdl.gold,test.time.by.22q.by.c.by.t.cfg)
    # }
    
    #working.mdl <- lme(fixed=Score~Time,random=~1|scanid,data=df.all,na.action = na.exclude) # to test this tryCatch
    ranef.test <- tryCatch(intervals(mdl.gold),error=function(err){return(FALSE)}) # check if data is singular with random effects after adding interactions
    if(!is.list(ranef.test)){
      while(!is.list(ranef.test)){ # delete last random effect and refit until random effect var-cov matrix can be computed
        cfg.gold$random <- cfg.gold$random[-length(cfg.gold$random)] 
        mdl.gold <- lme.ms(cfg.gold,df.all)
        ranef.test <- tryCatch(intervals(mdl.gold),error=function(err){return(FALSE)})
      }
    }
    results.PC$mdl.best.max.tsig <- t.poly.max.sig
    print(paste0('PC',PC,', best model maximum time order: ',results.PC$mdl.best.max.tsig))
    results.PC$mdl.best <- mdl.gold
    results.PC$cfg.best <- cfg.gold
    
    # get predicted values from best model
    results.PC$mdl.best$data$pred.best <- predict(results.PC$mdl.best)
    
    # get stats for interactions, i.e. group differences, from best model
    results.PC$mdl.best.stats <- summary(results.PC$mdl.best)$tTable
    # time-by-22q interactions
    coef.names <- rownames(results.PC$mdl.best.stats)
    find.ints <- grep('Is22q:Time*',coef.names)
    results.PC$p.ints.22q <- results.PC$mdl.best.stats[find.ints,'p-value',drop=FALSE] # p-values for interactions
    # time-by-correct/threat interactions, excluding 22q
    find.ints <- which(grepl('Time[0-9]:',coef.names) & !grepl('Is22q',coef.names))
    results.PC$p.ints.correct.threat <- results.PC$mdl.best.stats[find.ints,'p-value',drop=FALSE] # p-values for interactions
    
    results[[PC]] <- results.PC
}

sapply(results, function(X) X$p.ints.22q)
sapply(results, function(X) X$mdl.best.stats[rownames(X$p.ints.correct.threat),])

save(results, file = paste0(savedir,'lme_all_trials/',component_design,'_PCTimeCoursesLMEEffectsCorrectThreat.RData'))
