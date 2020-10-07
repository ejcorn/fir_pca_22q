# compare time parameters of each PC to social cognition and executive function

rm(list=ls())

args <- commandArgs(TRUE)
name_root <- args[1]
basedir <- args[2]
component_design <- args[3]

basedir <- '~/Dropbox/Cornblath_Bassett_Projects/BrainStates22q/fir_pca_22q/'
name_root <- 'CPCA_IDSchaefer200Z1xcp_6p_noFilter'
#basedir <- '/cbica/home/cornblae/ecornblath/fir_pca_22q/'
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
covariates <- c('scanage_months','sex','BrainSegVol','idemo_meanrelrms','handedness') # specify covariates

demo2 <- read.csv('data/n85_22q_deleted_demo_dx_cnb.csv',stringsAsFactors=F)
demo2 <- demo2[,setdiff(colnames(demo2),covariates)] # remove covariates from this df - only use one from processed file
idemo.performance <- read.csv(paste0(masterdir,'analyses/behavior/idemo/IDEmoAccuracy.csv'),stringsAsFactors = F,row.names = 1)
# only study 22q subjects - want to explain variation within this phenotype
demo <- demo[demo$study=='22q',]
demo.cog <- merge(demo2,demo,by='bblid')
ind <- order(match(demo.cog$bblid,demo$bblid)) # need to do this b.c
demo.cog <- demo.cog[ind,]
demo.cog <- cbind(demo.cog,idemo.performance[as.character(demo.cog$scanid),])
# demo.cog <- merge(demo,idemo.performance,by='scanid') # use this instead of the above 5 lines in order to use all subjects
rownames(demo.cog) <- demo.cog$scanid
# join and merge diagnoses, social cognition, execeff, etc.

grp.colors <- getGroupColors()

# load FIR betas
savedir <- paste0(masterdir,'analyses/fir/subject_fir_correct_incorrect_pca/cpc_timecourse/',component_design,'/')
stim.types <- list(threat=1,nonthreat=0)
response.types <- list(correct=1,incorrect=0)

stim.type <- stim.types[[1]]
response.type <- response.types[[1]]
load(file = paste0(savedir,'lme_all_trials/',component_design,'_PCTimeCoursesLMEEffectsCorrectThreat.RData'))

if(grepl('xcp_36p_despike',name_root)){
  ncomps <- 6
  results <- results[1:ncomps]
} else if(grepl('xcp_6p_noFilter',name_root)){
  ncomps <- 5
  results <- results[2:6] # remove comp 1 which is global signal
}

savedir <- paste0(savedir,'predict_cog_22q/')
dir.create(savedir,recursive = T)

#covariates <- c(covariates,'exe_eff','cpxres_eff') # add in some general cognition covariates
cog.vars <- c('allcorrect')
mdls <- list(lm); names(mdls) <- cog.vars

results.cog <- PC.summary <- list() # store all of the R^2 values and p-values in one matrix and post-hoc correct then plot
for(stim.type in names(stim.types)){
  for(response.type in names(response.types)){
    res.name <- paste0(stim.type,response.type) 
    stim.type.idx <- stim.types[[stim.type]]
    response.type.idx <- response.types[[response.type]]
    # load each response type
    results.bold <- results
    ncomps <- length(results.bold)
    fun <- function(x) max(abs(x))
    df.list <- lapply(results.bold, function(X) X$mdl.best$data[complete.cases(X$mdl.best$data) & X$mdl.best$data$Threat == stim.type.idx & X$mdl.best$data$Correct == response.type.idx,])
    PC.summary[[res.name]] <- lapply(df.list, function(df) aggregate(df,by=list(scanid=df$scanid),fun)[,c('scanid','Score'),drop=F])
    PC.summary[[res.name]] <- do.call(cbind,lapply(PC.summary[[res.name]], function(X) data.frame(X[setdiff(colnames(X),'scanid')],row.names=X$scanid)))
    colnames(PC.summary[[res.name]]) <- paste0('PC',1:ncomps,res.name)
  }
}
# make one big data frame to predict outcomes
df.all <- demo.cog[,c(covariates)]
rownames(df.all) <- demo.cog$scanid
for(PC.df in PC.summary){
  df.all <- cbind(df.all,PC.df[rownames(df.all),])
}
for(y in cog.vars){
  # concatenate outcome variable, scores (predictor of interest), and specified covariates
  df.pred <- cbind(y=demo.cog[,y],df.all)
  # exclude people with any NAs
  df.pred <- df.pred[complete.cases(df.pred),]
  # exclude outlier wrt outcome variable, if it's a continuous var:
  if(length(unique(df.pred$y))>2){
    df.pred <- df.pred[!outlier.mask(df.pred$y),]
  }
  resid.PCs <- residuals(mdls[[y]](formula=reformulate(termlabels = covariates,response='y',intercept = T),data=df.pred))
  # now add residuals to df
  df.pred <- cbind(df.pred,y.r=resid.PCs)
  # remove covariates and original dependent
  df.pred <- df.pred[,setdiff(colnames(df.pred),c('y',covariates))]
  # now fit a model to predict residuals from remaining variables (either score or random effects)
  lm.PCs <- lm.beta(mdls[[y]](y.r~.,data=df.pred))
  results.cog[[y]]$yhat <- fitted(lm.PCs)
  results.cog[[y]]$y <- lm.PCs$model$y.r
  results.cog[[y]]$f <- summary(lm.PCs)$fstatistic
  results.cog[[y]]$p <- get.ftest.pval(lm.PCs)
  results.cog[[y]]$rsq <- get.rsq(lm.PCs)
  results.cog[[y]]$coefs <- get.coef(lm.PCs)
  results.cog[[y]]$mdl <- lm.PCs
  results.cog[[y]]$df <- df.pred
}
PC.features <- setdiff(colnames(df.pred),'y.r')
names(PC.features) <- PC.features
c.tests <- lapply(PC.features,function(PC) unlist(cor.test(df.pred[,PC],df.pred$y.r)[c('estimate','p.value')]))
p.adjust(sapply(c.tests,function(X) X['p.value']),method = 'fdr')


