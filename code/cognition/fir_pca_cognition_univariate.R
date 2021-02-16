# compare time parameters of each PC to social cognition and executive function

rm(list=ls())

args <- commandArgs(TRUE)
name_root <- args[1]
basedir <- args[2]
component_design <- args[3]
fin <- 6
st <- 1

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
savedir <- paste0(masterdir,'analyses/fir/cpc_timecourse_fin',fin,'st',st,'/',component_design,'/')
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
cog.vars <- c('allcorrect','soccog_eff')
mdls <- list(lm,glm); names(mdls) <- cog.vars

results.cog <- PC.summary <- list() # store all of the R^2 values and p-values in one matrix and post-hoc correct then plot
for(stim.type in names(stim.types)){
  for(response.type in names(response.types)){
    res.name <- paste0(stim.type,response.type) 
    stim.type.idx <- stim.types[[stim.type]]
    response.type.idx <- response.types[[response.type]]
    # load each response type
    fun <- function(x) max(abs(x))
    df.list <- lapply(results, function(X) X$mdl.best$data[complete.cases(X$mdl.best$data) & X$mdl.best$data$Threat == stim.type.idx & X$mdl.best$data$Correct == response.type.idx,])
    PC.summary[[res.name]] <- lapply(df.list, function(df) aggregate(df,by=list(scanid=df$scanid),fun)[,c('scanid','Score'),drop=F])
    PC.summary[[res.name]] <- do.call(cbind,lapply(PC.summary[[res.name]], function(X) data.frame(X[setdiff(colnames(X),'scanid')],row.names=X$scanid)))
    colnames(PC.summary[[res.name]]) <- paste0('PC',1:ncomps,res.name)
  }
}

y <- 'allcorrect'
df.resid <- demo.cog
df.resid$y <- df.resid[,y]
if(length(unique(df.resid$y))>2){df.resid$y <- rank(df.resid$y)} # if y is continuous convert to rank
# exclude outlier wrt outcome variable, if it's a continuous var:
# if(length(unique(df.resid$y))>2){
#   df.resid <- df.resid[!outlier.mask(df.resid$y),]
# }
# residualize outcome variable
y.r <- residuals(lm(formula=reformulate(termlabels = covariates,response='y',intercept = T),data=df.resid)) 
if(length(unique(df.resid$y))>2){y.r <- rank(y.r)}  # if original y was continuous convert residuals to rank
r.mat <- p.mat <- matrix(NA,nrow=length(PC.summary),ncol=ncomps,dimnames = list(names(PC.summary),1:ncomps))
vals <- list()
for(res.name in names(PC.summary)){
  vals[[res.name]] <- list()
  # concatenate outcome variable, scores (predictor of interest), and specified covariates
  df.pred <- cbind(y.r=y.r,PC.summary[[res.name]][names(y.r),])
  for(PC in 1:ncomps){
    c.test <- cor.test(df.pred[,paste0('PC',PC,res.name)],df.pred$y.r)
    r.mat[res.name,PC] <- c.test$estimate
    p.mat[res.name,PC] <- c.test$p.value
    vals[[res.name]][[PC]] <- list(x=df.pred[,paste0('PC',PC,res.name)],y=df.pred$y.r)
  }
}
p.mat
p.mat <- list.posthoc.correct(list(p.mat),method = 'fdr')[[1]]
p.mat
plot(vals$nonthreatincorrect[[2]]$x,vals$nonthreatincorrect[[2]]$y)