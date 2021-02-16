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

masterdir <- paste0(basedir,'results/',name_root,'/')
savedir <- paste0(masterdir,'analyses/sc_vs_pc/')

dir.create(savedir,recursive = TRUE)

source(paste0(basedir,'code/miscfxns/packages.R'))
source(paste0(basedir,'code/plottingfxns/plottingfxns.R'))
source(paste0(basedir,'code/miscfxns/miscfxns.R'))

TR <- 3 # repetition time in seconds
nTR <- 204 # number of TRs in scan
demo <- read.csv(paste(basedir,'data/Demographics',name_root,'.csv',sep=''),stringsAsFactors=F)
rownames(demo) <- as.character(demo$scanid)
static.vars <- readMat(paste0(basedir,'data/Demographics',name_root,'.mat'))
Hippocampus <- 209:210
Amygdala <- 211:212 # index of amygdala and hippocampus
Accumbens <- 213:214
ROI.Structure <- Amygdala
covariates <- c('scanage_months','sex','BrainSegVol','idemo_meanrelrms','handedness','dti64MeanRelRMS') # specify covariates


#for(sc_data_name in c('QA_Pass','GFA_Pass','SIFT_radius2')){
sc_data_name <- 'GFA_Pass'
  matFile.SC <- readMat(paste0(basedir,'data/sc/StructuralConnectivity',sc_data_name,name_root,'.mat'))
  
  # difference in weighted degree
  
  SC <- matFile.SC$SC
  dti.missing <- matFile.SC$dti.missing
  ids <- matFile.SC$SC.IDs
  if(!identical(as.character(ids),as.character(demo$bblid))){print('ERROR: IDs do not match')}
  ROI.Strength <- do.call(rbind,setNames(lapply(1:dim(SC)[3],function(N) mean(SC[ROI.Structure,,N])),demo$scanid))


# load BOLD data
load.dir <- paste0(masterdir,'analyses/fir/cpc_timecourse_fin',fin,'st',st,'/',component_design,'/')
stim.types <- list(threat=1,nonthreat=0)
response.types <- list(correct=1,incorrect=0)
load(file = paste0(load.dir,'lme_all_trials/',component_design,'_PCTimeCoursesLMEEffectsCorrectThreat.RData'))
if(grepl('xcp_36p_despike',name_root)){
  ncomps <- 6
  results <- results[1:ncomps]
} else if(grepl('xcp_6p_noFilter',name_root)){
  ncomps <- 5
  results <- results[2:6] # remove comp 1 which is global signal
}
PC.summary <- list()
for(stim.type in names(stim.types)){
  for(response.type in names(response.types)){
    event.name <- paste0(stim.type,response.type) 
    stim.type.idx <- stim.types[[stim.type]]
    response.type.idx <- response.types[[response.type]]
    # load each response type
    fun <- function(x) max(x)-min(x)
    df.list <- lapply(results, function(X) X$mdl.best$data[complete.cases(X$mdl.best$data) & X$mdl.best$data$Threat == stim.type.idx & X$mdl.best$data$Correct == response.type.idx,])
    PC.summary[[event.name]] <- lapply(df.list, function(df) aggregate(df,by=list(scanid=df$scanid),fun)[,c('scanid','Score'),drop=F])
    PC.summary[[event.name]] <- do.call(cbind,lapply(PC.summary[[event.name]], function(X) data.frame(X[setdiff(colnames(X),'scanid')],row.names=X$scanid)))
    colnames(PC.summary[[event.name]]) <- paste0('PC',1:ncomps)
    PC.summary[[event.name]]$Event <- event.name
    PC.summary[[event.name]] <- cbind(PC.summary[[event.name]],demo[rownames(PC.summary[[event.name]]),covariates]) # get covariates associated with each subject's PC peak
    PC.summary[[event.name]]$ROIStrength <- ROI.Strength[rownames(PC.summary[[event.name]]),] # get amygdala strength asasociated with each subject
  }
}

m.res <- list()
p.mat <- matrix(NA,nrow=length(PC.summary),ncol=ncomps,dimnames = list(names(PC.summary),1:ncomps))
for(event.name in names(PC.summary)){
  m.res[[event.name]] <- list()
  for(PC in 1:ncomps){
    PC.summary[[event.name]]$y <- PC.summary[[event.name]][,paste0('PC',PC)]
    m.res[[event.name]][[paste0('PC',PC)]] <- summary(lm.beta(lm(y~scanage_months+sex+handedness+BrainSegVol +dti64MeanRelRMS+idemo_meanrelrms+ROIStrength,data=PC.summary[[event.name]])))$coef['ROIStrength',]
    p.mat[event.name,PC] <- m.res[[event.name]][[paste0('PC',PC)]]['Pr(>|t|)']
    #m.res[[event.name]][[paste0('PC',PC)]] <- unlist(cor.test(PC.summary[[event.name]]$y,PC.summary[[event.name]]$ROIStrength,use='pairwise.complete.obs',method = 'spearman')[c('estimate','p.value')])
  }
}
p.mat
plot(PC.summary[[event.name]]$ROIStrength,PC.summary[[event.name]]$PC5)
