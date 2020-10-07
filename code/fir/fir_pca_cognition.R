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
# demo <- demo[demo$study=='22q',]
# demo.cog <- merge(demo2,demo,by='bblid')
# ind <- order(match(demo.cog$bblid,demo$bblid)) # need to do this b.c
# demo.cog <- demo.cog[ind,]
# demo.cog <- cbind(demo.cog,idemo.performance[as.character(demo.cog$scanid),])
demo.cog <- merge(demo,idemo.performance,by='scanid') # use this instead of the above 5 lines in order to use all subjects
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
mdls <- list(lm,glm); names(mdls) <- cog.vars

results.cog <- list() # store all of the R^2 values and p-values in one matrix and post-hoc correct then plot
for(stim.type in names(stim.types)){
  for(response.type in names(response.types)){
    res.name <- paste0(stim.type,response.type) 
    stim.type.idx <- stim.types[[stim.type]]
    response.type.idx <- response.types[[response.type]]
    results.cog[[res.name]] <-list()
    # load each response type
    results.bold <- results
    ncomps <- length(results.bold)
    fun <- function(x) max(abs(x))
    df.list <- lapply(results.bold, function(X) X$mdl.best$data[complete.cases(X$mdl.best$data) & X$mdl.best$data$Threat == stim.type.idx & X$mdl.best$data$Correct == response.type.idx,])
    PC.summary <- lapply(df.list, function(df) aggregate(df,by=list(scanid=df$scanid),fun)[,c('scanid','Score'),drop=F])
    PC.summary <- do.call(cbind,lapply(PC.summary, function(X) data.frame(X[setdiff(colnames(X),'scanid')],row.names=X$scanid)))
    colnames(PC.summary) <- paste0('PC',1:ncomps)
    #PC.ranefs <- lapply(results.bold, function(X) as.data.frame(ranef(X$mdl.best)))
    #PC.ranefs <- lapply(1:length(PC.ranefs),function(PC) setNames(PC.ranefs[[PC]],paste0(names(PC.ranefs[[PC]]),'PC',PC)))
    #PC.summary <- do.call(cbind,PC.ranefs)
    
    for(y in cog.vars){
      # concatenate outcome variable, scores (predictor of interest), and specified covariates
      df.PCs <- cbind(y=demo.cog[,y],PC.summary[as.character(demo.cog$scanid),],demo.cog[,covariates])
      # exclude people with any NAs
      df.PCs <- df.PCs[complete.cases(df.PCs),]
      # exclude outlier wrt outcome variable, if it's a continuous var:
      if(length(unique(df.PCs$y))>2){
        df.PCs <- df.PCs[!outlier.mask(df.PCs$y),]
      }
      resid.PCs <- residuals(mdls[[y]](formula=reformulate(termlabels = covariates,response='y',intercept = T),data=df.PCs))
      # now add residuals to df
      df.PCs <- cbind(df.PCs,y.r=resid.PCs)
      # remove covariates and original dependent
      df.PCs <- df.PCs[,setdiff(colnames(df.PCs),c('y',covariates))]
      # now fit a model to predict residuals from remaining variables (either score or random effects)
      lm.PCs <- lm.beta(mdls[[y]](y.r~.,data=df.PCs))
      results.cog[[res.name]][[y]]$yhat <- fitted(lm.PCs)
      results.cog[[res.name]][[y]]$y <- lm.PCs$model$y.r
      results.cog[[res.name]][[y]]$f <- summary(lm.PCs)$fstatistic
      results.cog[[res.name]][[y]]$p <- get.ftest.pval(lm.PCs)
      results.cog[[res.name]][[y]]$rsq <- get.rsq(lm.PCs)
      results.cog[[res.name]][[y]]$coefs <- get.coef(lm.PCs)
      results.cog[[res.name]][[y]]$mdl <- lm.PCs
    }
  }
}

iter <- setNames(names(results.cog),names(results.cog))
pval.list.un <- lapply(iter, function(X) lapply(setNames(cog.vars,cog.vars), function(y) list(p=results.cog[[X]][[y]]$p)))
pval.list <- list.posthoc.correct(pval.list.un,'fdr')

f.list <- lapply(iter, function(X) lapply(setNames(cog.vars,cog.vars), function(y) list(f=results.cog[[X]][[y]]$f)))

p.list <- list()
col <- relist(flesh=brewer.pal(10,'BuGn')[-c(1:3,8:9)],skeleton = pval.list)
col <- relist(flesh=rep('#c85795',nrow(expand.grid(iter,cog.vars))),skeleton=pval.list)
ylab.pretty <- c(soccog_eff='Social Cognition Factor',allcorrect='Emotion ID Correct (%)')
for(p in iter){
    for(y in cog.vars){
      df.plt <- data.frame(y=results.cog[[p]][[y]]$y,yhat=results.cog[[p]][[y]]$yhat)
      rsq <- signif(results.cog[[p]][[y]]$rsq,2)
      r.lab <- paste0('R^2 == ',rsq)
      p.list[[p]][[y]] <- ggplot(df.plt) + geom_point(aes(x=yhat,y=y),color=col[[p]][[y]],stroke=0,alpha=0.5) + 
        geom_smooth(aes(x=yhat,y=y),fill=col[[p]][[y]],color=col[[p]][[y]],method='lm') + 
        xlab(paste('PC Peaks,',p,'\n(Fitted)')) + ylab(ylab.pretty[y]) +
        theme_classic() + theme(text = element_text(size = 6)) + 
        theme(plot.title = element_text(size=8,hjust=0.5,face = "bold")) +
        theme(plot.margin=grid::unit(c(0,0,0,0), "mm")) + theme(legend.position = 'none') +
        annotate(geom='text',x=Inf,y=-Inf,label=paste0('p[FDR] == ',signif(pval.list[[p]][[y]]$p,2)),parse=T,vjust=-1,hjust=1,size=2)+
      annotate(geom='text',x=Inf,y=-Inf,label=r.lab,parse=T,vjust=-2,hjust=1,size=2)
    }
}
p.all <- plot_grid(plotlist = Reduce(c,p.list),ncol=4)
ggsave(p.all,filename = paste0(savedir,'PredictCognitiveVariablesFromPCScores.pdf'),
       units = 'cm',height = 9,width = 18)
