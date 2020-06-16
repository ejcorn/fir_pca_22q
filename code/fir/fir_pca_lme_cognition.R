# compare time parameters of each PC to social cognition and executive function

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
#demo.cog <- merge(demo,idemo.performance,by='scanid') # use this instead of the above 5 lines in order to use all subjects
rownames(demo.cog) <- demo.cog$scanid
# join and merge diagnoses, social cognition, execeff, etc.

grp.colors <- getGroupColors()

# load FIR betas
component_design <- 'ThreatNonthreatAllStimuliStratified'
savedir <- paste0(masterdir,'analyses/fir/subject_fir_correct_incorrect_pca/cpc_timecourse/',component_design,'/')
dir.create(paste0(savedir,'lme'),recursive=T)
stim.types <- list(all='all',threat='threat',nonthreat='nonthreat')
response.types <- list(correct='correct',incorrect='incorrect')

stim.type <- stim.types[[1]]
response.type <- response.types[[1]]
results <- list()
load(file = paste0(savedir,'lme/',component_design,'_PCTimeCoursesLMEEffects.RData'))
savedir <- paste0(savedir,'lme/predict_cog_22q/')
dir.create(savedir,recursive = T)

#covariates <- c(covariates,'exe_eff','cpxres_eff') # add in some general cognition covariates
cog.vars <- c('allcorrect','soccog_eff')
mdls <- list(lm,lm,lm); names(mdls) <- cog.vars

results.cog <- list() # store all of the R^2 values and p-values in one matrix and post-hoc correct then plot
for(stim.type in c('all')){
  for(response.type in response.types){
    # load each response type
    results.bold <- results[[stim.type]][[response.type]]
    ncomps <- length(results.bold)
    results.dnames <- list(paste0('PC',1:length(results.bold)),cog.vars)
    # make results matrix:
    res.mat <- matrix(NA,nrow = length(results.bold),ncol=length(cog.vars),dimnames=results.dnames)
    results.cog[[paste0(stim.type,response.type)]]$betas <- results.cog[[paste0(stim.type,response.type)]]$p.values <- res.mat
    results.cog[[paste0(stim.type,response.type)]]$rsq <- results.cog[[paste0(stim.type,response.type)]]$f.p.values <- res.mat
    
    # extract some kind of summary metric of each PC's time course
    #PC.summary <- lapply(results.bold, function(X) as.data.frame(ranef(X$mdl.best))) # extract random effects
    df.list <- lapply(results.bold, function(X) X$mdl.best$data[complete.cases(X$mdl.best$data),])
    PC.summary <- lapply(df.list, function(df) aggregate(df,by=list(scanid=df$scanid),mean)[,c('scanid','Score'),drop=F])
    PC.summary <- lapply(PC.summary, function(X) data.frame(X[setdiff(colnames(X),'scanid')],row.names=X$scanid))
    #compare random effects within 22q to social cognitive deficits
    for(y in cog.vars){
      # concatenate outcome variable, scores (predictor of interest), and specified covariates
      df.PCs <- lapply(1:ncomps, function(PC) cbind(y=demo.cog[,y],PC.summary[[PC]][as.character(demo.cog$scanid),,drop=F],demo.cog[,covariates]))
      # exclude people with any NAs
      df.PCs <- lapply(df.PCs, function(X) X[complete.cases(X),])
      # exclude outlier wrt outcome variable
      #df.PCs <- lapply(df.PCs, function(X) X[!outlier.mask(X$y),])
      # replace variables with rank
      #for(PC in 1:ncomps){df.PCs[[PC]]$y <- rank(df.PCs[[PC]]$y)}
      # fit one model using all variables
      lm.PCs <- lapply(df.PCs,function(X) lm.beta(mdls[[y]](y~.,data=X)))
      # first regress out covariates -- shouldn't be necessary to fit this on every PC, but rarely missingness can differ by PCs
      resid.PCs <- lapply(df.PCs,function(X) residuals(mdls[[y]](formula=reformulate(termlabels = covariates,response='y',intercept = T),data=X)))
      # now add residuals to df
      df.PCs <- lapply(1:length(df.PCs), function(PC) cbind(df.PCs[[PC]],y.r=resid.PCs[[PC]]))
      # remove covariates and original dependent
      df.PCs <- lapply(df.PCs, function(X) X[,setdiff(colnames(X),c('y',covariates))])
      # now fit a model to predict residuals from remaining variables (either score or random effects)
      lm.PCs <- lapply(df.PCs,function(X) lm.beta(mdls[[y]](y.r~.,data=X)))

      coef.list <- lapply(lm.PCs,get.coef)
      results.cog[[paste0(stim.type,response.type)]][[y]] <- list()
      results.cog[[paste0(stim.type,response.type)]][[y]]$models <- lm.PCs
      results.cog[[paste0(stim.type,response.type)]][[y]]$coef.tables <- coef.list
      results.cog[[paste0(stim.type,response.type)]]$rsq[,y] <- sapply(lm.PCs,get.rsq)
      results.cog[[paste0(stim.type,response.type)]]$f.p.values[,y] <- sapply(lm.PCs,get.ftest.pval)
      results.cog[[paste0(stim.type,response.type)]]$betas[,y] <- sapply(coef.list,function(X) X['Score','Standardized'])
      results.cog[[paste0(stim.type,response.type)]]$p.values[,y] <- sapply(coef.list,function(X) X['Score','Pr(>|t|)'])
      #p.list <- lapply(lm.PCs, function(m) p.xy.flex(x=m$fitted.values,y=m$model$y,xlab = 'Fitted',ylab=y))
      p.list <- lapply(lm.PCs, function(m) p.xy.flex(x=get.partial.resids(m,'Score')$x,y=m$model$y,xlab = 'Score',ylab=y))
      p.all <- plot_grid(plotlist = p.list,nrow=1,align='hv')
      ggsave(plot = p.all, filename = paste0(savedir,stim.type,response.type,'Score_Predict',y,'.pdf'),
             width = 18,height = 4, units = "cm",useDingbats=F)
    }
    
  }
}

pval.list <- lapply(setNames(names(results.cog),names(results.cog)), function(X) list(f.p.values=results.cog[[X]]$f.p.values))
pval.list <- list.posthoc.correct(pval.list,'fdr')

p.list <- lapply(names(results.cog), function(X) imagesc(results.cog[[X]]$rsq *(p.adjust(pval.list[[X]]$f.p.values,'fdr') < 0.05),
                                                         cmap = 'Blues',ttl = X,caxis_name = expression(R^2)) + 
                   theme(axis.text.x = element_text(size=8,angle=90,hjust=1,vjust=0.5),
                         legend.position = 'bottom',legend.key.height = unit(0.1,'cm'),legend.key.width=unit(0.4,'cm'))) 
p.all <- plot_grid(plotlist = p.list,nrow=1,align='hv')
ggsave(plot = p.all, filename = paste0(savedir,'MeanScoreRsquare.pdf'),
       width = 7,height = 6, units = "cm",useDingbats=F)

pval.list <- lapply(setNames(names(results.cog),names(results.cog)), function(X) list(p.values=results.cog[[X]]$p.values))
pval.list <- list.posthoc.correct(pval.list,'fdr')
p.list <- lapply(names(results.cog), function(X) imagesc(results.cog[[X]]$betas *(p.adjust(pval.list[[X]]$p.values,'fdr') < 0.05),
          cmap = 'redblue_asymmetric',ttl = X,caxis_name = expression(beta)) + 
            theme(axis.text.x = element_text(size=8,angle=90,hjust=1,vjust=0.5),
                  legend.position = 'bottom',legend.key.height = unit(0.1,'cm'),legend.key.width=unit(0.4,'cm'))) 
p.all <- plot_grid(plotlist = p.list,nrow=1,align='hv')
ggsave(plot = p.all, filename = paste0(savedir,'MeanScoreBetas.pdf'),
       width = 7,height = 6, units = "cm",useDingbats=F)

