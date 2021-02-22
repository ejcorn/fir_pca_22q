# compare time parameters of each PC to social cognition and executive function

rm(list=ls())

args <- commandArgs(TRUE)
name_root <- args[1]
basedir <- args[2]
component_design <- args[3]
fin <- 6
st <- 0

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

# define keys for plot labels
nice.x.name.key <- c(threatcorrect='Threat Correct',threatincorrect='Threat Incorrect',nonthreatcorrect='Non-Threat Correct',nonthreatincorrect='Non-Threat Incorrect')
nice.y.name.key <- c(allcorrect='Emotion ID\nAccuracy (Rank)',overall_speed='Overall Speed')
nice.y.name.key.hist <- c(allcorrect = 'Emotion ID\nAccuracy (%)',overall_speed='Overall Speed')

subj.samples <- c('22q','AllSubjects','PNC')
subj.samples <- c('22q','PNC')
subj.sample <- '22q'
cog.vars <- c('allcorrect')
overlap.hist <- lapply(cog.vars,function(X) list()) # make histograms of performance for 22q and PNC samples used in this analysis
for(subj.sample in subj.samples){ # do this both within 22q sample only and within the full sample and within controls only
  demo <- read.csv(paste(basedir,'data/Demographics',name_root,'.csv',sep=''),stringsAsFactors=F)
  #cnb <- read.csv(paste0(basedir,'data/CNBMeds_',name_root,'.csv'),stringsAsFactors = F)
  #demo <- merge(demo,cnb[,c('scanid','cnb')],by='scanid')
  rownames(demo) <- as.character(demo$scanid)
  covariates <- c('scanage_months','sex','BrainSegVol','idemo_meanrelrms','handedness') # specify covariates
  demo2 <- read.csv('data/n85_22q_deleted_demo_dx_cnb.csv',stringsAsFactors=F)
  demo2 <- demo2[,setdiff(colnames(demo2),covariates)] # remove covariates from this df - only use one from processed file
  idemo.performance <- read.csv(paste0(masterdir,'analyses/behavior/idemo/IDEmoAccuracy.csv'),stringsAsFactors = F,row.names = 1)
  if(subj.sample == '22q'){
    # only study 22q subjects - want to explain variation within this phenotype
    demo <- demo[demo$study=='22q',]
    demo.cog <- merge(demo2,demo,by='bblid')
    # ind <- order(match(demo.cog$bblid,demo$bblid)) # not sure why this was in here - don't think it does anything
    # demo.cog <- demo.cog[ind,]
    demo.cog <- cbind(demo.cog,idemo.performance[as.character(demo.cog$scanid),])
  } else if(subj.sample == 'AllSubjects'){
    demo.cog <- merge(demo,idemo.performance,by='scanid') # use this instead of the above 5 lines in order to use all subjects
  } else if(subj.sample == 'PNC'){
    demo <- demo[demo$study=='pnc_sample',]
    demo.cog <- merge(demo,idemo.performance,by='scanid') # use this instead of the above 5 lines in order to use all subjects
  }
  rownames(demo.cog) <- demo.cog$scanid
  # join and merge diagnoses, social cognition, execeff, etc.
  
  grp.colors <- getGroupColors()
  
  # load FIR betas
  savedir <- paste0(masterdir,'analyses/fir/cpc_timecourse_fin',fin,'st',st,'/',component_design,'/')  
  
  stim.types <- list(threat=1,nonthreat=0)
  response.types <- list(correct=1,incorrect=0)
  results <- list()
  load(file = paste0(savedir,'lme_all_trials/',component_design,'_PCTimeCoursesLMEEffectsCorrectThreat.RData'))
  savedir <- paste0(savedir,'predict_cog_22q/',subj.sample,'/')
  dir.create(savedir,recursive = T)
  if(grepl('xcp_36p_despike',name_root)){
    ncomps <- 6
    results <- results[1:ncomps]
  } else if(grepl('xcp_6p_noFilter',name_root)){
    ncomps <- 5
    results <- results[2:6] # remove comp 1 which is global signal
  }
  
  if(subj.sample == '22q'){cog.vars <- c(cog.vars)}
  
  # Add executive function as covariate
  if(subj.sample=='22q'){covariates <- c(covariates)} ####
  
  results.cog <- list() # store all of the R^2 values and p-values in one matrix and post-hoc correct then plot
  # loop through independent variables: PC peak values during each stim-response pair
  stim.type <- 'threat'; response.type <- 'incorrect'
  for(stim.type in names(stim.types)){
    for(response.type in names(response.types)){
      res.name <- paste0(stim.type,response.type) 
      stim.type.idx <- stim.types[[stim.type]]
      response.type.idx <- response.types[[response.type]]
      stim.resp.key <- paste0(stim.type,response.type) # results list is indexed by combination of stimulus and response types
      # load each response type
      results.dnames <- list(paste0('PC',1:length(results)),cog.vars)
      # make results matrix, initialize with NaNs
      res.mat <- matrix(NA,nrow = length(results),ncol=length(cog.vars),dimnames=results.dnames)
      results.cog[[stim.resp.key]]$betas <- results.cog[[stim.resp.key]]$p.values <- res.mat
      results.cog[[stim.resp.key]]$rsq <- results.cog[[stim.resp.key]]$f.p.values <- res.mat
      results.cog[[stim.resp.key]]$df <- res.mat
      
      # extract some kind of summary metric of each PC's time course
      #PC.summary <- lapply(results, function(X) as.data.frame(ranef(X$mdl.best))) # extract random effects
      fun <- function(x) max(abs(x)) #mean(x[x<0]) #max(abs(x)) # specify function for quantifying PC time course
      df.list <- lapply(results, function(X) X$mdl.best$data[complete.cases(X$mdl.best$data) & X$mdl.best$data$Threat == stim.type.idx & X$mdl.best$data$Correct == response.type.idx,])
      PC.summary <- lapply(df.list, function(df) aggregate(df,by=list(scanid=df$scanid),fun)[,c('scanid','Score'),drop=F])
      PC.summary <- lapply(PC.summary, function(X) data.frame(X[setdiff(colnames(X),'scanid')],row.names=X$scanid))
       
      ### for testing ways of quantifying PC time courses  
      # fun <- function(x) x[which.max(abs(x))] #mean(x[x>0]) # specify function for quantifying PC time course
      # fun <- function(x) mean(x[x>0]) # specify function for quantifying PC time course
      # df.list <- lapply(results, function(X) X$mdl.best$data[complete.cases(X$mdl.best$data) & X$mdl.best$data$Threat == stim.type.idx & X$mdl.best$data$Correct == response.type.idx,])
      # PC.summary.wma <- lapply(df.list, function(df) aggregate(df,by=list(scanid=df$scanid),fun)[,c('scanid','Score'),drop=F])
      # PC.summary.wma <- lapply(PC.summary.wma, function(X) data.frame(X[setdiff(colnames(X),'scanid')],row.names=X$scanid))
      # 
      # # min.sid <- 6968;3086; max.sid <- 6565 ;9826 # PC2
      # min.sid <- 9826;8120; max.sid <- 8891
      # X <- df.list[[4]]
      # X[X$scanid==max.sid,]
      # plot(PC.summary.wma[[4]]$Score,PC.summary[[4]]$Score)
      ###
      
      for(y in cog.vars){ # loop through dependent variables (accuracy, etc. whatever is specified)
        # concatenate outcome variable, scores (predictor of interest in outer 2 loops), and specified covariates
        df.PCs <- lapply(1:ncomps, function(PC) cbind(y=demo.cog[,y],PC.summary[[PC]][as.character(demo.cog$scanid),,drop=F],demo.cog[,covariates]))
        # exclude people with any NAs
        df.PCs <- lapply(df.PCs, function(X) X[complete.cases(X),])
        # exclude outlier wrt outcome variable
        #df.PCs <- lapply(df.PCs, function(X) X[!outlier.mask(X$y),])
        # plot distribution of outcome variable to justify using rank
        # the ==0 makes sure this only happens once since outcome variable is always the same
        if(stim.type.idx == 0 & response.type.idx == 0){
          overlap.hist[[y]][[subj.sample]] <- df.PCs[[1]]$y
          p <- ggplot() + geom_histogram(aes(x=df.PCs[[1]]$y)) + theme_bw() +
            xlab(nice.y.name.key.hist[[y]]) + ylab('Count') + theme(text=element_text(size=8))
          ggsave(plot = p, filename = paste0(savedir,'Histogram',y,'.pdf'),
                 width = 3,height = 3, units = "cm",useDingbats=F)
        }
        # replace variables with rank
        
        for(PC in 1:ncomps){df.PCs[[PC]]$y <- rank(df.PCs[[PC]]$y)}
        if(stim.type.idx == 0 & response.type.idx == 0){
          p <- ggplot() + geom_histogram(aes(x=df.PCs[[1]]$y)) + theme_bw() +
            xlab(nice.y.name.key[[y]]) + ylab('Count') + theme(text=element_text(size=8))
          ggsave(plot = p, filename = paste0(savedir,'HistogramRank',y,'.pdf'),
                 width = 3,height = 3, units = "cm",useDingbats=F)
        }
        # fit one model using all variables for each PC
        lm.PCs <- lapply(df.PCs,function(X) lm.beta(lm(y~.,data=X)))
        # first regress out covariates -- shouldn't be necessary to fit this on every PC, but rarely missingness can differ by PCs
        # resid.PCs <- lapply(df.PCs,function(X) residuals(lm(formula=reformulate(termlabels = covariates,response='y',intercept = T),data=X)))
        # # now add residuals to df
        # df.PCs <- lapply(1:length(df.PCs), function(PC) cbind(df.PCs[[PC]],y.r=resid.PCs[[PC]]))
        # # remove covariates and original dependent
        # df.PCs <- lapply(df.PCs, function(X) X[,setdiff(colnames(X),c('y',covariates))])
        # # now fit a model to predict residuals from remaining variables (either score or random effects)
        # lm.PCs <- lapply(df.PCs,function(X) lm.beta(lm(y.r~.,data=X)))
        
        # loop through models for each PC and extract coefficients, p-vals, etc. 
        coef.list <- lapply(lm.PCs,get.coef)
        results.cog[[stim.resp.key]][[y]] <- list()
        results.cog[[stim.resp.key]][[y]]$models <- lm.PCs
        results.cog[[stim.resp.key]][[y]]$coef.tables <- coef.list
        results.cog[[stim.resp.key]]$rsq[,y] <- sapply(lm.PCs,get.rsq)
        results.cog[[stim.resp.key]]$f.p.values[,y] <- sapply(lm.PCs,get.ftest.pval)
        results.cog[[stim.resp.key]]$betas[,y] <- sapply(coef.list,function(X) X['Score','Standardized'])
        results.cog[[stim.resp.key]]$p.values[,y] <- sapply(coef.list,function(X) X['Score','Pr(>|t|)'])
        results.cog[[stim.resp.key]]$df[,y] <- sapply(lm.PCs,function(X) summary(X)$df[2])
        #p.list <- lapply(lm.PCs, function(m) p.xy.flex(x=m$fitted.values,y=m$model$y,xlab = 'Fitted',ylab=y))
        
        # Plot partial residuals of each PC peak val/score alongside y
        p.list <- lapply(lm.PCs, function(m) p.xy.flex(x=get.partial.resids(m,'Score')$x,y=m$model$y,xlab = 'Score',ylab=y,col = '#c85795'))
        p.all <- plot_grid(plotlist = p.list,nrow=1,align='hv')
        ggsave(plot = p.all, filename = paste0(savedir,stim.type,response.type,'Score_Predict',y,'.pdf'),
               width = 18,height = 4, units = "cm",useDingbats=F)
      }
      
    }
  }
  

  # Load all of the results from each stim-response combo for each pc
  # FDR correct over all of those tests
  # plot as a matrix
  names(results.cog) <- nice.x.name.key[names(results.cog)]
  pval.list <- lapply(setNames(names(results.cog),names(results.cog)), function(X) list(p.values=results.cog[[X]]$p.values))
  pval.list <- list.posthoc.correct(pval.list,'fdr')
  
  pval.list.1mat <- lapply(setNames(cog.vars,cog.vars), function(y) sapply(setNames(names(results.cog),names(results.cog)), function(X) results.cog[[X]]$p.values[,y,drop=F]))
  pval.list.1mat <- list.posthoc.correct(pval.list.1mat,'fdr')
  for(y in cog.vars){
    pvals <- pval.list.1mat[[y]]
    # extract matrix of betas
    betas <- sapply(setNames(names(results.cog),names(results.cog)), function(X) results.cog[[X]]$betas[,y,drop=F])
    # extract matrix of degrees of freedom
    df <- sapply(setNames(names(results.cog),names(results.cog)), function(X) results.cog[[X]]$df[,y,drop=F])
    # name matrix rows
    rownames(betas) <- rownames(pvals) <- rownames(df) <- paste0('PC',1:ncomps)
    
    # store those matrices in a list
    stat.list <- list(BetasStd=t(betas),Pfdr=t(pvals),df=t(df))
    # save stats
    lapply(names(stat.list), function(s) write.csv(stat.list[[s]],paste0(savedir,s,'_',y,'.csv')))
    # plot beta matrix with overlayed p values
    p <- imagesc(X=t(betas),cmap='custom1_asymmetric',overlay = t(p.signif.matrix(pvals)),overlay.text.sz = 2,
                 overlay.text.col = 'white',ttl = nice.y.name.key[y],caxis_name = expression(beta)) + nice_cbar() + 
      theme(legend.key.height = unit(0.4,'cm'))
    
    ggsave(plot = p, filename = paste0(savedir,'PeakScoreBetas_',y,'.pdf'),
           width = 9,height = 4.5, units = "cm",useDingbats=F)
    
    # plot strongest positive and strongest negative associations individually
    if(y =='allcorrect'){
      b.pos <- which(betas==max(betas),arr.ind=T)
      PC.pos <- b.pos[,'row']
      Event.pos <- colnames(betas)[b.pos[,'col']]
      m.pos <- results.cog[[Event.pos]][[y]]$models[PC.pos][[1]]
      p <- p.xy.flex(x=get.partial.resids(m.pos,'Score')$x,y=m.pos$model$y,xlab = paste0('PC',PC.pos,' Peak,\n',Event.pos),ylab=nice.y.name.key[y],
                     p=pvals[PC.pos,Event.pos],ptxt = 'p[FDR] == ',rtxt = 'beta ==',r=betas[PC.pos,Event.pos],parse = T,col='#c85795',
                     pt.sz=1,alpha = 0.5,ln.sz=0.75,pos='bottom.right')
      ggsave(plot = p, filename = paste0(savedir,'PeakScoreScatter_PC',PC.pos,Event.pos,'_',y,'.pdf'),
             width = 4.4,height = 4.4, units = "cm",useDingbats=F)
      b.neg <- which(betas==min(betas),arr.ind=T)
      PC.neg <- b.neg[,'row']
      Event.neg <- colnames(betas)[b.neg[,'col']]
      m.neg <- results.cog[[Event.neg]][[y]]$models[PC.neg][[1]]
      p <- p.xy.flex(x=get.partial.resids(m.neg,'Score')$x,y=m.neg$model$y,xlab = paste0('PC',PC.neg,' Peak,\n',Event.neg),ylab=nice.y.name.key[y],
                     p=pvals[PC.neg,Event.neg],ptxt = 'p[FDR] == ',rtxt = 'beta ==',r=betas[PC.neg,Event.neg],parse = T,col='#c85795',
                     pt.sz=1,alpha = 0.5,ln.sz=0.75,pos='top.right')
      ggsave(plot = p, filename = paste0(savedir,'PeakScoreScatter_PC',PC.neg,Event.neg,'_',y,'.pdf'),
             width = 4.4,height = 4.4, units = "cm",useDingbats=F)
  
    }
  }
}  

### analyze betas and compare between PNC and 22q
savedir <- paste0(masterdir,'analyses/fir/cpc_timecourse_fin',fin,'st',st,'/',component_design,'/predict_cog_22q/')  
q22.betas <- read.csv(paste0(savedir,'22q/BetasStd_allcorrect.csv'))[,-1] # remove labels
PNC.betas <- read.csv(paste0(savedir,'PNC/BetasStd_allcorrect.csv'))[,-1] # remove labels

p <- p.xy.flex(x=unlist(PNC.betas),y=unlist(q22.betas),xlab='PNC Betas',ylab='22q Betas',
               r.method='spearman',ttl = 'Relationship between PC Scores and Accuracy in 22q vs. PNC',
              col = 'violet')
ggsave(plot = p, filename = paste0(savedir,'PNCvs22qBetas_allcorrect.pdf'),
       width = 4.4,height = 4.4, units = "cm",useDingbats=F)

### compare distribution of accuracy between PNC and 22q samples
cols <- name(getGroupColors(),c('PNC','22q'))
p <- ggplot() + geom_histogram(aes(x=100*overlap.hist$allcorrect$PNC,fill = 'PNC',alpha=0.2)) +
  geom_histogram(aes(x=100*overlap.hist$allcorrect$`22q`,fill='22q',alpha=0.2)) + scale_fill_manual(values=cols,name='') + scale_alpha(NULL) +
  xlab('Emotion Identification Accuracy (%)') + ylab('Count') + theme_bw() + scale_y_continuous(expand=c(0,0)) + theme(text=element_text(size=8))
ggsave(plot = p, filename = paste0(savedir,'PNCvs22qBetas_allcorrect.pdf'),
       width = 4.4,height = 4.4, units = "cm",useDingbats=F)
