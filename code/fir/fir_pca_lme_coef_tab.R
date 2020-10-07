rm(list=ls())

args <- commandArgs(TRUE)
name_root <- args[1]
basedir <- args[2]
component_design <- args[3]

name_root <- 'CPCA_IDSchaefer200Z1xcp_6p_noFilter'
basedir <- '~/Dropbox/Cornblath_Bassett_Projects/BrainStates22q/fir_pca_22q/'
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

component_design <- 'ThreatNonthreatAllStimuliStratified'
savedir <- paste0(masterdir,'analyses/fir/subject_fir_correct_incorrect_pca/cpc_timecourse/',component_design,'/lme_all_trials/')
load(file = paste0(savedir,component_design,'_PCTimeCoursesLMEEffectsCorrectThreat.RData'))
savedir <- paste0(savedir,'coef_tables/')
dir.create(savedir,recursive = T)
source('code/statfxns/lme_msfxns.R') # functions for mixed effects model selection

## make coefficient tables for each models

ncomps <- length(results)
PC<-1
library(xtable)
latexify <- function(x,sci=TRUE){
  sci.mask <- abs(log(abs(x),base=10))>2
  x[sci.mask] <- sapply(x[sci.mask], function(y) format(y, scientific = sci))
  x[sci.mask] <- sanitize.numbers(x[sci.mask],type = "latex", math.style.exponents = TRUE)
  return(x)
}
signif.gated <- function(x,n=2){
  x[abs(log(abs(x),base=10))>2] <- signif(x[abs(log(abs(x),base=10))>2], n)
  return(x)
}
PC <- 2

if(grepl('xcp_36p_despike',name_root)){
  PC.offset <- 0 
} else if(grepl('xcp_6p_noFilter',name_root)){
  PC.offset <- -1 # offset PC index because of global signal
}

keys <- list('Is22q','scanage_months','sex','BrainSegVol','idemo_meanrelrms','handedness')
names(keys) <- c('22q Status','Age (months)','Sex (Female)','Total Brain Volume','Mean Head Motion','Handedness')
for(PC in 1:ncomps){
  m <- results[[PC]]$mdl.best
  #m.4 <- lme4.ms(results[[PC]]$cfg.best,m$data)
  m.s <- summary(m)
  c.t.fixed <- results[[PC]]$mdl.best.stats
  c.t.fixed[,which(!colnames(c.t.fixed) %in% 'DF')] <- signif(c.t.fixed[,which(!colnames(c.t.fixed) %in% 'DF')],2)
  colnames(c.t.fixed)[1] <- 'Estimate'
  c.t.latex <- as.data.frame(results[[PC]]$mdl.best.stats[,c('Value','DF','p-value')])
  c.t.latex <- as.data.frame(lapply(c.t.latex,function(x) signif(x,2)))
  c.t.latex$LaTeX <- paste0('$\\beta = $ ',latexify(c.t.latex$Value),', $p = $ ',
                                latexify(c.t.latex$`p.value`),', $df =$ ',latexify(c.t.latex$DF,sci=F))
  rownames(c.t.latex) <- rownames(c.t.fixed)
  write.csv(c.t.latex,file = paste0(savedir,'LaTeXPC',PC+PC.offset,'.csv'))
  #write.csv(c.t.fixed,file = paste0(savedir,'PC',PC,'FixedEffects.csv'),row.names = T)
  c.t.random <- intervals(m)$reStruct$scanid
  c.t.random <- c.t.random[grep('sd',rownames(c.t.random)),]
  rownames(c.t.random) <- results[[PC]]$cfg.best$random
  rownames(c.t.random)[1] <- 'Intercept'
  #write.csv(c.t.random,file = paste0(savedir,'PC',PC,'RandomEffects.csv'),row.names = T)
  # rework the tables to join into one
  c.t.fixed <- unname(rbind(c('Coefficient Name',colnames(c.t.fixed)),cbind(rownames(c.t.fixed),c.t.fixed)))
  # give variables nice names
  for(e in names(keys)){c.t.fixed[which(c.t.fixed[,1]==keys[e]),1] <- e} # for covariates usign keys variable
  c.t.fixed[,1] <- gsub('Time1','Time',c.t.fixed[,1])
  c.t.fixed[,1] <- gsub('Time([0-9])','Time^\\1',c.t.fixed[,1])
  c.t.random <- as.matrix(c.t.random)
  c.t.random <- unname(rbind(c('Coefficient Name',colnames(c.t.random)),cbind(rownames(c.t.random),c.t.random)))
  c.t.random[,1] <- gsub('Time1','Time',c.t.random[,1])
  c.t.random[,1] <- gsub('Time([0-9])','Time^\\1',c.t.random[,1])
  m.out <- rbind(c(paste0('PC',PC+PC.offset),'','Fixed Effects','','',''),c.t.fixed,
        rep('',ncol(c.t.fixed)),
        rep('',ncol(c.t.fixed)),
        c('','','Random Effects','','',''),
        cbind(c.t.random,'',''))
  colnames(m.out) <- rep('',ncol(m.out))
  write.csv(m.out,file = paste0(savedir,'PC',PC+PC.offset,'FixedAndRandomEffects.csv'),row.names = F)
}
