rm(list=ls())

args <- commandArgs(TRUE)
name_root <- args[1]
basedir <- args[2]
component_design <- args[3]
fin <- 6
st <- 1

# basedir <- '~/Dropbox/Cornblath_Bassett_Projects/BrainStates22q/fir_pca_22q/'
# name_root <- 'CPCA_IDSchaefer200Z1xcp_6p_noFilter'
# #basedir <- '/cbica/home/cornblae/ecornblath/fir_pca_22q/'
# component_design <- 'ThreatNonthreatAllStimuliStratified'

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
static.vars <- readMat(paste0(basedir,'data/Demographics',name_root,'.mat'))
Hippocampus <- 209:210
Amygdala <- 211:212 # index of amygdala and hippocampus
Accumbens <- 213:214
Amygdala.Hippocampus <- c(Amygdala,Hippocampus)

for(sc_data_name in c('QA_Pass','GFA_Pass','SIFT_radius2')){
  matFile.SC <- readMat(paste0(basedir,'data/sc/StructuralConnectivity',sc_data_name,name_root,'.mat'))
  
  # difference in weighted degree
  
  SC <- matFile.SC$SC
  dti.missing <- matFile.SC$dti.missing
  ids <- matFile.SC$SC.IDs
  strength <- do.call(rbind,setNames(lapply(1:dim(SC)[3],function(N) rowSums(SC[,,N])),ids))
  colnames(strength) <-paste0('Region',1:ncol(strength))
  
  if(!identical(rownames(strength),as.character(demo$bblid))){print('ERROR: IDs do not match')}
  
  # compare weighted degree between groups while controlling for relevant variables
  demo$is22q <- demo$study == '22q' # make binary variable for 22q
  m.list <- lapply(colnames(strength), function(S) lm.beta(lm(s~scanage_months+sex+handedness+BrainSegVol+is22q,data=cbind(demo,s=strength[,S]))))
  get.lm.pval <- function(m,coef.name='is22qTRUE'){return(summary(m)$coef[coef.name,'Pr(>|t|)'])}
  c.list <- do.call(rbind,lapply(m.list,coefficients))
  p.list <- p.adjust(sapply(m.list,get.lm.pval),method='fdr')
  
  # plot significant differences on brain
  # fname <- paste0(savedir,'22qMinusPNC_RegionalStrength',sc_data_name,'.mat')
  # nodeData <- unname(c.list[1:200,'is22qTRUE',drop=FALSE])#*(p.list[1:200] < 0.05))
  # writeMat(fname,nodeData=nodeData,plotTitles=c('Strength, 22q-PNC',' '))
  # system(paste0('source ~/.bash_profile; source activate pyforge ; python code/visualize/brainvis_fir.py ',name_root,' ',static.vars$atlasScale,' ',fname))

  # now weight strength by each component of FIR
  
  savedir <- paste0(masterdir,'analyses/sc_vs_pc/',component_design,'/')
  dir.create(savedir,recursive = T)
  CPC <- readMat(paste0(masterdir,'analyses/fir/cpc_timecourse_fin',fin,'st',st,'/',
                        component_design,'/pncvs22qcoeff/FIRGroup',component_design,'_CPCAComponentsBootstrappedThreshold.mat'))
  CPC.coeff <- CPC$nodeDataAll
  if(grepl('xcp_36p_despike',name_root)){
    n.pcs <- 6
    PCs.label <- PCs.idx <- 1:n.pcs
  } else if(grepl('xcp_6p_noFilter',name_root)){
    n.pcs <- 5
    PCs.idx <- 2:6 # remove comp 0 which is global signal
    PCs.label <- 1:n.pcs # for plots, call PC1 PC0 b/c it's global signal
    names(PCs.label) <- PCs.idx
  }
  #CPC.weighted.strength <- strength %*% (CPC.coeff^2) # weight strength by CPC maps
  #CPC.weighted.strength <- sapply(1:n.pcs, function(PC) rowMeans(strength[,CPC.coeff[,PC]!=0]))
  CPC.weighted.strength <- matrix(NA,nrow=length(ids),ncol=n.pcs,dimnames=list(ids,paste0('CPCDegreeWeight',PCs.idx)))
  for(PC.i in PCs.idx){
    for(id in ids){
      PC.mask <- CPC.coeff[,PC.i]!=0 # here just average connectivity within areas that are up or down
      PC.mask1 <- CPC.coeff[,PC.i]>0; PC.mask2 <- CPC.coeff[,PC.i]<0; # here average connections between up areas and down areas
      CPC.weighted.strength[as.character(id),paste0('CPCDegreeWeight',PC.i)] <- mean(SC[PC.mask1,PC.mask2,which(ids==id)])
      CPC.weighted.strength[as.character(id),paste0('CPCDegreeWeight',PC.i)] <- mean(SC[Amygdala,PC.mask,which(ids==id)])
    }
  }
  colnames(CPC.weighted.strength) <- paste0('CPCDegreeWeight',PCs.label)
  # compute global strength to use as a covariate
  global.strength <- do.call(rbind,setNames(lapply(1:dim(SC)[3],function(N) mean(SC[Amygdala,,N])),ids))
  
  CPC.weighted.strength <- cbind(CPC.weighted.strength, GlobalStrength=global.strength[,1])
  
  # regress weighted strength/connectivity on  covariates + global strength
  m.list <- lapply(paste0('CPCDegreeWeight',1:n.pcs), function(P) lm(s~scanage_months+sex+handedness+BrainSegVol+is22q +dti64MeanRelRMS,
                                                                   data=cbind(demo,s=CPC.weighted.strength[,P],GlobalStrength=CPC.weighted.strength[,'GlobalStrength'])))
  # add analysis of global strength alone
  m.list <- c(m.list,list(lm(GlobalStrength~scanage_months+sex+handedness+BrainSegVol+is22q+dti64MeanRelRMS,data=cbind(demo,GlobalStrength=CPC.weighted.strength[,'GlobalStrength']))))
  # extract betas for 22q status
  get.lm.coef <- function(m,coef.name='is22qTRUE'){return(summary(m)$coef[coef.name,'Estimate'])}
  c.list <- do.call(rbind,lapply(m.list,get.lm.coef))
  p.list <- p.adjust(sapply(m.list,get.lm.pval),method='fdr')
  
  # compute partial residuals for each regression: coefficient*X + residuals
  df.plt <- do.call(rbind,lapply(1:length(m.list),function(j) data.frame(j=paste('PC',j),y=residuals(m.list[[j]]) + c.list[j]*m.list[[j]]$model$is22q,x=m.list[[j]]$model$is22q,stringsAsFactors=F)))
  df.plt$j[df.plt$j == paste('PC',n.pcs+1)] <- 'Global Strength'
  # plot
  
  x.tick.labs <- c(paste('PC',1:n.pcs),'Global Strength')
  #p.start <- ggplot(df.plt) + geom_jitter(aes(x=j,y=y,color=x),position=position_jitterdodge())
  p.start <- ggplot(df.plt) + geom_boxplot(aes(x=j,y=y,fill=x),position=position_dodge(),size=0.25,outlier.size = 0.25)
  p <- p.start + scale_fill_manual(breaks=c(TRUE,FALSE),labels=c('22q','PNC'),values=getGroupColors()) +
    annotate(geom='text',x=x.tick.labs,y=1.1*max(df.plt$y),label=p.signif.matrix(as.matrix(p.list)),size=2.5,color='red')+
    scale_x_discrete(limits=x.tick.labs)+ theme_classic() + theme(axis.text.x=element_text(size=8,angle=90,hjust=1,vjust=0.5))+
    theme(text=element_text(size=8),legend.title = element_blank(),legend.key.size = unit(0.1,'cm')) + xlab('') +ylab('PC-Weighted Strength')
  p
  ggsave(p,filename = paste0(savedir,'CPCWeightedStrength_GlobalStrength',sc_data_name,'.pdf'),
         units = 'cm',height = 6,width = 6)
  p
}  
  
