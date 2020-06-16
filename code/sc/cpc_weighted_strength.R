rm(list=ls())

args <- commandArgs(TRUE)
name_root <- args[1]
numClusters <- as.numeric(args[2])
basedir <- args[3]

name_root <- 'ScanIDSchaefer200Z1_22q_PreQC_XCP36pdespike_us_100reps'
basedir <- '/data/tesla-data/ecornblath/brain_states_22q/'
basedir <- '~/Dropbox/Cornblath_Bassett_Projects/BrainStates22q/brain_states_22q/'

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

sc_data_name <- 'StreamlineCountPass'
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
m.list <- lapply(colnames(strength), function(S) lm(s~scanage_months+sex+handedness+BrainSegVol+is22q,data=cbind(demo,s=strength[,S])))
get.lm.pval <- function(m,coef.name='is22qTRUE'){return(summary(m)$coef[coef.name,'Pr(>|t|)'])}
c.list <- do.call(rbind,lapply(m.list,coefficients))
p.list <- p.adjust(sapply(m.list,get.lm.pval),method='fdr')

# plot significant differences on brain
# fname <- paste0(savedir,'22qMinusPNC_RegionalStrength',sc_data_name,'.mat')
# nodeData <- unname(c.list[,'is22qTRUE',drop=FALSE]*(p.list < 0.05))
# writeMat(fname,nodeData=nodeData,plotTitles=c('Strength, 22q-PNC',' '))
# system(paste0('source ~/.bash_profile; source activate pyforge ; python code/assesscluster/brainvis_fir.py ',name_root,' ',static.vars$atlasScale,' ',fname))

# now weight strength by each component of FIR

component_design <- 'ThreatNonthreatAllStimuliStratified'
savedir <- paste0(masterdir,'analyses/sc_vs_pc/',component_design,'/')
dir.create(savedir,recursive = T)
CPC <- readMat(paste0(masterdir,'analyses/fir/subject_fir_correct_incorrect_pca/cpc_timecourse/',
                      component_design,'/pncvs22qcoeff/FIRGroup',component_design,'_CPCAComponentsBootstrappedThreshold.mat'))
CPC.coeff <- CPC$nodeData
n.pcs <- ncol(CPC.coeff)
#CPC.weighted.strength <- strength %*% (CPC.coeff^2) # weight strength by CPC maps
#CPC.weighted.strength <- sapply(1:n.pcs, function(PC) rowMeans(strength[,CPC.coeff[,PC]!=0]))
CPC.weighted.strength <- matrix(NA,nrow=length(ids),ncol=n.pcs,dimnames=list(ids,NULL))
for(PC.i in 1:n.pcs){
  for(id in ids){
    PC.mask <- CPC.coeff[,PC.i]!=0 # here just average connectivity within areas that are up or down
    PC.mask1 <- CPC.coeff[,PC.i]>0; PC.mask2 <- CPC.coeff[,PC.i]<0; # here average connections between up areas and down areas
    CPC.weighted.strength[as.character(id),PC.i] <- mean(SC[PC.mask1,PC.mask2,which(ids==id)])
  }
}
colnames(CPC.weighted.strength) <- paste0('CPCDegreeWeight',1:n.pcs)
# compute global strength to use as a covariate
global.strength <- do.call(rbind,setNames(lapply(1:dim(SC)[3],function(N) mean(SC[,,N])),ids))

CPC.weighted.strength <- cbind(CPC.weighted.strength, GlobalStrength=global.strength[,1])

# regress weighted strength/connectivity on  covariates + global strength
m.list <- lapply(paste0('CPCDegreeWeight',1:n.pcs), function(P) lm(s~scanage_months+sex+handedness+BrainSegVol+is22q +GlobalStrength+dti64MeanRelRMS,
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


