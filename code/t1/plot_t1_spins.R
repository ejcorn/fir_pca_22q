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
grp.colors <- getGroupColors()

# specify component design and load PCs
component_design <- 'ThreatNonthreatAllStimuliStratified'
CPC <- readMat(paste0(masterdir,'analyses/fir/subject_fir_correct_incorrect_pca/cpc_timecourse/',
                      component_design,'/pncvs22qcoeff/FIRGroup',component_design,'_CPCAComponentsBootstrappedThreshold.mat'))
PCs.oi <- c(1:6) # these are components with altered temporal profiles in 22q
PC.maps <- CPC$nodeData[,PCs.oi]

# load
savedir <- paste0(masterdir,'analyses/t1/pc_spin/',component_design,'/')
dir.create(savedir,recursive = T)
matData <- readMat(paste0(masterdir,'analyses/t1/CTSAPermSchaefer200.mat'))
t1.maps <- list(CT=matData$ct.data,SA=matData$sa.data)
t1.maps$SA <- t1.maps$SA*-1
t1.maps.perm <- list(CT=lapply(1:dim(matData$ct.data.perm)[2],function(y) matData$ct.data.perm[,y]),
                     SA=lapply(1:dim(matData$sa.data.perm)[2],function(y) -1*matData$sa.data.perm[,y]))


# do areas with significant loadings on PCs tend to
# have more structural differences?
# surface area: small values more prominent, indicate increased SA in 22q -- multiply this by -1
# ct: large values more promininet, indicate decreased CT in 22q

# For all of the statistically significant areas of each PC, get the mean beta value
dfun <- function(PCs,t1){return(sapply(1:ncol(PCs), function(PC) mean(abs(t1)[PCs[,PC]!=0],na.rm=T))) }

results <- list()
metric.long.names <- list(CT='Cortical Thickness',SA='Surface Area')
for(t1.metric in names(t1.maps)){
  t1.map.j <- t1.maps[[t1.metric]] # iterate through t1 maps
  t1.maps.perm.j <- t1.maps.perm[[t1.metric]]
  d.obs <- dfun(PC.maps,t1.map.j) # get metric for each PC observed
  d.null <- sapply(t1.maps.perm.j, function(t1.perm) dfun(PC.maps,t1.perm)) # get metric for null PCs
  results[[t1.metric]]$p.spin <- sapply(1:ncol(PC.maps), function(PC) mean(d.obs[PC] < d.null[PC,])) # test if structural disruption is more extreme in direction of 22q-TD diff
  results[[t1.metric]]$d.obs <- setNames(d.obs,paste0('PC',PCs.oi))
  results[[t1.metric]]$d.null <- matrix.to.df(t(d.null),dnames = list(NULL,paste0('PC',PCs.oi)))
  # now get the actual beta values for morph difference within each PC to just do t-tests to see if they are >0
  results[[t1.metric]]$beta.dist <- lapply(1:ncol(PC.maps), function(PC) data.frame(name=paste0('PC',PCs.oi[PC]),values=t1.map.j[PC.maps[,PC]!= 0]))
  results[[t1.metric]]$beta.dist.p <- lapply(results[[t1.metric]]$beta.dist, function(X) t.test(X$values)$p.value)
  results[[t1.metric]]$name <- metric.long.names[[t1.metric]]
}
results$CT$p.spin
results$SA$p.spin
# make plots
pval.list <- list.posthoc.correct(lapply(results, function(X) X$p.spin),'fdr')
for(t1.metric in names(t1.maps)){results[[t1.metric]]$p.spin.fdr <- pval.list[[t1.metric]]}
p.list <- lapply(results, function(X) ggplot() + geom_boxplot(data=collapse.columns(X$d.null),aes(x=names,y=values),alpha=0.6,fill="#0B775E",outlier.stroke = 0,outlier.size = 0.5) +
         geom_point(aes(x=names(X$d.obs),y=X$d.obs),color="#E1BD6D",shape='diamond',size=2) +
           annotate(geom='text',x=names(X$d.obs),y=max(X$d.null),label=ifelse(X$p.spin.fdr<0.05,yes='*',no=''),size=4,color="#F2300F")+
           annotate(geom='text',x=names(X$d.obs),y=max(X$d.null),vjust=1,label=ifelse(X$p.spin<0.05,yes='*',no=''),size=4,color="#35274A")+
           theme_classic()+ggtitle(X$metric)+theme(plot.title = element_text(hjust=0.5))+
           xlab('') + ylab(paste0(X$name,' (MAD)')) + theme(text=element_text(size=8),
          axis.text.x=element_text(angle=90,hjust=1,vjust=0.5)))
p.all <- plot_grid(plotlist=p.list,nrow=1)
ggsave(plot = p.all, filename = paste0(savedir,'PC_CTSA_MAD.pdf'),
       width = 6,height = 5, units = "cm",useDingbats=F)


