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
source(paste0(basedir,'code/miscfxns/miscfxns.R'))
masterdir <- paste(basedir,'results/',name_root,'/',sep='')
savedir <- paste0(masterdir,'analyses/fir/cpc_timecourse_fin',fin,'st',st,'/',component_design,'/pncvs22qcoeff/bootexplained/')

grpnames <- list(PNC='HCs',q22='22q',all='Group')
df.list <- list()
for(coeff.grp in names(grpnames)){
  df.list[[coeff.grp]] <- list()
  for(ts.grp in names(grpnames)[-3]){
    coeff.ts.name <- paste0(grpnames[[ts.grp]],' Data -> ',grpnames[[coeff.grp]],' Coefficients')
    df.boot <- collapse.columns(read.csv(paste0(savedir,coeff.grp,'Coeff_',ts.grp,'TS_Bootstrap.csv')))
    # add Coeff.TS for label and Coeff.Sort in order to sort them by coefficient more easily later on
    df.list[[coeff.grp]][[ts.grp]] <- data.frame(df.boot,Coeff.TS=coeff.ts.name,Coeff.Sort=paste0(coeff.grp,ts.grp),stringsAsFactors = F)
  }
}
df.plot <- do.call(rbind,do.call(rbind,df.list))

ncomps.plot <- 6
PCs.plot <- 2:ncomps.plot
# sort fill by which components are being used
coef.lvls <- unique(df.plot$Coeff.TS)[order(unique(df.plot$Coeff.Sort))]
df.plot$Coeff.TS <- factor(df.plot$Coeff.TS,levels = coef.lvls,ordered = T)
# make bars for median values
df.medians <- sapply(coef.lvls,function(x) 
  sapply(PCs.plot, function(PC) median(df.plot$values[df.plot$Coeff.TS==x & df.plot$names ==paste0('PC',PC)])))
df.medians <- data.frame(df.medians,PC=paste0('PC',PCs.plot),check.names = F)
df.medians <- collapse.columns(df.medians,cnames=coef.lvls,groupby = 'PC')
df.medians$names <- factor(df.medians$names,levels=coef.lvls,ordered=T)
  
p <- ggplot() + geom_col(data=df.medians,aes(x=PC,y=values,fill=names),position = position_dodge(1),alpha=0.5)+
  geom_boxplot(data=df.plot,aes(x=names,y=values,fill=Coeff.TS),position = position_dodge(1),outlier.size = 0.05,size=0.05,outlier.stroke = 0) +
  scale_x_discrete(limits=paste0('PC',PCs.plot),labels=paste0('PC',1:ncomps.plot)) + 
  scale_y_continuous(limits=c(0,NA),expand=c(0,0)) + theme_classic() +
  scale_fill_brewer(palette = 'Paired',name='',limits=coef.lvls,breaks=coef.lvls) + xlab('') + ylab('Variance Explained (%)') +
  theme(text=element_text(size=8),legend.key.size = unit(0.1,'cm'),legend.position = 'bottom')
ggsave(plot = p, filename = paste0(savedir,'VarianceExplainedBootstrappedSpecificity.pdf'),
       width = 9,height = 6, units = "cm",useDingbats=F)
# make scree plot
data.dir <- paste0(masterdir,'analyses/fir/cpc_timecourse_fin',fin,'st',st,'/',component_design,'/')
explained.mat <- readMat(paste0(data.dir,'GroupCPCAComponentsExplained.mat'))$explained
ncomps.scree <- 50 #length(explained.mat)
p <- ggplot() + geom_line(aes(x=1:ncomps.scree,y=explained.mat[1:ncomps.scree])) + geom_vline(xintercept = 6,linetype ='dashed') +
  theme_bw() + xlab('PC #')  + ylab(expression(R^2)) + theme(text=element_text(size=8))
ggsave(plot = p, filename = paste0(savedir,'ScreePlot.pdf'),
       width = 9,height = 9, units = "cm",useDingbats=F)

