rm(list=ls())
args <- commandArgs(TRUE)
basedir <- args[1]
name_root <- args[2]

PNC.cohort <- 'CPCA_IDSchaefer200Z1xcp_36p_despike'
basedir <- '~/Dropbox/Cornblath_Bassett_Projects/BrainStates22q/fir_pca_22q/'
#basedir <- '/cbica/home/cornblae/ecornblath/brain_states_22q/'
setwd(basedir)

masterdir <- paste0(basedir,'results/',name_root,'/')
savedir <- paste0(masterdir,'analyses/sample/')
dir.create(savedir,recursive = TRUE)

source(paste0(basedir,'code/miscfxns/packages.R'))
source(paste(basedir,'code/plottingfxns/GeomSplitViolin.R',sep=''))
source(paste0(basedir,'code/plottingfxns/plottingfxns.R'))
source(paste0(basedir,'code/miscfxns/miscfxns.R'))

demo <- read.csv(paste(basedir,'data/Demographics',name_root,'.csv',sep=''),stringsAsFactors=F)
#demo <- read.csv(paste0(basedir,'data/PNC_',PNC.cohort,'_sample.csv'),stringsAsFactors=F)
RNcolors <- getGroupColors()

p1 <- p.single.jitter(df=demo,yname='scanage',grpname = 'study',
                ylabel = 'Age',xlabel = '',xticklabels = c('22q','PNC'),
                cols=RNcolors,alpha=0.6) + theme(legend.position = 'none')
p2 <- p.single.jitter(df=demo,yname='idemo_meanrelrms',grpname = 'study',
                      ylabel = 'Mean FD',xlabel = '',xticklabels = c('22q','PNC'),
                      cols=RNcolors,alpha=0.6) + theme(legend.position = 'none')
demo$sex <- as.factor(demo$sex)
levels(demo$sex) <- c('male','female')
p3 <- ggplot(demo) + geom_bar(aes(x=study,fill=sex),position = position_dodge()) +
  scale_fill_brewer(palette='Blues')+
  xlab('')+theme_classic() + theme(text = element_text(size = 8),legend.position = 'none')

demo$BrainSegVol <- demo$BrainSegVol/1000 # convert to cm^3
p3 <- p.single.jitter(df=demo,yname='BrainSegVol',grpname = 'study',
                      ylabel = expression(BrainSegVolcm^3),xlabel = '',xticklabels = c('22q','PNC'),
                      cols=RNcolors,alpha=0.6) + theme(legend.position = 'none')

demo$race <- as.factor(demo$race)
p4 <- ggplot(demo) + geom_bar(aes(x=study,fill=race),position = position_dodge()) +
  scale_fill_brewer(palette='Blues')+ ylab('Count')+
  xlab('')+theme_classic() + theme(text = element_text(size = 8),legend.position = 'none')
ggsave(plot=p4+theme(legend.position = 'bottom',legend.key.size = unit(0.1,'cm')),filename = paste0(savedir,'RaceLegend.pdf'))
p.all <- plot_grid(plotlist=list(p1,p2,p3,p4),nrow=2)
ggsave(plot = p.all, filename = paste0(savedir,'SampleCharacteristics',name_root,'.pdf'),
       height = 9,width =9, units = "cm")

sum.22q <- summary(demo$scanage[demo$study=='22q'])
sum.PNC <- summary(demo$scanage[demo$study=='pnc_sample'])
tab <- as.data.frame(rbind(sum.22q,sum.PNC))
rownames(tab) <- c('22q11.2DS','PNC')
tab <- round(tab,digits = 2)
library(xtable)
xtable(tab,caption = 'Table 1. Sample demographics.',label = 'table:table1')

demo$sex <- as.numeric(demo$sex)
data.22q <- demo[demo$study=='22q',]
data.PNC <- demo[demo$study=='pnc_sample',]
q22 <- c(paste(signif(mean(data.22q$scanage),3),'±',signif(sd(data.22q$scanage),2)),
           paste(signif(100*mean(data.22q$sex == 1),3),'%',sep=''),
           paste(signif(100*mean(data.22q$race == 1),3),'%',sep=''),
           paste(signif(100*mean(data.22q$race == 2),3),'%',sep=''),
           paste(signif(100*mean(!data.22q$race %in% c(1,2)),3),'%',sep=''),
           paste(signif(mean(data.22q$BrainSegVol/1000),3),'±',signif(sd(data.22q$BrainSegVol/1000),2))
)
PNC <- c(paste(signif(mean(data.PNC$scanage),3),'±',signif(sd(data.PNC$scanage),2)),
             paste(signif(100*mean(data.PNC$sex == 1),3),'%',sep=''),
             paste(signif(100*mean(data.PNC$race == 1),3),'%',sep=''),
             paste(signif(100*mean(data.PNC$race == 2),3),'%',sep=''),
             paste(signif(100*mean(!data.PNC$race %in% c(1,2)),3),'%',sep=''),
             paste(signif(mean(data.PNC$BrainSegVol/1000),3),'±',signif(sd(data.PNC$BrainSegVol/1000),2))
)
tab <- data.frame(`22q11.2DS`=q22,PNC = PNC,check.names = F)
rownames(tab) <- c('Age (y)','Male','White','African American','Other Race','Total Brain Volume (cm^3)')
xtable(tab,caption = 'Table 1. Sample demographics.',label = 'table:table1')
