rm(list=ls())

args <- commandArgs(TRUE)
name_root <- args[1]
basedir <- args[2]
name_root <- 'ScanIDSchaefer200Z1_22q_PreQC_XCP36pdespike_us_100reps'
basedir <- '~/Dropbox/Cornblath_Bassett_Projects/BrainStates22q/brain_states_22q/'

setwd(basedir)

source(paste0(basedir,'code/miscfxns/packages.R'))
library(rstatix) # can't get this to install on cfn right now

masterdir <- paste(basedir,'results/',name_root,'/',sep='')

source(paste(basedir,'code/plottingfxns/GeomSplitViolin.R',sep=''))
source(paste0(basedir,'code/plottingfxns/plottingfxns.R'))
source(paste0(basedir,'code/miscfxns/miscfxns.R'))
source(paste0(basedir,'code/statfxns/statfxns.R'))

demo <- read.csv(paste(basedir,'data/Demographics',name_root,'.csv',sep=''),stringsAsFactors=F)
rownames(demo) <- as.character(demo$scanid)
grp.colors <- getGroupColors()

# load FIR betas
component_design <- 'ThreatNonthreatAllStimuliStratified'
ncomps <- 6
fin <- 6
TR <- 3
t.names <- as.character(TR*1:fin)
savedir <- paste0(masterdir,'analyses/fir/subject_fir_correct_incorrect_pca/cpc_timecourse/',component_design,'/')
dir.create(paste0(savedir,'anova'),recursive=T)
stim.types <- list(all='all',threat='threat',nonthreat='nonthreat')
response.types <- list(correct='correct',incorrect='incorrect')

matData <- lapply(stim.types, function(S)
  lapply(response.types, function(R) readMat(paste0(savedir,S,R,'FIRBetas_CPCScores.mat'))))

covariates <- c('scanage','sex','BrainSegVol','idemo_meanrelrms','handedness')
stim.type <- stim.types[[1]]
response.type <- response.types[[1]]
PC <-1
results <- list()
for(stim.type in stim.types){
  results[[stim.type]] <- list()
  for(response.type in response.types){
    results[[stim.type]][[response.type]] <- list()
    PC.scores <- matData[[stim.type]][[response.type]]$betaFIR.reshape
    scanids.mat <- as.numeric(unname(unlist(matData[[stim.type]][[response.type]]$subjInd.scanID)))
    p.data.all <- p.pred.all <- list()
    for(PC in 1:ncomps){
      results.PC <- list()
      results.PC$PC <- paste0('PC',PC)
      PC.scores.j <- setNames(as.data.frame(PC.scores[,,PC]),t.names)
      PC.scores.j <- cbind(PC.scores.j,scanid=scanids.mat,Is22q=as.numeric(demo[as.character(scanids.mat),'study']=='22q'),stringsAsFactors=F)
      PC.scores.j <- cbind(PC.scores.j,demo[as.character(scanids.mat),covariates])
      PC.df <- collapse.columns(PC.scores.j,cnames = t.names,groupby = c('scanid','Is22q',covariates))
      colnames(PC.df)[1:2] <- c('Score','Time')
      PC.df$Time <- as.numeric(PC.df$Time)
      
      # do a two-way repeated measures anova
      df.aov <- PC.df[complete.cases(PC.df),]
      df.aov$Time <- as.factor(df.aov$Time)
      df.aov$Is22q <- as.factor(df.aov$Is22q)
      df.aov$scanid <- as.factor(df.aov$scanid)
      df.aov$sex <- as.factor(df.aov$sex)
      results.PC$shapiro <- rstatix::shapiro_test(data = rstatix::group_by(.data=df.aov,Is22q,Time),Score)
      results.PC$res.aov <- rstatix::anova_test(
        data = df.aov, dv = Score, wid = scanid,
        within = Time,between = Is22q,covariate = c(scanage,handedness,sex,BrainSegVol,idemo_meanrelrms)
      )
      results.PC$res.aov.t <- rstatix::get_anova_table(results.PC$res.aov)
      results.PC$df.aov <- df.aov
      
      # plot raw data of PC time courses
      PC.df.mean.se <- PC.df %>%
        group_by(Is22q, Time) %>% get_summary_stats(Score,type='mean_se')
      p.data.all[[PC]] <- ggplot(data = PC.df.mean.se) + geom_line(aes(x = Time,y = mean, color=as.character(Is22q)),alpha=0.5,size=0.5) +
        geom_errorbar(aes(x=Time,ymin=mean-2*se,ymax=mean+2*se, color=as.character(Is22q)),alpha=0.5,size=0.5,width=0.8)+
        ggtitle(paste0('PC',PC,':'))+ xlab('Time (s)') + ylab('Estimated HDR')+
        scale_color_manual(values=grp.colors,limits=c('0','1'),guide='none')+
        theme_bw() + theme(plot.title = element_text(hjust=0.5),text=element_text(size=8))
      results[[stim.type]][[response.type]][[PC]] <- results.PC
    }
    results.sr <- results[[stim.type]][[response.type]]
    # correct group comparisons for multiple comparisons
    p.vals <- p.adjust(sapply(results.sr, function(X) X$res.aov.t[X$res.aov.t$Effect=='Is22q:Time','p']),'fdr')
    p.labels <- paste('p = ',signif(p.vals,2))
    p.data.all <- lapply(1:length(p.data.all), function(PC) p.data.all[[PC]] + annotate(geom='text',x=Inf,y=-Inf,hjust=1,vjust=-0.5,label=p.labels[[PC]],size=1.5))
    
    p.data <- plot_grid(plotlist=p.data.all,nrow=1,rel_widths = rep(1,ncomps),rel_heights = rep(1,ncomps))
    ggsave(plot = p.data, filename = paste0(savedir,'anova/',stim.type,response.type,'_PCTimeCoursesMeanSE.pdf'),
           width = 18,height = 4, units = "cm",useDingbats=F)
  }
}
aov.all <- list()
for(PC in 1:ncomps){
  # concatenate PC scores for all responses and analyze threat-non-threat differences for each PC
  all.df.aov <- do.call(rbind,lapply(c('threat','nonthreat'), function(S)
    do.call(rbind,lapply(c('incorrect'), function(R) cbind(results[[S]][[R]][[PC]]$df.aov,Stim=S,Response=R)))))
  all.df.aov$Stim <- as.factor(all.df.aov$Stim)
  aov.all[[PC]] <- rstatix::anova_test(
    data = all.df.aov, dv = Score, wid = scanid,
    within = c(Time,Stim),between = Is22q#,covariate = c(scanage,handedness,sex,BrainSegVol,idemo_meanrelrms)
  )
}