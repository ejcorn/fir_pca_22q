rm(list=ls())
name_root <- 'ScanIDSchaefer200Z1_22q_100reps'
#basedir <- '/data/tesla-data/ecornblath/brain_states_22q/'
basedir <- '~/Dropbox/Cornblath_Bassett_Projects/BrainStates22q/brain_states_22q/'
setwd(basedir)

masterdir <- paste0(basedir,'results/',name_root,'/')
savedir <- paste0(masterdir,'analyses/behavior/idemo/')
dir.create(savedir,recursive = TRUE)

source(paste0(basedir,'code/miscfxns/packages.R'))
source(paste(basedir,'code/plottingfxns/GeomSplitViolin.R',sep=''))
source(paste0(basedir,'code/plottingfxns/plottingfxns.R'))
source(paste0(basedir,'code/miscfxns/miscfxns.R'))

demo <- read.csv(paste(basedir,'data/Demographics',name_root,'.csv',sep=''),stringsAsFactors=F)
RNcolors <- getGroupColors()
demo$bblid <- as.character(demo$bblid)
demo$bblid[nchar(demo$bblid)==5] <- paste0('0',demo$bblid[nchar(demo$bblid) == 5]) # add leading 0 to bblid
stimulus.types <- c('all','anger','fear','happy','neutral','sad')
stimulus.types.scores <- c('all','anger','fear','happy','noe','sad')
#stimulus.types <- c('threat','nonthreat')
rownames(demo) <- as.character(demo$scanid) # only access demo by referring to scanids

# load stimulus timing data
stim.all <- read.csv(paste0(basedir,'data/task/idemo/stimulus/all_3col.txt'),header=F,sep='\t') # read tab delimited file
stim.by.type <- lapply(stimulus.types, function(S) read.csv(paste0(basedir,'data/task/idemo/stimulus/',S,'_3col.txt'),header=F,sep='\t'))

# declare path roots to pnc scores for each subject *in your sample*
pnc.path.roots <- paste0(basedir,'data/task/idemo/pnc_scores/scores/',demo$bblid[demo$study == 'pnc_sample'],'_00',demo$scanid[demo$study == 'pnc_sample'],'_all_')
names(pnc.path.roots) <- demo$scanid[demo$study == 'pnc_sample']

# declare path roots to 22q scores for each subject
bblid_scanid <- paste0(demo$bblid[demo$study == '22q'],'_',demo$scanid[demo$study == '22q']) # start with bblid and scanids in sample
bblid_scanid <- substr(bblid_scanid,start=2,stop = nchar(bblid_scanid)) # remove leading 0 
fnames.scores <- lapply(bblid_scanid, function(id) Sys.glob(paste0(basedir,'data/scores2/',id,'/scores/iDemo2.10/*correct.csv'))) # list all of the files to figure out the prefixes... some subjects have 000000 as their bblid
q22.path.roots <- sapply(fnames.scores, function(f) paste0(unique(substr(f,start=1,stop=132)),'_all_')) # truncate to prefix
#q22.path.roots <- paste0(basedir,'data/scores2/',bblid_scanid,'/scores/iDemo2.10/',demo$bblid[demo$study == '22q'],'_00',demo$scanid[demo$study == '22q'],'_all_')

names(q22.path.roots) <- demo$scanid[demo$study == '22q']

all.path.roots <- c(q22.path.roots,pnc.path.roots) # IMPORTANT: this is done in same order as demo$study: 22q first then PNC
scores.exists <- lapply(stimulus.types, function(S) file.exists(paste0(all.path.roots,S,'_correct.csv')))

loadresponses.by.type <- function(stimulus.types,all.path.roots,response.type){
  # INPUTS:
  # stimulus.types: list of prefixes for facial emotion stimuli corresponding to file names of task response data
  # all.path.roots: Nx1 vector
  correct.by.type <- list()
  for(S in stimulus.types){
    correct.by.type[[S]] <- list()
    for(N in names(all.path.roots)){
      P <- all.path.roots[N]
      if(file.exists(paste0(P,S,'_',response.type,'.csv'))){
        correct.by.type[[S]][[N]] <- read.csv(paste0(P,S,'_',response.type,'.csv'),header=F,sep='\t')	
        if(correct.by.type[[S]][[N]][1] ==0){correct.by.type[[S]][[N]] <- data.frame(V1=numeric(),V2=numeric(),V3=numeric())} # if no responses logged, make df empty
      } else{correct.by.type[[S]][[N]] <- 'missing file'}
    }
  }
  return(correct.by.type)
}
countresponse.rows <- function(response.by.type){
  
  stimulus.types <- names(response.by.type)
  all.names <- names(response.by.type[[stimulus.types[1]]])
  df <- data.frame(scanid = all.names)
  rownames(df) <- all.names
  for(S in stimulus.types){
    df[,S] <- NA
    for(N in all.names){
      if(!is.character(response.by.type[[S]][[N]])){df[N,S] <- nrow(response.by.type[[S]][[N]])} 
    }
  }
  return(df)
}

correct.by.type <- loadresponses.by.type(stimulus.types.scores,all.path.roots,'correct')
incorrect.by.type <- loadresponses.by.type(stimulus.types.scores,all.path.roots,'incorrect')
nr.by.type <- loadresponses.by.type(stimulus.types.scores,all.path.roots,'nr')

df.correct <- countresponse.rows(correct.by.type)
df.incorrect <- countresponse.rows(incorrect.by.type)
df.nr <- countresponse.rows(nr.by.type)

n.items <- sapply(stim.by.type,nrow) # number of items for each stim type
names(n.items) <- stimulus.types.scores # name vector elements for lapply
# convert counts to proportions
for(cname in stimulus.types.scores){
  df.correct[,cname] <- df.correct[,cname]/n.items[cname]
  df.incorrect[,cname] <- df.incorrect[,cname]/n.items[cname]
  df.nr[,cname] <- df.nr[,cname]/n.items[cname]
}
save(df.correct,df.incorrect,df.nr,stimulus.types.scores,file = paste0(savedir,'IDEmoTaskPerformance',name_root,'.RData'))

#########################
### tabulate and plot ###
#########################

# add up all response types: correct, incorrect, no response. should be 1
rowSums(cbind(df.correct$all,df.incorrect$all,df.nr$all),na.rm=TRUE)

# many pnc subjects are missing nr and incorrect data
# on CFN if you run: ls -alh /data/jux/BBL/studies/22q/22qBassett/idemo/pnc/scores
# the files with permissions -rw-rw---- correspond to subjects whose incorrect and nr data didn't copy over

# plot 22q vs. pnc task performance for each emotion for each metric (correct, incorrect, no response)
df.list <- list(Correct=df.correct,Incorrect=df.incorrect,NoResponse=df.nr)

for(df.name in names(df.list)){
  df <- df.list[[df.name]]
  df.plot <- collapse.columns(df,stimulus.types.scores)
  df.plot$metric <- df.name
  df.plot$group <- rep(demo$study,length(stimulus.types.scores))
  p <- p.group.jitter(df = df.plot,yname = 'values',xname = 'names',
                      grpname = 'group',cols = RNcolors,
                      ylabel = df.name,xlabel = 'Stimulus Type',alpha=0.6,ylim=c(0,1))
  ggsave(plot = p, filename = paste0(savedir,'PNCvs22qIDEmo_',df.name,name_root,'.pdf'),height = 5,width =18, units = "cm")
  
}
