rm(list=ls())

args <- commandArgs(TRUE)
name_root <- args[1]
numClusters <- as.numeric(args[2])
basedir <- args[3]

name_root <- 'ScanIDSchaefer200Z1_22q_PreQC_100reps'
basedir <- '/data/tesla-data/ecornblath/brain_states_22q/'
basedir <- '~/Dropbox/Cornblath_Bassett_Projects/BrainStates22q/brain_states_22q/'

setwd(basedir)

masterdir <- paste0(basedir,'results/',name_root,'/')
savedir <- paste0(masterdir,'analyses/behavior/idemo/')
dir.create(savedir,recursive = TRUE)

source(paste0(basedir,'code/miscfxns/packages.R'))
source(paste0(basedir,'code/plottingfxns/plottingfxns.R'))
source(paste0(basedir,'code/miscfxns/miscfxns.R'))

TR <- 3 # repetition time in seconds
nTR <- 204 # number of TRs in scan
demo <- read.csv(paste(basedir,'data/Demographics',name_root,'.csv',sep=''),stringsAsFactors=F)
RNcolors <- getGroupColors()
demo$bblid <- as.character(demo$bblid)
demo$bblid[nchar(demo$bblid)==5] <- paste0('0',demo$bblid[nchar(demo$bblid) == 5]) # add leading 0 to bblid
stimulus.types <- c('all','anger','fear','happy','neutral','sad')
stimulus.types.scores <- c('all','anger','fear','happy','noe','sad')
#stimulus.types <- c('threat','nonthreat')
rownames(demo) <- as.character(demo$scanid) # only access demo by referring to scanids
pnc.scanids <- as.character(demo$scanid[demo$study == 'pnc_sample'])
q22.scanids <- as.character(demo$scanid[demo$study == '22q'])

# load stimulus timing data
stim.all <- read.csv(paste0(basedir,'data/task/idemo/stimulus/all_3col.txt'),header=F,sep='\t') # read tab delimited file
stim.by.type <- lapply(stimulus.types, function(S) read.csv(paste0(basedir,'data/task/idemo/stimulus/',S,'_3col.txt'),header=F,sep='\t'))
names(stim.by.type) <- stimulus.types.scores

# declare path roots to pnc scores for each subject *in your sample*
pnc.path.roots <- paste0(basedir,'data/task/idemo/pnc_scores/scores/',demo[pnc.scanids,'bblid'],'_00',pnc.scanids,'_all_')
pnc.fnames <- lapply(pnc.path.roots,function(p) Sys.glob(paste0(p,'*.csv')))

names(pnc.path.roots) <- pnc.scanids
# declare path roots to 22q scores for each subject
bblid_scanid <- paste0(demo[q22.scanids,'bblid'],'_',q22.scanids) # start with bblid and scanids in sample
bblid_scanid <- substr(bblid_scanid,start=2,stop = nchar(bblid_scanid)) # remove leading 0 
fnames.scores <- lapply(bblid_scanid, function(id) Sys.glob(paste0(basedir,'data/scores2/',id,'/scores/iDemo2.10/*all_correct_corrected.csv'))) # list all of the files to figure out the prefixes... some subjects have 000000 as their bblid
q22.path.roots <- sapply(fnames.scores, function(f) paste0(unique(substr(f,start=1,stop=132)),'_all_')) # truncate to prefix
#q22.path.roots <- paste0(basedir,'data/scores2/',bblid_scanid,'/scores/iDemo2.10/',demo$bblid[demo$study == '22q'],'_00',demo$scanid[demo$study == '22q'],'_all_')

names(q22.path.roots) <- demo$scanid[demo$study == '22q']

all.path.roots <- c(q22.path.roots,pnc.path.roots) # IMPORTANT: this is done in same order as demo$study: 22q first then PNC
scores.exists <- lapply(stimulus.types, function(S) file.exists(paste0(all.path.roots,S,'_incorrect_corrected.csv')))

# my questions:
# for pnc, nobody has the "corrected" files so how do i match those to time series
# for 22q, 2 subjects are missing corrected files ("17631_7746" and "18099_8729")

load.corrected.responses.by.type <- function(stimulus.types,all.path.roots,response.type){
  # INPUTS:
  # stimulus.types: list of prefixes for facial emotion stimuli corresponding to file names of task response data
  # all.path.roots: Nx1 vector of file path roots for each subject's behavioral data. 
  # response.type: specify correct, incorrect, or no response files
  #
  # this script will append stimulus.types to all.path.roots and load data
  # OUTPUT:
  # list of response data, first level is stimulus type, second level within each stimulus type
  # is a subject
  correct.by.type <- list()
  for(S in stimulus.types){
    correct.by.type[[S]] <- list()
    for(N in names(all.path.roots)){
      P <- all.path.roots[N]
      fname <- paste0(P,S,'_',response.type,'_corrected.csv')
      #print(fname)
      if(file.exists(fname)){
        correct.by.type[[S]][[N]] <- read.csv(fname,header=F,sep=' ')
       # print(correct.by.type[[S]][[N]])
        #if(correct.by.type[[S]][[N]][1] ==0){correct.by.type[[S]][[N]] <- data.frame(V1=numeric(),V2=numeric(),V3=numeric())} # if no responses logged, make df empty
      } else{correct.by.type[[S]][[N]] <- 'missing file'}
    }
  }
  return(correct.by.type)
}

correct.by.type <- load.corrected.responses.by.type(stimulus.types.scores,all.path.roots,'correct')
for(stim.type in stimulus.types.scores){
  stim.timing <- round(stim.by.type[[stim.type]]$V1) # find time in seconds where stimulus type is presented
  # construct indicator vector of incorrect responses for each subject
  timing.correct <- timing.incorrect.nr <- list() # hold indicator vectors
  for(x.name in names(correct.by.type[[stim.type]])){
    i.name <- paste0('s',x.name) # matlab doesn't like numeric field names
    X <- correct.by.type[[stim.type]][[x.name]]
    if(is.data.frame(X)){
      correct.timing <- round(X$V1) # find time in seconds where response is correct
      timing.correct[[i.name]] <- correct.timing
      # place ones in positions where stimulus is presented AND response is NOT correct
      timing.incorrect.nr[[i.name]] <- stim.timing[!stim.timing %in% correct.timing] 
    } else{timing.correct[[i.name]] <- NaN; timing.incorrect.nr[[i.name]] <- NaN}
  }
  print('writing mat data')
  writeMat(con = paste0(savedir,'CorrectResponseTiming_',stim.type,'.mat'),indicator=timing.correct,scanids=names(timing.correct))
  writeMat(con = paste0(savedir,'IncorrectNRTiming_',stim.type,'.mat'),indicator=timing.incorrect.nr,scanids=names(timing.incorrect.nr))
  print('done writing mat data')
}

