# this script loads behavioral data for 22q and PN
# it accounts for the fact that:
# some subjects have _correct.csv files and some have _correct_corrected.csv files
# (differences is that *corrected.csv is already adjusted by 18s for first 6 volumes)

rm(list=ls())

args <- commandArgs(TRUE)
name_root <- args[1]
numClusters <- as.numeric(args[2])
basedir <- args[3]

name_root <- 'ScanIDSchaefer200Z1_22q_PreQC_XCP1normfcon_100reps'
basedir <- '/data/tesla-data/ecornblath/fir_pca_22q/'
basedir <- '~/Dropbox/Cornblath_Bassett_Projects/BrainStates22q/fir_pca_22q/'

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
pnc.path.roots <- paste0(basedir,'data/task/idemo/pnc_scores/scores_plus_missing040920/',demo[pnc.scanids,'bblid'],'_00',pnc.scanids,'_all_')
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
# separately process correct.csv and corrected.csv files
corrected.scores.exists <- lapply(stimulus.types.scores, function(S) file.exists(paste0(all.path.roots,S,'_correct_corrected.csv')))
uncorrected.scores.exists <- lapply(stimulus.types.scores, function(S) file.exists(paste0(all.path.roots,S,'_correct.csv')))

load.corrected.responses.by.type <- function(stimulus.types,all.path.roots,response.type,ndisc=6,TR=3){
  # INPUTS:
  # stimulus.types: list of prefixes for facial emotion stimuli corresponding to file names of task response data
  # all.path.roots: Nx1 vector of file path roots for each subject's behavioral data. 
  # response.type: specify correct, incorrect, or no response files
  # ndisc: number of volumes discarded at beginning of scan (only matters for using uncorrected files)
  # TR: scan TR (only matters for using uncorrected files)
  #
  # this script will loop through stimulus types, load files for each subject
  # based on path roots in all.path.roots, and attempt to load and correct an uncorrected file if corrected unavailable
  #
  # OUTPUT:
  # list of response data, first level is stimulus type, second level within each stimulus type
  # is a subject
  correct.by.type <- list()
  for(S in stimulus.types){
    correct.by.type[[S]] <- list()
    for(N in names(all.path.roots)){
      P <- all.path.roots[N]
      fname.corrected <- paste0(P,S,'_',response.type,'_corrected.csv')
      fname.uncorrected <- paste0(P,S,'_',response.type,'.csv')
      #print(fname)
      if(file.exists(fname.corrected)){
        correct.by.type[[S]][[N]] <- read.csv(fname.corrected,header=F,sep=' ')
        # print(correct.by.type[[S]][[N]])
        #if(correct.by.type[[S]][[N]][1] ==0){correct.by.type[[S]][[N]] <- data.frame(V1=numeric(),V2=numeric(),V3=numeric())} # if no responses logged, make df empty
      } else if(file.exists(fname.uncorrected) & !file.exists(fname.corrected)){ # if no corrected file, load uncorrected file
        print(paste0('subject ',N,' ',S,' -- using uncorrected file -- subtracting ',TR*ndisc,'s (',ndisc,' vol) from start'))
        df.tmp <- read.csv(fname.uncorrected,header=F,sep='\t') # uncorrected files are tab separated
        df.tmp$V1 <- df.tmp$V1 - TR*ndisc+0.001 # correct by subtracting off timing of discarded warmup volumes
        correct.by.type[[S]][[N]] <- df.tmp
      } else{correct.by.type[[S]][[N]] <- 'missing file'}
    }
  }
  return(correct.by.type)
}

assign.response.timing <- function(X){
  # INPUTS:
  # X: data for a subject from load.corrected.responses.by.type for a given stimulus and response type (correct, incorrect)
  #
  # OUTPUTS:
  # return timing of responses to stimuli
  
  if(is.data.frame(X)){
    response.indices <- round(X$V1[X$V2 >0]) # find time in seconds where there is a stimulus (V2 column says duration of stimulus)
    if(length(response.indices) ==0){response.indices <- NaN}
    return(response.indices)
  } else{return(NaN)} # note that missing files and files that are empty reflecting none of a particular response will both be NaN. this doesn't matter for all, bc everyone got at least 1
} # but it does matter for some of the emotions where it's possible for subjects to miss all of a single emotion

# load timing of each type of response
correct.by.type <- load.corrected.responses.by.type(stimulus.types.scores,all.path.roots,'correct')
nr.by.type <- load.corrected.responses.by.type(stimulus.types.scores,all.path.roots,'nr')
incorrect.by.type <- load.corrected.responses.by.type(stimulus.types.scores,all.path.roots,'incorrect')

for(stim.type in 'noe'){
  stim.timing <- round(stim.by.type[[stim.type]]$V1) # find time in seconds where stimulus type is presented
  # construct indicator vector of incorrect responses for each subject
  timing.correct <- timing.incorrect.nr <- timing.incorrect <- timing.nr <- list() # hold indicator vectors
  for(x.name in names(correct.by.type[[stim.type]])){ # see line 120 - all identical
    i.name <- paste0('s',x.name)
    timing.correct[[i.name]] <- assign.response.timing(correct.by.type[[stim.type]][[x.name]])
    if(x.name == 8419){timing.correct[[i.name]] <- timing.correct[[i.name]] + 6}
    timing.incorrect.nr[[i.name]] <- NaN # only add this if there is data on correct
    if(!is.na(timing.correct[[i.name]][1])){timing.incorrect.nr[[i.name]] <- stim.timing[!stim.timing %in% timing.correct[[i.name]]]} # get indices of where stimulus is presented AND response is NOT correct (NR +incorrect)
    if(length(timing.incorrect.nr[[i.name]]) ==0){timing.incorrect.nr[[i.name]] <- NaN}
    
    timing.incorrect[[i.name]] <- assign.response.timing(incorrect.by.type[[stim.type]][[x.name]])
    if(x.name == 8419){timing.incorrect[[i.name]] <- timing.incorrect[[i.name]] + 6}
    timing.nr[[i.name]] <- assign.response.timing(nr.by.type[[stim.type]][[x.name]])
  }
}

# is NOT correct = incorrect + NR?
i.names <- paste0('s',names(correct.by.type[[stim.type]]))
timing.incorrect.nr.reconstruct <- sapply(i.names, function(i.name) sort(c(timing.incorrect[[i.name]],timing.nr[[i.name]])))
timing.incorrect.nr.reconstruct[sapply(timing.incorrect.nr.reconstruct,length) ==0] <- NaN

diffmask <- sapply(i.names, function(i.name) !identical(timing.incorrect.nr[[i.name]],timing.incorrect.nr.reconstruct[[i.name]]))
timing.incorrect.nr[diffmask]
timing.incorrect.nr.reconstruct[diffmask]

Sys.glob(paste0(all.path.roots['8419'],'*'))
read.csv("/Users/Eli/Dropbox/Cornblath_Bassett_Projects/BrainStates22q/brain_states_22q/data/scores2/15409_5175/scores/iDemo2.10/015409_005175_all_sad_nr_corrected.csv",header=F)

# 8419 has only 59 stimuli, AND the stimulus timing in the corrected files doesn't line up with the expected stimulus timing. appears to be off by two volumes
