# this script loads behavioral data for 22q and PN
# it accounts for the fact that:
# some subjects have _correct.csv files and some have _correct_corrected.csv files
# (differences is that *corrected.csv is already adjusted by 18s for first 6 volumes)

rm(list=ls())

args <- commandArgs(TRUE)
name_root <- args[1]
basedir <- args[2]

#name_root <- 'ScanIDSchaefer200Z1_22q_PreQC_XCP36pdespike_us_100reps'
# basedir <- '/data/tesla-data/ecornblath/brain_states_22q/'
#basedir <- '~/Dropbox/Cornblath_Bassett_Projects/BrainStates22q/brain_states_22q/'

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
#q22.path.roots <- sapply(fnames.scores, function(f) paste0(dirname(f),'/',unique(substr(f,start=1,stop=132)),'_all_')) # truncate to prefix
q22.path.roots <- paste0(basedir,'data/scores2/',bblid_scanid,'/scores/iDemo2.10/',demo$bblid[demo$study == '22q'],'_00',demo$scanid[demo$study == '22q'],'_all_')

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
      } else{correct.by.type[[S]][[N]] <- 'missing'}
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
  } else{return('missing')} # if there are no responses but file exists return NaN, if there are no files at all return 'missing'
} 

infer.response.timing <- function(correct.timing,incorrect.timing,nr.timing,stims){
  # INPUTS:
  # correct, incorrect, and nr timing
  # all possible stimuli
  # 
  # OUTPUTS:
  # if only 1 of the 3 are missing, infer the missing one by exclusion and replace
  timings <- list(correct.timing,incorrect.timing,nr.timing)
  missing.mask <- sapply(timings, function(X) X == 'missing')
  if(sum(missing.mask) == 1){
    timings[missing.mask] <- stims[!stims %in% unlist(timings[!missing.mask])]
    print('successfully inferred missing response timing')
    return(timings)
  } else { return(timings)}
}
# load timing of each type of response
correct.by.type <- load.corrected.responses.by.type(stimulus.types.scores,all.path.roots,'correct')
nr.by.type <- load.corrected.responses.by.type(stimulus.types.scores,all.path.roots,'nr')
incorrect.by.type <- load.corrected.responses.by.type(stimulus.types.scores,all.path.roots,'incorrect')

# aggregate performance on task -- start by initializing performance dataframe
idemo.performance <- data.frame(scanid=demo$scanid,bblid=demo$bblid) 
rownames(idemo.performance) <- as.character(idemo.performance$scanid)
for(resp.type in c('correct','incorrect','nr')){
  for(stim.type in c(stimulus.types.scores,'threat','nonthreat')){
    idemo.performance[,paste0(stim.type,resp.type)] <- NA
  }
}
perf.quant <- function(timing){
  # given performance timing, return number of responses accounting for missings, excludeds, and nas
  if(any(timing == 'missing')){return(NaN)
    } else {return(sum(!is.na(timing)))}
}

#identical(names(correct.by.type),names(incorrect.by.type)); identical(names(nr.by.type),names(incorrect.by.type))
nr.threshold <- 0.3 * nrow(stim.all) # ** Exclude subjects who didn't respond to > 30% of items
timing.correct <- timing.incorrect.nr <- timing.incorrect <- timing.nr <- list() # hold indicator vectors
for(stim.type in stimulus.types.scores){
  stim.timing <- round(stim.by.type[[stim.type]]$V1) # find time in seconds where stimulus type is presented
  timing.correct[[stim.type]] <- timing.incorrect.nr[[stim.type]] <- timing.incorrect[[stim.type]] <- timing.nr[[stim.type]] <- list() # hold indicator vectors
  for(x.name in names(correct.by.type[[stim.type]])){ # see line 141 - all identical
    i.name <- paste0('s',x.name)
    timing.correct[[stim.type]][[i.name]] <- assign.response.timing(correct.by.type[[stim.type]][[x.name]])
    # get indices of where stimulus is presented AND response is NOT correct (NR +incorrect)
    #timing.incorrect.nr.rec[[stim.type]][[i.name]] <- NaN # only add this if there is data on correct
    #if(!any(is.na(timing.correct[[stim.type]][[i.name]]))){timing.incorrect.nr.rec[[stim.type]][[i.name]] <-sort(stim.timing[!stim.timing %in% timing.correct[[stim.type]][[i.name]]])} # incorrect + nr is just "not correct"
    
    timing.incorrect[[stim.type]][[i.name]] <- assign.response.timing(incorrect.by.type[[stim.type]][[x.name]])
    timing.nr[[stim.type]][[i.name]] <- assign.response.timing(nr.by.type[[stim.type]][[x.name]])
    idemo.performance[x.name,paste0(stim.type,'correct')] <- perf.quant(timing.correct[[stim.type]][[i.name]]) 
    idemo.performance[x.name,paste0(stim.type,'incorrect')] <- perf.quant(timing.incorrect[[stim.type]][[i.name]])
    idemo.performance[x.name,paste0(stim.type,'nr')] <- perf.quant(timing.nr[[stim.type]][[i.name]])
    # if 1 individual timing file is missing, attempt to infer it from the other two types of responses
    if(sum('missing' %in% c(timing.correct[[stim.type]][[i.name]],timing.incorrect[[stim.type]][[i.name]],timing.nr[[stim.type]][[i.name]])) == 1){ 
      list[timing.correct[[stim.type]][[i.name]],timing.incorrect[[stim.type]][[i.name]],timing.nr[[stim.type]][[i.name]]] <-
        infer.response.timing(timing.correct[[stim.type]][[i.name]],timing.incorrect[[stim.type]][[i.name]],timing.nr[[stim.type]][[i.name]])
    }
    
    timing.incorrect.nr[[stim.type]][[i.name]] <- sort(c(timing.incorrect[[stim.type]][[i.name]],timing.nr[[stim.type]][[i.name]])) # merge incorrect and nr by concatenating them. remove NAs in case they only have inc, nrs, but not both
    if(any(c('missing','excluded') %in% timing.incorrect.nr[[stim.type]][[i.name]])){timing.incorrect.nr[[stim.type]][[i.name]] <- unique(timing.incorrect.nr[[stim.type]][[i.name]])}
    if(length(timing.incorrect.nr[[stim.type]][[i.name]]) ==0){timing.incorrect.nr[[stim.type]][[i.name]] <- NaN} # if after removing NAs there is nothing, then replace with 1 NA
    # ** Exclude subjects who didn't respond to > 30% of all items
    if(length(timing.nr$all[[i.name]]) > nr.threshold){timing.incorrect[[stim.type]][[i.name]] <- timing.correct[[stim.type]][[i.name]] <- timing.nr[[stim.type]][[i.name]] <- 'excluded'}
  }
  # exclude subject 8419, who seems to have a different order of stimulus presentation (see check_extract_TRlocked_responses_v2.R)
  timing.correct[[stim.type]][['s8419']] <- timing.incorrect.nr[[stim.type]][['s8419']] <- timing.incorrect[[stim.type]][['s8419']] <- timing.nr[[stim.type]][['s8419']] <- 'excluded'
  # exclude subject 8484 who has files but no data in them
  timing.correct[[stim.type]][['s8484']] <- timing.incorrect.nr[[stim.type]][['s8484']] <- timing.incorrect[[stim.type]][['s8484']] <- timing.nr[[stim.type]][['s8419']] <- 'excluded'
  
  print(paste0('writing mat data: ',stim.type))
  writeMat(con = paste0(savedir,'CorrectResponseTiming_',stim.type,'.mat'),indicator=timing.correct[[stim.type]],scanids=names(timing.correct[[stim.type]]))
  writeMat(con = paste0(savedir,'IncorrectNRTiming_',stim.type,'.mat'),indicator=timing.incorrect.nr[[stim.type]],scanids=names(timing.incorrect.nr[[stim.type]]))
  writeMat(con = paste0(savedir,'IncorrectTiming_',stim.type,'.mat'),indicator=timing.incorrect[[stim.type]],scanids=names(timing.incorrect[[stim.type]]))
  writeMat(con = paste0(savedir,'NRTiming_',stim.type,'.mat'),indicator=timing.nr[[stim.type]],scanids=names(timing.nr[[stim.type]]))
  print(paste0('done writing mat data: ',stim.type))
}

# generate threat and non threat contrasts by concatenating. threat = anger + fear, non-threat = sad + happy
cat.timing <- function(timing,i.names,stim.types){
  # INPUTS:
  # timing: list whose elements contain stimulus types, the elements of which contain subject timing of idemo task responses
  # i.names: s#### scan ids, names of 2nd layer of timing
  # stim.types: stimulus types to concatenate, correspond to names of 1st layer of timing
  #
  # OUTPUTS:
  # timing.cat: list of subject timings for stimulus types in stim.types
  timing.cat <- sapply(i.names, function(i.name) sort(unlist(unname(sapply(stim.types, function(s) timing[[s]][[i.name]]) ))))
  timing.cat[sapply(timing.cat,length) ==0] <- NaN # if all NaNs, just reduce it to 1 NaN
  timing.cat[sapply(timing.cat,function(X) 'missing' %in% X)] <- 'missing' # if missing just make it 1 missing
  timing.cat[sapply(timing.cat,function(X) 'excluded' %in% X)] <- 'excluded' # if excluded just make it 1 missing
  if(!identical(sapply(timing.cat,length),sapply(timing.cat,function(x) length(unique(x))))){print('warning: duplicates introduced')} # check for duplicates being introduced
  return(timing.cat)
}

i.names <- paste0('s',names(correct.by.type[['all']]))
threat <- c('anger','fear')
timing.correct.threat <- cat.timing(timing.correct,i.names,threat)
timing.incorrect.threat <- cat.timing(timing.incorrect,i.names,threat)
timing.nr.threat <- cat.timing(timing.nr,i.names,threat)
timing.incorrect.nr.threat <- cat.timing(timing.incorrect.nr,i.names,threat)
print('writing mat data threat')
writeMat(con = paste0(savedir,'CorrectResponseTiming_threat.mat'),indicator=timing.correct.threat,scanids=names(timing.correct.threat))
writeMat(con = paste0(savedir,'IncorrectNRTiming_threat.mat'),indicator=timing.incorrect.nr.threat,scanids=names(timing.incorrect.nr.threat))
writeMat(con = paste0(savedir,'IncorrectTiming_threat.mat'),indicator=timing.incorrect.threat,scanids=names(timing.incorrect.threat))
writeMat(con = paste0(savedir,'NRTiming_threat.mat'),indicator=timing.nr.threat,scanids=names(timing.nr.threat))
print('done writing mat data threat')

nonthreat <- c('happy','sad','noe')
timing.correct.nonthreat <- cat.timing(timing.correct,i.names,nonthreat)
timing.incorrect.nonthreat <- cat.timing(timing.incorrect,i.names,nonthreat)
timing.nr.nonthreat <- cat.timing(timing.nr,i.names,nonthreat)
timing.incorrect.nr.nonthreat <- cat.timing(timing.incorrect.nr,i.names,nonthreat)
print('writing mat data nonthreat')
writeMat(con = paste0(savedir,'CorrectResponseTiming_nonthreat.mat'),indicator=timing.correct.nonthreat,scanids=names(timing.correct.nonthreat))
writeMat(con = paste0(savedir,'IncorrectNRTiming_nonthreat.mat'),indicator=timing.incorrect.nr.nonthreat,scanids=names(timing.incorrect.nr.nonthreat))
writeMat(con = paste0(savedir,'IncorrectTiming_nonthreat.mat'),indicator=timing.incorrect.nonthreat,scanids=names(timing.incorrect.nonthreat))
writeMat(con = paste0(savedir,'NRTiming_nonthreat.mat'),indicator=timing.nr.nonthreat,scanids=names(timing.nr.nonthreat))
print('done writing mat data nonthreat')

# add threat correct and incorrect to performance
for(resp.type in c('correct','incorrect','nr')){
  idemo.performance[,paste0('threat',resp.type)] <- rowSums(idemo.performance[,paste0(threat,resp.type)])
  idemo.performance[,paste0('nonthreat',resp.type)] <- rowSums(idemo.performance[,paste0(nonthreat,resp.type)])
}
# normalize all columns to number of stimuli of each type
stim.counts <- sapply(stim.by.type, function(X) nrow(X))
stim.counts['threat'] <- sum(stim.counts[threat])
stim.counts['nonthreat'] <- sum(stim.counts[nonthreat])

for(stim.type in names(stim.counts)){
  idemo.performance[,grepl(stim.type,colnames(idemo.performance))] <- idemo.performance[,grepl(stim.type,colnames(idemo.performance))] / stim.counts[stim.type]
}

write.csv(x=idemo.performance,file = paste0(savedir,'IDEmoAccuracy.csv'))
