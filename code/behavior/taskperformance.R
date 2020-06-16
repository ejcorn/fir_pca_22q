args <- commandArgs(TRUE)
name_root <- args[1]
numClusters <- as.numeric(args[2])
basedir <- args[3]

name_root <- 'ScanIDSchaefer200Z1_22q_100reps'
basedir <- '/data/tesla-data/ecornblath/brain_states_22q/'

masterdir <- paste(basedir,'results/',name_root,'/',sep='')

source(paste0(basedir,'code/miscfxns/packages.R'))
source(paste(basedir,'code/plottingfxns/GeomSplitViolin.R',sep=''))
source(paste0(basedir,'code/plottingfxns/plottingfxns.R'))

demo <- read.csv(paste(basedir,'data/Demographics',name_root,'.csv',sep=''),stringsAsFactors=F)
demo$bblid <- as.character(demo$bblid)
demo$bblid[nchar(demo$bblid)==5] <- paste0('0',demo$bblid[nchar(demo$bblid) == 5]) # add leading 0 to bblid
stimulus.types <- c('anger','fear','happy','neutral','sad')
#stimulus.types <- c('threat','nonthreat')

stim.all <- read.csv(paste0(basedir,'data/task/idemo/stimulus/all_3col.txt'),header=F,sep='\t') # read tab delimited file
stim.by.type <- lapply(stimulus.types, function(S) read.csv(paste0(basedir,'data/task/idemo/stimulus/',S,'_3col.txt'),header=F,sep='\t'))
pnc.path.roots <- paste0(basedir,'data/task/idemo/pnc_scores/scores/',demo$bblid[demo$study == 'pnc_sample'],'_00',demo$scanid[demo$study == 'pnc_sample'],'_all_')
names(pnc.path.roots) <- demo$scanid[demo$study == 'pnc_sample']
q22.path.roots <- paste0(basedir,'data/task/idemo/22q_scores/scores/',demo$bblid[demo$study == '22q'],'_00',demo$scanid[demo$study == '22q'],'_all_')
names(q22.path.roots) <- demo$scanid[demo$study == '22q']
all.path.roots <- c(pnc.path.roots,q22.path.roots)

# check if data exists for pnc subjects in my sample
pnc.exists <- lapply(stimulus.types, function(S) file.exists(paste0(pnc.path.roots,S,'_correct.csv')))
do.call('cbind',pnc.exists)
#task.data.exists <- file.exists(paste0(q22.path.roots,'anger_correct.csv')) # check which files exist... if subject has data for one they have it for all
# most of the bblids in my demo don't have files in the scores file
# so I can't compare scores to brain data in this sample
# run 
# task.data.exists <- file.exists(paste0(q22.path.roots,'anger_correct.csv')) # check which files exist... if subject has data for one they have it for all
# q22.path.roots[!task.data.exists]

# then do ls /data/tesla-data/ecornblath/brain_states_22q/data/task/idemo/22q_scores/scores/*anger_correct.csv | grep 17565 (or any bblid in the above list)

#all.path.roots <- all.path.roots[task.data.exists]

loadresponses.by.type <- function(stimulus.types,all.path.roots,response.type){
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

pnc.correct.by.type <- loadresponses.by.type(stimulus.types,pnc.path.roots,'correct')
pnc.df.correct <- countresponse.rows(pnc.correct.by.type)

correct.by.type <- loadresponses.by.type(stimulus.types,all.path.roots,'correct')
incorrect.by.type <- loadresponses.by.type(stimulus.types,all.path.roots,'incorrect')

df.correct <- countresponse.rows(correct.by.type)
df.incorrect <- countresponse.rows(incorrect.by.type)
