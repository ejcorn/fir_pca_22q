# this script processes the demographic data that was given to me by Kosha
# produces a master combined demographic file called data/PNC*_sample.csv
# the script ProcessData22q.m excludes match pairs where the 22q subject
# in the pair have no idemo data

rm(list=ls())
args <- commandArgs(TRUE)
basedir <- args[1]
#basedir <- paste0(getwd(),'/')
datadir <- '/cbica/home/mahadeva/22q_project/'

# load data
for(PNC.cohort in c('Clean','PreQC')){ # for samples matched post ('Clean') or pre ('PreQC') quality control exclusions
	source(paste(basedir,'code/miscfxns/miscfxns.R',sep=''))
	demo_22q <- read.csv(paste(basedir,'data/Demo22q.csv',sep=''))
	dx_22q <- read.csv(paste(basedir,'data/n85_22q_deleted_demo_dx_cnb.csv',sep=''))
	demo_PNC <- read.csv(paste(basedir,'data/DemoPNC_',PNC.cohort,'.csv',sep=''))
	
	demo_PNC$study <- 'pnc_sample'
	# reorder 22q dx data and delete extra subjects
	reorder.idx <- match(demo_22q$bblid,dx_22q$bblid)
	dx_22q <- dx_22q[reorder.idx,]
	print(identical(dx_22q$bblid,demo_22q$bblid))
	demo_22q <- cbind(demo_22q,dx_prodromal=dx_22q$dx_prodromal,dx_psychosis = dx_22q$dx_psychosis)
	# mark which PNC subjects are matched with prodromal 22q subjects
	demo_PNC <- cbind(demo_PNC,dx_prodromal=dx_22q$dx_prodromal,dx_psychosis = dx_22q$dx_psychosis)

	# rename some columns to match PNC and 22q samples
	demo_22q <- rename.column(demo_22q,'dti64_tsnr','dti64Tsnr')
	demo_22q <- rename.column(demo_22q,'dti64_meanrelrms','dti64MeanRelRMS')
	demo_22q <- rename.column(demo_22q,'dti64_exclude','dti64Exclude')

	combined.columns <- c(intersect(names(demo_22q),names(demo_PNC)))
	# add these columns that only exist for 22q
	additional.columns <- c('restbold_acquired','dti64_acquired','idemo_acquired','idemo_exclude')
	for(ac in additional.columns){demo_PNC[ac] <- NA}
	demo_PNC$idemo_exclude <- 0

	combined.columns <- c(combined.columns,additional.columns)
	demo_combined <- rbind(demo_22q[,combined.columns],demo_PNC[,combined.columns])

	# for pnc subjects, add scanning session label
	scan.bbl.ses <- read.csv(paste0(basedir,'data/pnc_labels_tinashe.csv')) # this links bblids and scan ids with the fmriprep/BIDS output 'ses-1','ses-2','ses-3'
	demo_combined$pnc_session <- as.numeric(sapply(demo_combined$scanid, function(id) ifelse(id %in% scan.bbl.ses$scanID,yes=substr(scan.bbl.ses$sessionID[scan.bbl.ses$scanID == id],4,4),no=NA))) # get session id if pnc subject, otherwise na
	# add in freesurfer brain volume data --- CAUTION: must be sure that scan sessions link up with this
	#demo_combined$fs_name <- ifelse(demo_combined$study=='pnc_sample',yes='pnc_sample2',no='22q_deleted')
	
	fs.files <- sapply(demo_combined$bblid, function(id) ifelse(id %in% demo_PNC$bblid,  # declare paths to aseg stats for pnc and 22q
		#yes=Sys.glob(paste0('/cbica/home/mahadeva/22q_project/freesurfer53/pnc/',id,'/*',demo_combined$scanid[demo_combined$bblid==id],'/stats/aseg.stats')), # freesurfer from BBL
		yes=Sys.glob(paste0(datadir,'pnc_sample2/derivatives_uniqueSession/freesurfer/sub-',id,'PNC',demo_combined$pnc_session[demo_combined$bblid==id],'/stats/aseg.stats')), # use freesurfer from Arun's fmriprep reconstruction
		no=paste0(datadir,'22q_deleted/derivatives/freesurfer/sub-',id,'/stats/aseg.stats')))
	fs.files <- setNames(fs.files,demo_combined$bblid)
	# for 22q subject 16095, arun just ran this subject as a separate xcp run so it wouldn't average T1s over scanning sessions years apart
	fs.files['16095'] <- paste0(datadir,"22q_deleted/derivatives_sub-16095/freesurfer/sub-16095/stats/aseg.stats")
	#fs.files <- setNames(paste0('/cbica/home/mahadeva/22q_project/',demo_combined$fs_name,'/derivatives/freesurfer/sub-',demo_combined$bblid,'/stats/aseg.stats'),demo_combined$bblid)	 # after rerun of fmriprep use this
	aseg.stats <- sapply(names(fs.files), function(id) ifelse(file.exists(fs.files[[id]]),yes=read.delim(fs.files[[id]],stringsAsFactors=F),no=NA)) 
	aseg.brainsegvol <- sapply(aseg.stats,function(x) x[grep('BrainSegVol, Brain Segmentation Volume,',x)]) # extract line containing brain seg vol
	BrainSegVol <- sapply(aseg.brainsegvol, function(x) ifelse(length(x) > 0,yes=as.numeric(strsplit(x,', ')[[1]][4]),no=NA)) # get brainsegvol if recon all exists otherwise na	
	demo_combined$BrainSegVol <- BrainSegVol
	#demo_combined <- remove.column(demo_combined,'fs_name')
	
	# add handedness data
	hand.PNC <- read.csv(paste0(basedir,'data/n1601_dataRelease_MASTER_demogs.csv'))
	hand.22q <- read.csv(paste0(basedir,'data/22q_handedness.csv'))
	hand.all <- sapply(demo_combined$bblid,function(id) ifelse(id %in% hand.PNC$bblid,yes=hand.PNC$handednessv2[hand.PNC$bblid==id],no=hand.22q$hand[hand.22q$bblid==id]))
	print(demo_combined$bblid)
	print(hand.PNC[,c('bblid','handednessv2')])
	print(hand.all)
	hand.all[hand.all == 9] <- NA
	demo_combined$handedness <- hand.all

	print(colnames(demo_combined))
	write.csv(demo_combined,file=paste0(basedir,'data/PNC_',PNC.cohort,'_sample.csv'),row.names=F,na = 'NA')
	if(PNC.cohort == 'Clean'){rm(list=setdiff(ls(),c('basedir','datadir')))}
}

# pnc.clean <- read.csv(paste0(basedir,'data/PNC_Clean_sample.csv'))
# pnc.preqc <- read.csv(paste0(basedir,'data/PNC_PreQC_sample.csv'))
# pnc.all <- rbind(pnc.clean,pnc.preqc)[,c('bblid','scanid','pnc_session')]
# pnc.all <- pnc.all[!is.na(pnc.all$pnc_session),]
# pnc.all <- pnc.all[!duplicated(pnc.all),]
# write.csv(pnc.all,file=paste0(basedir,'data/PNC_unique_sessions.csv'),row.names=F)