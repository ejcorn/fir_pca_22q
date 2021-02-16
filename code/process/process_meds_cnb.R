rm(list=ls())
args <- commandArgs(TRUE)
basedir <- args[1]
name_root <- args[2]

name_root <- 'CPCA_IDSchaefer200Z1xcp_6p_noFilter'
basedir <- '~/Dropbox/Cornblath_Bassett_Projects/BrainStates22q/fir_pca_22q/'
# #basedir <- '/cbica/home/cornblae/ecornblath/fir_pca_22q/'
setwd(basedir)

masterdir <- paste0(basedir,'results/',name_root,'/')
savedir <- paste0(masterdir,'analyses/sample/')
dir.create(savedir,recursive = TRUE)

source(paste0(basedir,'code/miscfxns/packages.R'))
source(paste(basedir,'code/plottingfxns/GeomSplitViolin.R',sep=''))
source(paste0(basedir,'code/plottingfxns/plottingfxns.R'))
source(paste0(basedir,'code/miscfxns/miscfxns.R'))
source(paste0(basedir,'code/statfxns/statfxns.R'))

demo <- read.csv(paste(basedir,'data/Demographics',name_root,'.csv',sep=''),stringsAsFactors=F)
demo$cnb <- NA
#demo <- read.csv(paste0(basedir,'data/PNC_',PNC.cohort,'_sample.csv'),stringsAsFactors=F)

# get cnb data for 22q subjects
iq <- read.csv(paste0(basedir,'data/iq_data_22q.csv'))
cnb <- read.csv(paste0(basedir,'data/cnb.csv'),stringsAsFactors = F)
# first need to select only CNB data that is within 1 year of scan
# convert dates to date format
cnb$cnb_date <- as.Date.character(cnb$cnb_date,format = '%Y%m%d')
cnb$DOSCAN <-  as.Date.character(cnb$DOSCAN,format = '%d-%b-%y')
cnb$CnbMinusScan <- cnb$cnb_date - cnb$DOSCAN
# now identify CNBs that are far away from date of scan
cnb$RemoteCNBFlag <- is.na(cnb$CnbMinusScan) | abs(cnb$CnbMinusScan)>365

for(id in demo$scanid[demo$study=='22q']){
  if(!cnb$RemoteCNBFlag[cnb$SCANID == id]){
    demo[demo$scanid == id,'cnb'] <- cnb$CIQ_ACCURACY[cnb$SCANID == id]
  }
}

cnb.pnc <- read.csv(paste0(basedir,'data/cnb_pnc_out2.csv'))
# loop through pnc subjects, find cnb closest to scan date, no further than 1 year away from scan
for(id in demo$bblid[demo$study=='pnc_sample']){
  subj.cnb <- cnb.pnc[cnb.pnc$bblid == id,] # get cnbs for 1 subject
  subj.cnb$CnbMinusScan <- subj.cnb$cnbAgeMonths - demo$scanage_months[demo$bblid==id]
  if(min(abs(subj.cnb$CnbMinusScan))<=12){ # if there are any cnbs within 1 year of studied scan
    cnb.idx <- which.min(abs(subj.cnb$CnbMinusScan))# get the cnb closest to scan
    if(length(cnb.idx)>1){print('two CNBs equally close to scan')} # make sure only one CNB
    demo$cnb[demo$bblid == id] <- subj.cnb$cnb_iq_ar[cnb.idx] 
  }
}

# analyze iq data
t.test(demo$cnb[demo$study=='22q'],demo$cnb[demo$study=='pnc_sample'])
iq.cnb <- merge(cnb,iq,by.x = 'bblid',by.y = 'bblid')
cor(iq.cnb[,c('CIQ_ACCURACY','viq','piq_pdi','fsiq_mdi')],use='pairwise.complete.obs')
summary(lm(cnb~study,data=demo)) # R2 = 60% for pnc vs 22q IQ/cnb overall accuracy

# get medication use data for each 22q subject
med <- read.csv(paste0(basedir,'data/medicine_22q.csv'),stringsAsFactors = F)
med$DOSCAN <- as.Date.character(med$DOSCAN,format = '%d-%b-%y')
med$DOENT <- as.Date.character(med$DOENT,format = '%d-%b-%y')
med$DOMED_START <- as.Date.character(med$DOMED_START,format = '%d-%b-%y')
med$DOMED_END <- as.Date.character(med$DOMED_END,format = '%d-%b-%y')
# find antipsychotics - manually reviewed and these keywords get all antipsychotics
d12ag <- c('Antipsychotic','ANTIPSYCHOTIC') 
d12ag.mask <- Reduce('|',lapply(d12ag, function(X) grepl(X,med$PHARMCLASS)))
med$StartMedBeforeScan <- (med$DOMED_START-med$DOSCAN) <0 
# **** if no dates just assume they were on med prior to scan
med$StartMedBeforeScan[is.na(med$StartMedBeforeScan)] <- TRUE
#View(med[d2ag.mask & med$StartMedBeforeScan,])
d12ag.ids <- unique(med$BBLID[d12ag.mask & med$StartMedBeforeScan])
demo$d12ag <- demo$bblid %in% d12ag.ids
# save cnb.iq.med csv
cnb.iq <- demo[,c('scanid','bblid','cnb','d12ag')]
write.csv(x=cnb.iq,file = paste0('data/CNBMeds_',name_root,'.csv'),row.names = F)
