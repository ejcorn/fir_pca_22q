dates <- read.table('data/medications/dates.txt',stringsAsFactors = F)
retain.mask <- grep('22q_deleted/1 *',dates$V1) # find only imaging files
# get ID + 1-2 lines down which has folder name with dates in it
dates.ids <- data.frame(bblid=dates[retain.mask,],
                        img.date.1 = dates[retain.mask+1,],
                        img.date.2 = dates[retain.mask+2,],
                        stringsAsFactors = F)
# truncate ids to extract just id
dates.ids$bblid <- substr(dates.ids$bblid,start= nchar(dates.ids$bblid)-5,stop=nchar(dates.ids$bblid)-1)
# truncate dates to get just dates
dates.ids$img.date.1 <- substr(dates.ids$img.date.1,start = 1,stop=8)
dates.ids$img.date.2 <- substr(dates.ids$img.date.2,start = 1,stop=8)
# reformat dates
dates.ids$img.date.1 <- as.Date.character(dates.ids$img.date.1,format = '%Y%m%d')
dates.ids$img.date.2 <- as.Date.character(dates.ids$img.date.2,format = '%Y%m%d')
write.csv(dates.ids,file = 'data/medications/bblid_scandates.csv')
