setwd("/Volumes/promise_p/ADNI/RS/derivatives/conn/Denoised_timeseries/conn")

library(RNifti)
library(oro.nifti)
library(readr)

totalmotion <- read_csv("/Users/aiyingzhang/ADNI/wholedata/totalmotion_add.csv",col_types = cols(...1 = col_skip()))

detailmotion <- read_csv("/Users/aiyingzhang/ADNI/wholedata/detailmotion_add.csv", col_names = FALSE)

n <- dim(totalmotion)[1]
sub_filtered <- c() 
motion <- list()

for (h in 1:n){
  num_badmotion <- totalmotion[h,]$num
  time <- totalmotion[h,]$time
  subid <- totalmotion[h,]$id
  if(num_badmotion < time*0.3){
    sub_filtered <- c(sub_filtered,subid)
    mt <- rep(0,time)
    motionfile = t(detailmotion[detailmotion$X1==subid,c(2:time)])
    ind = which(motionfile==1)
    if(length(ind)>0 ){
      for(m in 1:length(ind)){
        if(motionfile[ind[m]+1]==1 && motionfile[ind[m]+2]==1 && ind[m]+2 < time){
          mt[ind[m]+1] = 1
          mt[ind[m]+2] = 1
        }
      }
    }    
    motion[[subid]] <- mt
  }
  
}

subfilt <- data.frame(sub_filtered)
subfilt$id <- substr(sub_filtered,1,8)


subinfo <- read_csv("/Users/aiyingzhang/ADNI/wholedata/subinfo_add.csv",col_types = cols(subid = col_character()))
subinfo <- subinfo[,-1]


id_match <- match(subfilt$id,subinfo$subid)
subselected <- subinfo[id_match,]
subselected$subject <- sub_filtered

write.csv(subselected,"/Users/aiyingzhang/ADNI/wholedata/filteredsubinfo_add.csv")

subselected <- read_csv("/Users/aiyingzhang/ADNI/wholedata/filteredsubinfo_add.csv",col_types = cols(subid = col_character()))

n_ft <- dim(subselected)[1]
for(hf in 1:n_ft){
  subject <- subselected[hf,]$subject 
  motion_ind <- which(motion[[subject]]==1)
  path <- paste("/Volumes/promise_p/ADNI/RS/derivatives/conn/Denoised_timeseries/conn/sub-",subject,"/func",sep='')
  a  <- list.files(path, pattern = 'dsub*')
  file <- paste(path,'/',a,sep = '')
  
  rsfmri= readNifti(file)
  rsfmri[,,,motion_ind] =NaN
  
  TR = 3
  x = dim(rsfmri)[1]
  y = dim(rsfmri)[2]
  z = dim(rsfmri)[3]
  t = dim(rsfmri)[4]
  # Estimate intrinsic neural timescale
  INTMap = array(0, dim = c(x,y,z))
  First_Zero_Crossing  = array(0, dim=c(x,y,z))
  for(i in 1:x){
    for(j in 1:y){
      for(k in 1:z){
        timecourse = rsfmri[i,j,k,5:t] # ignore pre-steady-state frames (~10 seconds / ~15 seconds)
        ac = acf(timecourse, lag.max = length(timecourse)-1, plot = FALSE,na.action = na.pass)$acf
        zc = which(ac<0)[1]
        if(all(timecourse==0) || length(zc)==0 || is.na(zc)){
          INTMap[i,j,k] = NaN
          First_Zero_Crossing[i,j,k] = NaN
        }else{
          INTMap[i,j,k] = sum(ac[2:(zc-1)]) * TR
          First_Zero_Crossing[i,j,k] = zc
        }
      }
    }
  }
  INTMap[gm_mask==0] = 0
  image <- nifti(INTMap)
  fname =  paste('/Users/aiyingzhang/ADNI/wholedata/INTMap/',subject,'_INTMap',sep = '')
  writeNIfTI(image, fname, verbose=TRUE)
  
}