##################################
# HCP Data
##################################

rm(list=ls())

parcel<-"Power264"
run<-"tfMRI_MOTOR_RL"

setwd(paste0("/Users/yizhao/Dropbox/MyFolder/Biostat-JHU/2019/HCP/Data/S1200/",parcel,"/",run))

dir.data<-paste0("/Volumes/ESD-USB3/HCP/",parcel,"/output")

#################################################
# Power 264 ROI module information

#=======================================
# version 1
ROI.info.v1<-read.csv("/Users/yizhao/Dropbox/MyFolder/Biostat-JHU/2019/HCP/Data/S1200/Power264/mni_module.csv")
ROI.info.v1$ID<-as.numeric(ROI.info.v1$Module)
rownames(ROI.info.v1)<-paste0("ROI",1:nrow(ROI.info.v1))

cc.ROI.v1<-colors()[c(136,564,500,469,50,200,460,654,146,430,464)]
#=======================================

#=======================================
# version 2
ROI.info.v2<-read.csv("/Users/yizhao/Dropbox/MyFolder/Neuroimaging/ROI/Power_264/ROI_module_new/Power264_module.csv")
ROI.info.v2$ID<-as.numeric(ROI.info.v2$module)
rownames(ROI.info.v2)<-paste0("ROI",1:nrow(ROI.info.v2))

cc.ROI.v2<-colors()[read.csv("/Users/yizhao/Dropbox/MyFolder/Neuroimaging/ROI/Power_264/ROI_module_new/Power264_module_color.csv")$R.colors]
#=======================================
#################################################

#################################################
subject<-as.numeric(list.files(dir.data))

# time course
ts<-vector("list",length=length(subject))
names(ts)<-subject
ts.indi<-rep(0,length(subject))
for(i in 1:length(subject))
{
  dir.tmp<-paste0(dir.data,"/",subject[i],"/",run)
  
  if(file.exists(dir.tmp))
  {
    ts.indi[i]<-1
    
    # time series
    env.tmp<-new.env()
    load(paste0(dir.tmp,"/ts.RData"),env.tmp)
    ts[[i]]<-env.tmp$dat
    colnames(env.tmp$dat)<-rownames(ROI.info.v1)
  }
}

# motion parameter
motion<-vector("list",length=length(subject))
names(motion)<-subject
motion.dt=FD<-motion
motion.indi<-rep(0,length(subject))
for(i in 1:length(subject))
{
  dir.tmp<-paste0(dir.data,"/",subject[i],"/",run)
  
  if(file.exists(dir.tmp)&file.exists(paste0(dir.tmp,"/motion.RData")))
  {
    motion.indi[i]<-1
    
    # motion parameter
    env.tmp<-new.env()
    load(paste0(dir.tmp,"/motion.RData"),env.tmp)
    
    motion[[i]]<-env.tmp$motion
    motion.dt[[i]]<-env.tmp$motion_dt
    FD[[i]]<-env.tmp$FD
  }
}

# event onsets
# event type: 
# cue
# lf and rf
# lh and rh
# t
TR<-0.73
onset<-vector("list",length=length(subject))
names(onset)<-subject
ev.files<-c("cue","lf","rf","lh","rh","t")
for(i in 1:length(subject))
{
  dir.tmp<-paste0(dir.data,"/",subject[i],"/",run,"/EVs")
  
  if(file.exists(dir.tmp))
  {
    event<-NULL
    for(j in 1:length(ev.files))
    {
      out.tmp<-read.table(paste0(dir.tmp,"/",ev.files[j],".txt"))
      event<-rbind(event,data.frame(event=rep(ev.files[j],nrow(out.tmp)),out.tmp[,1:2],out.tmp[1:2]/TR))
    }
    colnames(event)<-c("event","start_s","duration_s","start_scan","duration_scan")
    
    onset[[i]]<-event
  }
}
#################################################

save.image("Data.RData")
