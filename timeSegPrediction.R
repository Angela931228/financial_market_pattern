library("zoo")
library("tseries")
library("ggplot2")
library("EMD")
library("xlsx")
library("e1071")
library("hydroGOF")


setwd("/Users/Angela/Desktop/market_state")
hsi_price_raw<- read.xlsx("./data/hk/hsi_data.xlsx",1)
#hsi_price <- hsi_price_raw[c((nrow(hsi_price_raw)-1000):nrow(hsi_price_raw)),]


hsi_price_training <- hsi_price_raw[c(1:7000),]
hsi_price<- hsi_price_training
calAverageFrequency<- function(dt){
  f<- {}
  f<- hsi_price_raw[1,1]
  diff_dt<- diff(dt)
  for(i in 1:(length(diff_dt)-1)){
    if(diff_dt[i]*diff_dt[i+1]<0){
      f<- rbind(f, hsi_price_raw[i+1,1])
    }
  }
  
  fDuration<- diff(f)
  return(mean(fDuration))
}
calSegmentsCount<- function(dt){
  diff_dt<- diff(dt)
  cnt<-0
  for(i in 1:(length(diff_dt)-1)){
    if(diff_dt[i]*diff_dt[i+1]<0){
      cnt<- cnt+1
    }
  }
  return(cnt)
}
createSegments <- function(dt){
  f<- {}
  f<- hsi_price_raw[1,1]
  pr<- dt[1]
  diff_dt<- diff(dt)
  for(i in 1:(length(diff_dt)-1)){
    if(diff_dt[i]*diff_dt[i+1]<0){
      f<- rbind(f, hsi_price_raw[i+1,1])
      pr<- rbind(pr, dt[i+1])
    }
  }
  
  
  Duration<- diff(f)
  changeRate<- round(diff(log(pr)),digits = 3)
  result <- cbind(Duration, changeRate)
  return(result)
}
createStatesFromSegments <- function(segments,nSeg){
  i <- 1
  statesList <- {}
  while(i < (nrow(segments)-nSeg)){
    state<-{}
    for(k in 1:nSeg){
      state<- cbind(state,segments[i+k,1])
      state<- cbind(state,segments[i+k,2])
    }
    # i <- i + nSeg # sliding nSeg
    i<- i + 1
    statesList<- rbind(statesList,state)
  }
  return(statesList)
} 

hsi_emd <- emd(as.numeric(as.vector(hsi_price[,2])),hsi_price[,1],boundary="wave")
hsi_emd_res <- as.numeric(as.vector(hsi_price[,2]))- hsi_emd$imf[,1] 

calAverageFrequency(hsi_emd_res) 
calSegmentsCount(hsi_emd_res) #1434

hsi_segments <- createSegments((hsi_emd_res))

set.seed(20)
nSeg<-5
nCluster<- 10

hsi_states <-  createStatesFromSegments(hsi_segments,5)
hsi_state_mean<- apply(hsi_states,2,mean)
hsi_state_sd<- apply(hsi_states,2,sd)
hsi_states_s<-{}
for(i in 1:length(hsi_state_mean)){
  hsi_states_s<- cbind(hsi_states_s,(hsi_states[,i] -hsi_state_mean[i])/hsi_state_sd[i])
}
statesCluster <- kmeans(hsi_states_s, nCluster, nstart = 100)

getStatePoints<- function(state){
  statePoints<- {}
  statePoint<- c(0,0)
  i<-1
  statePoints<- rbind(statePoints,statePoint)
  while(i  < ((length(state)))){
    statePoint <-  cbind( statePoint[1] + state[i], (statePoint[2]+1)*(1+state[i+1])-1)
    i<- i+2
    statePoints<- rbind(statePoints,statePoint)
  }
  return(statePoints)
}

originCluster<-{}

for(i in 1:nCluster){
  originCluster<-rbind(originCluster,(statesCluster$centers[i,])*hsi_state_sd+ hsi_state_mean)
}
par(mfrow = c(5,4), mar=c(4.5,3,2,3))

#plot pattern
pdf(paste("statePattern", nCluster, ".pdf",sep=''))
for(i in 1:nCluster){
  statePoint1 <- getStatePoints(originCluster[i,])
  plot(statePoint1,type = 'l',main= i)
}
dev.off()

#plot patterns belong to cluster
pdf(paste("statesPlot", nCluster,'clusters', ".pdf",sep=''))
for(l in 1:nCluster){
  for(i in 1:200){
    if(statesCluster$cluster[i] == l){
      statePoint1 <- getStatePoints(hsi_states[i,])
      plot(statePoint1,type = 'l',main= statesCluster$cluster[i])
    }
  }
}
dev.off()



h_norm<- function(dt){
  result_norm<- (dt -hsi_state_mean[c(1:(nSeg*2-2))])/hsi_state_sd[c(1:(nSeg*2-2))]
  return(result_norm)
}

checkClusterIndex <- function(dt, dc){
  minRms <- 999999
  index<- 1
  seg_state<-{}
  for(k in 1:nrow(dt)){
    seg_state<- cbind(seg_state,dt[k,1])
    seg_state<- cbind(seg_state,dt[k,2])
  }
  for(i in 1:nCluster){
    rms <- rmse(as.vector(h_norm(seg_state)),as.vector(dc[i,c(1:(nSeg*2-2))]))
    if(rms<minRms){
      minRms <- rms
      index<-i
    }
  }
  return(index)
}
checkRT<- function(a,b,c){
  return( (abs(log(a)-log(b))>abs(c)) &&((a-b)*c>0) )
}
#Calculate return
m<- 7001
rts<-{}
while(m < nrow(hsi_price_raw)){
  hsi_price_t<- hsi_price_raw[c((m-7000):m),]
  hsi_emd_t <- emd(as.numeric(as.vector(hsi_price_t[,2])),hsi_price_t[,1],boundary="wave")
  hsi_emd_res_t <- as.numeric(as.vector(hsi_price_t[,2]))- hsi_emd_t$imf[,1] 
  #hsi_segments <- createSegments((hsi_emd_res))
  hsi_emd_res_der_t <- diff(hsi_emd_res_t)
  if(hsi_emd_res_der_t[length(hsi_emd_res_der_t)]*hsi_emd_res_der_t[length(hsi_emd_res_der_t)-1]<0){
    hsi_segments <- createSegments((hsi_emd_res_t))
    idx<- checkClusterIndex(hsi_segments[c((nrow(hsi_segments)-nSeg+2):nrow(hsi_segments)),],statesCluster$centers)
    dr_<- originCluster[idx,nSeg*2-1]
    cr_<- originCluster[idx,nSeg*2]
    h<-m+1
    hsi_price_h<- hsi_price_raw[c((h-7000):h),]
    hsi_emd_h <- emd(as.numeric(as.vector(hsi_price_h[,2])),hsi_price_h[,1],boundary="wave")
    hsi_emd_res_h <- as.numeric(as.vector(hsi_price_h[,2]))- hsi_emd_h$imf[,1] 
    hsi_emd_res_der_h <- diff(hsi_emd_res_h)
    while(!checkRT(hsi_price_raw[h,2],hsi_price_raw[m,2],cr_)&&( (difftime(hsi_price_raw[h,1],hsi_price_raw[m,1]))<dr_)&&(hsi_emd_res_der_h[length(hsi_emd_res_der_h)]*hsi_emd_res_der_h[length(hsi_emd_res_der_h)-1]>0)){
       h<-h+1
    }
    if(h<nrow(hsi_price_raw)){
      if((hsi_price_raw[h,2]-hsi_price_raw[m,2])*cr_>0){
        rt<-abs(log(hsi_price_raw[h,2])-log(hsi_price_raw[m,2]))
      }else{
        rt<- -abs(log(hsi_price_raw[h,2])-log(hsi_price_raw[m,2]))
      }
      print(rt)
      rts<- rbind(rts,rt)
    }
    m<- h
  }else{
    m<-m+1
  }
}


