library(data.table); library(lubridate); library(tidyverse); library(RcppRoll);library(stringr); library(readr)
library(readxl)
library(mgcv)
theme_set(theme_bw())
h<-theme(plot.title = element_text(hjust = 0.5,size=15))

# workdir <- "C:/Users/Lars Reiter/Onedrive/Arbejde/Narhval-ToGo/Finalized/AccDepth_1hourMean"    # filepath laptop
workdir <- "C:/Users/lenny/OneDrive/Arbejde/Narhval-ToGo/(1) Artikel til Tidsskrift/"    # filepath desktop
setwd(workdir) # set working directory to where your files are


# confidence interval function (t-test statistic)
tt <- function(x,sigma,n) x+qt(c(0.025,0.975),df=n-1)*sigma/sqrt(n)


# load data
df0 <- read.table('Finalhours_GW.txt', sep="\t",header=TRUE)

# Threshold for normal behaviour 
th <- 10 #  40 hours threshold normal behaviour

df <- df0 %>% select(hours,vedba,jerk,Ind) # ignore depth column
ddf <- df0 %>% select(hours,Depth,Ind) # include depth

# activity
dt <- as.data.table(df) # data table format
dt <- na.omit(dt) # remove NA rows
# depth
ddt <- as.data.table(ddf)
ddt <- na.omit(ddt)

# long time means (for each individual)
means <- aggregate(dt[, 2:3], by=list(Ind=dt$Ind,aboveth=dt$hours>th), FUN=mean) # individual means after threshold hours
means <- rename(means[means$aboveth==TRUE,],vedbaPooled=vedba,jerkPooled=jerk) # rename vedba/jerk
meansdepth <- aggregate(ddt[,2],by=list(aboveth=ddt$hours>th,Ind=ddt$Ind),FUN=mean) # individual depth mean after threshold hours
meansdepth <- rename(meansdepth[meansdepth$aboveth==TRUE,],depthPooled = Depth) # individual depth mean after threshold hours
maxdepth <- aggregate(ddt[,2],by=list(aboveth=ddt$hours>th,Ind=ddt$Ind),FUN=max) # individual depth mean after threshold hours
maxdepth <- rename(maxdepth[maxdepth$aboveth==TRUE,],depthmaxPooled = Depth) # individual depth mean after threshold hours
meansdepth <- merge(meansdepth,maxdepth,by="Ind")

# add long time mean column(s) for each dataframe
dt <- merge(dt, means,by="Ind") # for vedba/jerk df
ddt <- merge(ddt,meansdepth,by="Ind") # for depth df

dt <- mutate(dt,deltaV=vedba-vedbaPooled,deltaJ=jerk-jerkPooled)
ddt <- mutate(ddt,delta.depthmean=Depth-depthPooled,delta.depthmax=Depth-depthmaxPooled)

dt$aboveth<-NULL
dt$vedba<-NULL
dt$jerk<-NULL
dt$jerkPooled <-NULL
dt$vedbaPooled <-NULL
ddt$aboveth.x <-NULL
ddt$depthPooled <- NULL
ddt$depthmaxPooled <- NULL
ddt$Depth <- NULL

### FITTING ###
m.V <- gam(deltaV ~ s(hours),data=dt)
m.J <- gam(deltaJ ~ s(hours),data=dt)
m.Dmean <- gam(delta.depthmean ~ s(hours),data=ddt)
m.Dmax <- gam(delta.depthmean ~ s(hours),data=ddt)

dt$predV <- predict(m.V)
dt$predJ <- predict(m.J)
ddt$Dmean <- predict(m.Dmean)
ddt$Dmax <- predict(m.Dmax)


#### DATA PREPARATION: ACTIVITY ####

hrsa <- sort(aggregate(dt$hours,list(dt$Ind),max)$x) # activity table



# Active hours for all whales, and number of occurences of each number of hours
mt <- table(hrsa)
mx <-unlist(lapply(names(mt),strtoi)) # get (unique) hours
my <- as.vector(mt) # get number of occurences of unique hours
# adjust so the last hours (of a single whale) is left out (no CI possible)
my <- head(my,-1) # remove last whale
mx <- head(mx,-1) # remove last whale
my[end(my)[1]] <- my[end(my)[1]]  + 1 # add extra occurence to 2nd last number of hours (for 2nd last whale). +1 because after 2nd last whale we had 1 occurences


# create lists (for storage of dataframes)
bothlist <- list() 
aggbothlist <- list()

# chop male dataframe into dataframe segments (with variable active whales)
bothlist[[1]] <- subset(dt,hours<mx[1]) # first dataframe (for all whales)
for(i in 2:length(my)){
  bothlist[[i]] <- setdiff(subset(dt,hours<mx[i]),subset(dt,hours<mx[i-1])) # complete list of intermediate dataframes
}
# mean over all whales and compute variable confidence bond
degfreedom <- sum(my) # active whales
for(i in 1:length(my)){
  aggbothlist[[i]] <- bothlist[[i]] %>%
    group_by(hours) %>%
    summarise(vedbam = mean(deltaV), 
              vedbastd = sd(deltaV),
              jerkm=mean(deltaJ),
              jerkstd=sd(deltaJ),
              civedbaL=tt(mean(deltaV),sd(deltaV),degfreedom)[1],
              civedbaU=tt(mean(deltaV),sd(deltaV),degfreedom)[2],
              cijerkL=tt(mean(deltaJ),sd(deltaJ),degfreedom)[1],
              cijerkU=tt(mean(deltaJ),sd(deltaJ),degfreedom)[2])
  degfreedom <- degfreedom - my[i]
  print(degfreedom)
  
}
aggboth <- as.data.table(rbindlist(aggbothlist, use.names = TRUE, fill = TRUE)) # bind all data segments together


#### DATA PREPARATION: DEPTH ####

hrsd <- sort(aggregate(ddt$hours,list(ddt$Ind),max)$x) # depth all

# Active hours for males, and number of occurences of each number of hours (for males)
mt <- table(hrsd)
mx <-unlist(lapply(names(mt),strtoi)) # get (unique) hours
my <- as.vector(mt) # get number of occurences of unique hours
# adjust so the last hours (of a single whale) is left out (no CI possible)
my <- head(my,-2) # remove last 2 whales
mx <- head(mx,-2) # remove last 2 whales
my[end(my)[1]] <- my[end(my)[1]]  + 3 # add 3 occurences to 3rd last number of hours (+3 is the number of occurences of whales beyond these hours)

# create lists (for storage of dataframes)
depthlist <- list() 
aggdepthlist <- list()

# chop male dataframe into dataframe segments (with variable active whales)
depthlist[[1]] <- subset(ddt,hours<mx[1]) # first dataframe (for all whales)
for(i in 2:length(my)){
  depthlist[[i]] <- setdiff(subset(ddt,hours<mx[i]),subset(ddt,hours<mx[i-1])) # complete list of intermediate dataframes
}
# TODO: mean over males and females and compute variable confidence bounds
degfreedom <- sum(my) # active whales
for(i in 1:length(my)){
  aggdepthlist[[i]] <- depthlist[[i]] %>%
    group_by(hours) %>%
    summarise(depthm = mean(delta.depthmean),
              depthmax = max(delta.depthmax),
              depthstd = sd(delta.depthmean),
              cidepthL=tt( mean(delta.depthmean),sd(delta.depthmean),degfreedom)[1],
              cidepthU=tt( mean(delta.depthmean),sd(delta.depthmean),degfreedom)[2],
              cidepthmL=tt(max(delta.depthmax),sd(delta.depthmax),degfreedom)[1],
              cidepthmU=tt(max(delta.depthmax),sd(delta.depthmax),degfreedom)[2])
  degfreedom <- degfreedom - my[i]
  
}
aggdepth <- as.data.table(rbindlist(aggdepthlist, use.names = TRUE, fill = TRUE)) # bind all data segments together


#### PLOTTING ####


# VEDBA
ggplot(aggboth,aes(hours,vedbam))+
  geom_line()+
  geom_line(aes(y=civedbaU),linetype="dashed")+
  geom_line(aes(y=civedbaL),linetype="dashed")+
  geom_jitter(alpha=I(0.1),color="gray0")+
  geom_hline(yintercept=0,color="red")+
  xlim(c(0,20))

# JERK
ggplot(aggboth,aes(hours,jerkm))+
  geom_line()+
  geom_line(aes(y=cijerkU),linetype="dashed")+
  geom_line(aes(y=cijerkL),linetype="dashed")+
  geom_jitter(alpha=I(0.1),color="gray0")+
  geom_hline(yintercept=0,color="red")+
  xlim(c(0,20))

# depth MEAN
ggplot(aggdepth,aes(hours,depthm))+
  geom_line()+
  geom_line(aes(y=cidepthU),linetype="dashed")+
  geom_line(aes(y=cidepthL),linetype="dashed")+
  geom_jitter(alpha=I(0.1),color="gray0")+
  geom_hline(yintercept=0,color="red")+
  xlim(c(0,20))

# depth MAX
ggplot(aggdepth,aes(hours,depthmax))+
  geom_line()+
  geom_line(aes(y=cidepthmU),linetype="dashed")+
  geom_line(aes(y=cidepthmL),linetype="dashed")+
  geom_jitter(alpha=I(0.1),color="gray0")+
  geom_hline(yintercept=0,color="red")+
  xlim(c(0,20))

