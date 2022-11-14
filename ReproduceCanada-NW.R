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
df0 <- read.table('Finalhours-New.txt', sep="\t",header=TRUE)

# Threshold for normal behaviour 
th <- 40 #  40 hours threshold normal behaviour

HTth <- 60 # intersect whales at 60 minutes
# handling times
HTdf <- read_excel("HT.xlsx")
whales <- unique(df0$Ind)
maleNames <- c('Asgeir','Balder','Helge','Thor','Kyrri','Siggi','Nemo','Frederik','Eske','Bjarne')
HTdf <- HTdf[HTdf$Ind %in% whales,]
colnames(HTdf)<-c("Ind","handlingtime","size")
HTdf <- mutate(HTdf,HT=ifelse(HTdf$handlingtime < HTth, "Short handlingtime","Long handlingtime"),
               sizeId = ifelse(HTdf$size < quantile(HTdf$size,0.25),1,ifelse(HTdf$size<quantile(HTdf$size,0.75),2,3)),
               Sex = ifelse(HTdf$Ind %in% maleNames,"Male","Female"))

HTdf$handlingtime<-NULL
HTdf$size <- NULL

df <- df0 %>% select(hours,vedba,jerk,Ind) # ignore depth column
ddf <- df0 %>% select(hours,Depth,Ind) # include depth

# activity
dt <- as.data.table(df) # data table format
dt <- na.omit(dt) # remove NA rows
dt<-merge(dt,HTdf,by="Ind")
# depth
ddt <- as.data.table(ddf)
ddt <- na.omit(ddt)
ddt<-merge(ddt,HTdf,by="Ind")

# long time means (for each individual)
means <- aggregate(dt[, 3:4], by=list(Ind=dt$Ind,aboveth=dt$hours>th), FUN=mean) # individual means after threshold hours
means <- rename(means[means$aboveth==TRUE,],vedbaPooled=vedba,jerkPooled=jerk) # rename vedba/jerk
meansdepth <- aggregate(ddt[,3],by=list(aboveth=ddt$hours>th,Ind=ddt$Ind),FUN=mean) # individual depth mean after threshold hours
meansdepth <- rename(meansdepth[meansdepth$aboveth==TRUE,],depthPooled = Depth) # individual depth mean after threshold hours
maxdepth <- aggregate(ddt[,3],by=list(aboveth=ddt$hours>th,Ind=ddt$Ind),FUN=max) # individual depth mean after threshold hours
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

#### FITTING ####

#m.V <- gam(deltaV ~ s(hours,by=factor(Sex))+HT,data=dt)
#m.J <- gam(deltaJ ~ s(hours,by=factor(Sex))+HT,data=dt)
#m.Dmean <- gam(delta.depthmean ~ s(hours,by=factor(HT))+Sex,data=ddt)
#m.Dmax <- gam(delta.depthmean ~ s(hours,by=factor(HT))+Sex,data=ddt)

m.V <- gam(deltaV ~ s(hours,by=factor(HT)),data=dt)
m.J <- gam(deltaJ ~ s(hours,by=factor(HT)),data=dt)
m.Dmean <- gam(delta.depthmean ~ s(hours,by=factor(HT)),data=ddt)
m.Dmax <- gam(delta.depthmean ~ s(hours,by=factor(HT)),data=ddt)

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
    group_by(hours,HT) %>%
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
aggboth <- head(aggboth,-9)


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
    group_by(hours,HT) %>%
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
  geom_point(data=dt,aes(y=deltaV),alpha=I(0.1),color="gray0")+
  geom_hline(yintercept=0,color="red")+
  xlim(c(0,100))+
  facet_wrap(.~HT)
ggsave("MRvedbaNW.pdf",width=10,height=5.5,dpi=600)

# JERK
ggplot(aggboth,aes(hours,jerkm))+
  geom_line()+
  geom_line(aes(y=cijerkU),linetype="dashed")+
  geom_line(aes(y=cijerkL),linetype="dashed")+
  geom_jitter(data=dt,aes(y=deltaJ),alpha=I(0.1),color="gray0")+
  geom_hline(yintercept=0,color="red")+
  xlim(c(0,100))+
  facet_wrap(.~HT)
ggsave("MRjerkNW.pdf",width=10,height=5.5,dpi=600)

# depth MEAN
ggplot(aggdepth,aes(hours,depthm))+
  geom_line()+
  geom_line(aes(y=cidepthU),linetype="dashed")+
  geom_line(aes(y=cidepthL),linetype="dashed")+
  geom_jitter(alpha=I(0.2),color="gray0")+
  geom_hline(yintercept=0,color="red")+
  xlim(c(0,100))+
  facet_wrap(.~HT)
ggsave("MRdepthmNW.pdf",width=10,height=5.5,dpi=600)

# depth MAX
ggplot(aggdepth,aes(hours,depthmax))+
  geom_line()+
  geom_line(aes(y=cidepthmU),linetype="dashed")+
  geom_line(aes(y=cidepthmL),linetype="dashed")+
  geom_jitter(alpha=I(0.2),color="gray0")+
  geom_hline(yintercept=0,color="red")+
  xlim(c(0,100))+
  facet_wrap(.~HT)
ggsave("MRdepthmaxNW.pdf",width=10,height=5.5,dpi=600)

# VedBA
aggboth[aggboth$HT=="Long handlingtime",]$vedbam<0 # 34
aggboth[aggboth$HT=="Long handlingtime",]$civedbaL<0 # 37 
aggboth[aggboth$HT=="Long handlingtime",]$civedbaU<0 # 7

aggboth[aggboth$HT=="Short handlingtime",]$vedbam<0 # 17
aggboth[aggboth$HT=="Short handlingtime",]$civedbaL<0 # 64 
aggboth[aggboth$HT=="Short handlingtime",]$civedbaU<0 # 0
# Jerk
aggboth[aggboth$HT=="Long handlingtime",]$jerkm<0 # 7
aggboth[aggboth$HT=="Long handlingtime",]$cijerkL<0 # 37
aggboth[aggboth$HT=="Long handlingtime",]$cijerkU<0 # 0

aggboth[aggboth$HT=="Short handlingtime",]$jerkm<0 # 15
aggboth[aggboth$HT=="Short handlingtime",]$cijerkL<0 # 32
aggboth[aggboth$HT=="Short handlingtime",]$cijerkU<0 # 0
