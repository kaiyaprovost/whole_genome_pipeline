df = read.table("/Users/kprovost/Dropbox (AMNH)/Dissertation/temp_sweep_isl_rank.txt",
                sep="\t",header=T,stringsAsFactors = F)
df = df[order(df$rank),]

df_nofstout = df[df$Species!="T. curvirostre",]

## fst 
{
plot(df$fst,df$Islands)
mod=lm(df$Islands~df$fst)
abline(mod,col="red")
summary(mod) ## not sig

plot(df$fst,df$island.site)
mod=lm(df$island.site~df$fst)
abline(mod,col="red")
summary(mod) ## not sig 

plot(df$fst,df$Islands..S.)
mod=lm(df$Islands..S.~df$fst)
abline(mod,col="red")
summary(mod) ## not sig 

plot(df$fst,df$islands.s.site)
mod=lm(df$islands.s.site~df$fst)
abline(mod,col="red")
summary(mod) ## not sig 


plot(df$fst,df$Sweeps)
mod=lm(df$Sweeps~df$fst)
abline(mod,col="red")
summary(mod) ## not sig but positive rsq = 0.15, arsq=0.04

plot(df$fst,df$sweep.site)
mod=lm(df$sweep.site~df$fst)
abline(mod,col="red")
summary(mod) ## not sig 

plot(df$fst,df$Sweeps..S.)
mod=lm(df$Sweeps..S.~df$fst)
abline(mod,col="red")
summary(mod) ## not sig  but pos rsq = 0.18 adj 0.07

plot(df$fst,df$sweeps.s.site)
mod=lm(df$sweeps.s.site~df$fst)
abline(mod,col="red")
summary(mod) ## not sig 
}

## rank 
{
  plot(df$rank,df$Islands)
  mod=lm(df$Islands~df$rank)
  abline(mod,col="red")
  summary(mod) ## not sig
  
  plot(df$rank,df$island.site)
  mod=lm(df$island.site~df$rank)
  abline(mod,col="red")
  summary(mod) ## not sig but pos rsq 0.11 adj 0.005
  
  plot(df$rank,df$Islands..S.)
  mod=lm(df$Islands..S.~df$rank)
  abline(mod,col="red")
  summary(mod) ## not sig 
  
  plot(df$rank,df$islands.s.site)
  mod=lm(df$islands.s.site~df$rank)
  abline(mod,col="red")
  summary(mod) ## not sig 
  
  
  plot(df$rank,df$Sweeps)
  mod=lm(df$Sweeps~df$rank)
  abline(mod,col="red")
  summary(mod) ## ALMOST sig but positive rsq = 0.32, arsq=0.2
  
  plot(df$rank,df$sweep.site)
  mod=lm(df$sweep.site~df$rank)
  abline(mod,col="red")
  summary(mod) ## not sig 
  
  plot(df$rank,df$Sweeps..S.)
  mod=lm(df$Sweeps..S.~df$rank)
  abline(mod,col="red")
  summary(mod) ## almost sig  but pos rsq = 0.36 adj 0.28
  
  plot(df$rank,df$sweeps.s.site)
  mod=lm(df$sweeps.s.site~df$rank)
  abline(mod,col="red")
  summary(mod) ## not sig but pos rsq = 0.19 adj = 0.09
}

## str 
{
  plot(df$str,df$Islands)
  mod=lm(df$Islands~df$str)
  abline(mod,col="red")
  summary(mod) ## not sig
  
  plot(df$str,df$island.site)
  mod=lm(df$island.site~df$str)
  abline(mod,col="red")
  summary(mod) ## not sig 
  
  plot(df$str,df$Islands..S.)
  mod=lm(df$Islands..S.~df$str)
  abline(mod,col="red")
  summary(mod) ## not sig 
  
  plot(df$str,df$islands.s.site)
  mod=lm(df$islands.s.site~df$str)
  abline(mod,col="red")
  summary(mod) ## not sig 
  
  
  plot(df$str,df$Sweeps)
  mod=lm(df$Sweeps~df$str)
  abline(mod,col="red")
  summary(mod) ## not sig but positive rsq = 0.15, arsq=0.04
  
  plot(df$str,df$sweep.site)
  mod=lm(df$sweep.site~df$str)
  abline(mod,col="red")
  summary(mod) ## not sig 
  
  plot(df$str,df$Sweeps..S.)
  mod=lm(df$Sweeps..S.~df$str)
  abline(mod,col="red")
  summary(mod) ## not sig  but pos rsq = 0.18 adj 0.07
  
  plot(df$str,df$sweeps.s.site)
  mod=lm(df$sweeps.s.site~df$str)
  abline(mod,col="red")
  summary(mod) ## not sig 
}


## Structure.across.CFB 
{
  boxplot(df$Islands~df$Structure.across.CFB)
  mod=lm(df$Islands~df$Structure.across.CFB)
  abline(mod,col="red")
  summary(mod) ## not sig
  
  boxplot(df$island.site~df$Structure.across.CFB)
  mod=lm(df$island.site~df$Structure.across.CFB)
  abline(mod,col="red")
  summary(mod) ## not sig 
  
  boxplot(df$Islands..S.~df$Structure.across.CFB)
  mod=lm(df$Islands..S.~df$Structure.across.CFB)
  abline(mod,col="red")
  summary(mod) ## not sig 
  
  boxplot(df$islands.s.site~df$Structure.across.CFB)
  mod=lm(df$islands.s.site~df$Structure.across.CFB)
  abline(mod,col="red")
  summary(mod) ## not sig 
  
  
  boxplot(df$Sweeps~df$Structure.across.CFB)
  mod=lm(df$Sweeps~df$Structure.across.CFB)
  abline(mod,col="red")
  summary(mod) ## not sig but positive rsq = 0.15, arsq=0.04
  
  boxplot(df$sweep.site~df$Structure.across.CFB)
  mod=lm(df$sweep.site~df$Structure.across.CFB)
  abline(mod,col="red")
  summary(mod) ## not sig 
  
  boxplot(df$Sweeps..S.~df$Structure.across.CFB)
  mod=lm(df$Sweeps..S.~df$Structure.across.CFB)
  abline(mod,col="red")
  summary(mod) ## not sig  but pos rsq = 0.18 adj 0.07
  
  boxplot(df$sweeps.s.site~df$Structure.across.CFB)
  mod=lm(df$sweeps.s.site~df$Structure.across.CFB)
  abline(mod,col="red")
  summary(mod) ## not sig 
}
par(mfrow=c(2,2))
barplot(df$island.site,names=df$Species,
        col=as.numeric(as.factor(df$str)),
        las=2)
barplot(df$islands.s.site,names=df$Species,
        col=as.numeric(as.factor(df$str)),
        las=2)
barplot(df$sweeps.s.site,names=df$Species,
        col=as.numeric(as.factor(df$str)),
        las=2)
barplot(df$sweep.site,names=df$Species,
        col=as.numeric(as.factor(df$str)),
        las=2)




