labels = c(6,21,120,1000,2000,4000)
old = c(36,126,579,5400,5400*2,5400*4)
new = c(4.5,5,7.3,59,432,4320)

x = barplot(rbind(old,new),beside=T,log="y",
        names=labels,xlab="x1000 generations",
        ylab="Seconds per Run")

## calculate simulations finished per hour
hours=(0:4800)
old_h = old/(60*60)
new_h = new/(60*60)

old_per = 1/old_h
new_per = 1/new_h

old_all = old_h*800
new_all = new_h*800

df_o = NULL
df_n = NULL

for(hr in hours){
  row_o = floor(old_per * hr)
  row_n = floor(new_per * hr)
  
  row_o[row_o >=800] = 800
  row_n[row_n >=800] = 800
  
  if(is.null(df_o)) {
    df_o = row_o
    df_n = row_n
  } else {
    df_o = rbind(df_o,row_o)
    df_n = rbind(df_n,row_n)
  }
  
}

colnames(df_o) = labels
rownames(df_o) = hours
head(df_o)

df_o = unique(df_o)
old_completion=rowSums(df_o)
newhours_o = rownames(df_o)


colnames(df_n) = labels
rownames(df_n) = hours
head(df_n)
df_n = unique(df_n)

new_completion=rowSums(df_n)
newhours_n = rownames(df_n)

old_comp_small = rowSums(df_o[,1:3])
new_comp_small = rowSums(df_n[,1:3])

plot(as.numeric(newhours_o)*10,(old_completion)*10,type="l",
     ylab="Simulations Finished",xlim=c(0,1300),
     ylim=c(0,25000),
     xlab="Hours Needed")
points(as.numeric(newhours_n)*10,(new_completion)*10,col="red",type="l")
abline(v=129*10,lty=3)
abline(v=2*10,lty=3)
