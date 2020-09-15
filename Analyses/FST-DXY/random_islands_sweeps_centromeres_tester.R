setseed(12042012)

empirical_sweeps=2690
empirical_islands=83
genome_prop_cent=0.15
sweep_prop_cent=0.304
island_prop_cent=0.304
sweep_reg_prop_cent = 0.278
island_reg_prop_cent = 0.235


run_bernoulli = function(number,proportion){
  sum(Rlab::rbern(number,proportion))/number
}

run_bernoulli(number=empirical_islands,proportion=genome_prop_cent)
run_bernoulli(number=empirical_sweeps,proportion=genome_prop_cent)

island=sapply(1:100000,FUN = function(x){run_bernoulli(number=empirical_islands,proportion=genome_prop_cent)})
sweep=sapply(1:100000,FUN = function(x){run_bernoulli(number=empirical_sweeps,proportion=genome_prop_cent)})

png("islands_sweeps_proportions_in_centromeres_simulated.png",
    height=300,
    width=600)
par(mar=c(4,4,0,0))
hist(island,xlim=c(0,0.35),breaks=seq(0,0.4,0.01),ylim=c(0,45000),col=rgb(0,1,1,0.3),
     main="",ylab="Count",xlab="Proportion in Centromeres")
hist(sweep,xlim=c(0,0.35),breaks=seq(0,0.4,0.01),ylim=c(0,45000),col=rgb(1,0,1,0.3),add=T)
abline(v=sweep_prop_cent,col="purple",lwd=2,lty=2)
abline(v=sweep_reg_prop_cent,col="magenta",lwd=2,lty=3)
abline(v=island_reg_prop_cent,col="cyan",lwd=2,lty=3)
dev.off()

island_reg_above=sum(island>=island_reg_prop_cent)/100000 ## 1894
island_above=sum(island>=island_prop_cent)/100000 ## 5
sweep_reg_above=sum(sweep>=sweep_reg_prop_cent)/100000 ## 0
sweep_above=sum(sweep>=sweep_prop_cent)/100000 ## 0 

