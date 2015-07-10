set.seed(23020306)

## Simple trial run. Initial population equally distributed between
## genotypes. Constant environment model parameters.
sim.len <- 50
init.pop <- 500
sim.out <- doSim(mParFixEf, sim.length=sim.len, init.pop.size=init.pop)
simRunSum <- summariseSimRun(sim.out)
pG <- (simRunSum$ntGG + 0.5*simRunSum$ntGT)/simRunSum$ntT
plot(pG,type="l")

## Simulation with G allele only
sim.out <- doSim(mParFixEf, sim.length=500, init.pop.size=init.pop, init.G = "GG")
simRunSum <- summariseSimRun(sim.out)
plot(simRunSum$ntT, type="l", ylim=c(0,600))
lines(simRunSum$ntF, col="red")
lines(simRunSum$ntM, col="blue")

## Simulation with T allele only
sim.out <- doSim(mParFixEf, sim.length=500, init.pop.size=init.pop, init.G = "TT")
simRunSum <- summariseSimRun(sim.out)
plot(simRunSum$ntT, type="l", ylim=c(0,600))
lines(simRunSum$ntF, col="red")
lines(simRunSum$ntM, col="blue")

## Simulate invasion of superior into inferior genotype
numReps <- 10000
resG <- "TT"
mPar <- list(C=mParFixEf, S=mParStore)
invProp <- c(0.01, 0.02, 0.05)
invRate <- list()

for (m in seq_along(mPar)) {
    for (i in invProp) {
        label <- paste0(i,names(mPar)[[m]])
        invRate[[label]] <- simInv(mPar[[m]], inv.prop=i, num.reps=numReps, res=resG)

        print(paste(switch(names(mPar)[[m]], C="Constant", S="Stochastic"),
                    "environment, proportion of invader:", i))

        ## Mean invader log growth rate
        print(summary(invRate[[label]], na.rm=TRUE))

        ## Proportion positive vs negative growth rate
        print(paste("Invasion rate positive", sum(invRate[[label]]>0)/numReps*100, "%"))
    }
}
