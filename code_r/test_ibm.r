source("./code_r/test_ibm_fun.r")

set.seed(23020306)

sim.len <- 50
init.pop <- 500
sim.out <- doSim(mParFixEf, sim.length=sim.len, init.pop.size=init.pop)

simRunSum <- summariseSimRun(sim.out)

pG <- (simRunSum$ntGG + 0.5*simRunSum$ntGT)/simRunSum$ntT
plot(pG,type="l")

sim.out <- doSim(mParFixEf, sim.length=500, init.pop.size=init.pop, init.G = "GG")
simRunSum <- summariseSimRun(sim.out)
plot(simRunSum$ntT, type="l")
lines(simRunSum$ntF, col="red")
lines(simRunSum$ntM, col="blue")

sim.out <- doSim(mParFixEf, sim.length=500, init.pop.size=init.pop, init.G = "TT")
simRunSum <- summariseSimRun(sim.out)
plot(simRunSum$ntT, type="l")
lines(simRunSum$ntF, col="red")
lines(simRunSum$ntM, col="blue")

resG <- "TT"
sim.out <- doSim(mParFixEf, sim.length=500, init.pop.size=init.pop, init.G = resG)
z <- sim.out[[500]]

## introduce invader
inv.prop <- 0.05
invN <- inv.prop*sum(unlist(z))
inv <- sample(unlist(z), invN, replace=TRUE)
invS <- sample(dimnames(z)$S, invN, replace=TRUE)
index <- integer(0)
Arange <- dimnames(z)$A
for (i in seq_along(Arange)) {
    if (!is.null(z[["F",Arange[i],resG]])) {
        index <- c(index,i)
    }
}
invA <- sample(Arange[index], invN, replace=TRUE)
for (i in seq_along(inv)) {
    z[[invS[i],invA[i],switch(resG, GG="TT", TT="GG")]] <- inv[i]
}

sim.inv <- doSim(mParFixEf, sim.length=500, init.state=z)
simRunSum <- summariseSimRun(sim.inv)

pG <- (simRunSum$ntGG + 0.5*simRunSum$ntGT)/simRunSum$ntT
plot(pG,type="l")

n <- 1000
invRate <- numeric(n)
for (i in seq_len(n)) {
    invRate[i] <- getInvaderGrRt(mParFixEf, "TT")
}

getInvaderGrRt(mParStore, "TT")
