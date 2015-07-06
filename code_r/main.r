rm(list=ls(all=TRUE))

if (Sys.getenv("LOGNAME") == "Tamora" || Sys.getenv("LOGNAME") == "tamorajames") {
    setwd("~/Projects/MRes/soay_ibm")
}

source("./code_r/load.r")
source("./code_r/imported/demog_fun.r")
source("./code_r/imported/popgen_fun.r")
source("./code_r/ibm_fun.r")

set.seed(23020306)

sim.len <- 250
init.pop <- 700
sim.out <- doSim(mParFixEf, sim.length=sim.len, init.pop.size=init.pop)

simRunSum <- summariseSimRun(sim.out)

pG <- (simRunSum$ntGG + 0.5*simRunSum$ntGT)/simRunSum$ntT
plot(pG,type="l")
