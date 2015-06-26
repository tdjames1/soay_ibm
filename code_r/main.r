##rm(list=ls(all=TRUE))

if (Sys.getenv("LOGNAME") == "Tamora" || Sys.getenv("LOGNAME") == "tamorajames") {
    setwd("~/Projects/MRes/soay_ibm")
}

set.seed(23020306)

source("./code_r/load.r")
source("./code_r/ibm_fun.r")

sim.len <- 40
z <- doSim(mParFixEf, sim.length=sim.len)
