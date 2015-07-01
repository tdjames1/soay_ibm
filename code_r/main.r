rm(list=ls(all=TRUE))

if (Sys.getenv("LOGNAME") == "Tamora" || Sys.getenv("LOGNAME") == "tamorajames") {
    setwd("~/Projects/MRes/soay_ibm")
}

source("./code_r/load.r")
source("./code_r/imported/demog_fun.r")
source("./code_r/imported/popgen_fun.r")
source("./code_r/ibm_fun.r")

set.seed(23020306)

sim.len <- 50
z <- doSim(mParFixEf, sim.length=sim.len)
