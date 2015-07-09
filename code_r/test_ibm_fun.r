## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Functions for testing the Soay sheep coat colour IBM
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

getInvaderGrRt <- function(mPar, resG, resT=50, invT=50, init.pop=500, inv.prop=0.01) {
    ## Calculate average marginal growth rate for an invading allele
    ## introduced at low frequency into an existing population
    ## initially consisting of individuals of a single genotype.

    ## Args:
    ##   mPar: Model parameters. Numeric vector or matrix.
    ##   resG: Resident genotype.
    ##   resT: Number of timesteps to iterate with resident only.
    ##   invT: Number of timesteps to iterate with invader.
    ##   init.pop: Initial population size.
    ##   inv.prop: Invader proportion.

    ## Returns:
    ##   Mean of log invader growth rate.
    invG <- switch(resG, GG="TT", TT="GG")
    if (is.null(invG)) {
        stop (paste ("Invalid resident genotype:", resG))
    }

    ## Simulate resident population
    sim.res <- doSim(mPar, sim.length=resT, init.pop.size=init.pop, init.G = resG)
    res.z <- sim.res[[resT]]

    ## Introduce invader - sample sizes from resident population, assign to random sex/age
    invN <- inv.prop*sum(unlist(res.z))
    validA <- sapply(dimnames(res.z)$A, function(x) !is.null(z[["F", x, resG]]))
    invA <- sample(dimnames(res.z)$A[validA], invN, replace=TRUE)
    invS <- sample(dimnames(res.z)$S, invN, replace=TRUE)
    invZ <- sample(unlist(res.z), invN, replace=TRUE)
    for (i in seq_len(invN)) {
        res.z[[invS[i], invA[i], invG]] <- invZ[i]
    }
    sim.inv <- doSim(mPar, sim.length=invT, init.state=res.z)
    simRunSum <- summariseSimRun(sim.inv)
    invNt <- with(simRunSum, switch(invG, GG=ntGG+ntGT/2, TT=ntTT+ntGT/2))

    ## Return arithmetic mean of invader growth rate
    ## N.B. this returns NaN if the invasion fails to establish at first step
    return(mean(diff(log(invNt[invNt>0]))))
}
