## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Functions for the Soay sheep coat colour IBM
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

library(arm) ## invlogit

initState <- function (m.par, init.pop.size = 500, maxA = 15) {
    ## Initialise individual state array for the population.

    ## Args:
    ##   init.pop.size: Initial size of the population.
    ##   QUERY any other initialisation choices - age, sex, genotype?

    ## Returns:
    ##   Individual state array containing the size of each individual
    ##   in the population, structured by age, sex and coat colour
    ##   genotype.

    ## Population structure - age, sex, genotype
    Aset <- as.character(seq.int(0, maxA))
    Gset <- c("GG","GT","TT")
    Sset <- c("F","M")

    sets <- list(S=Sset, A=Aset, G=Gset)

    ## Initialise state array (uses util function from IPM code)
    z <- mk.flist(list(S=Sset, A=Aset, G=Gset))

    ## Following strategy in supplementary material to Rees et
    ## al. 2014, which is to get a random selection of new recruits
    ## from normal distribution based on mean parent size:

    ## z <- rnorm(init.pop.size, mean = m.par["rcsz.int"]
    ##            + m.par["rcsz.z"] * 3.2, sd = m.par["rcsz.sd"])

    ## Offspring size model developed for the Soay coat colour IPM
    ## takes into account sex, twin status, maternal age, population
    ## density and observation year:
    ## capWgt~1+sex+capWgtMum+isTwnMat+ageMum*obsY+Ntm1+(1|obsYf)

    ## With centred standardized density and year this reduces to:
    ## mu = mPar["sz.off.(Intercept)"]+mPar["sz.off.capWgtMum"]*x+
    ## mPar["sz.off.ageMum"]*A+mPar["sz.off.sexM"]+mPar["sz.off.isTwnMat"]

    ## FIXME Need mean values mean.sz.F, mean.age.F for the following,
    ## these values are made up.
    mean.sz.F <- 3.2
    mean.age.F <- 5

    ## Assume equal numbers of offspring of each sex and genotype.
    n <- init.pop.size/6

    ## Get proportion of twins in the new population based on twinning rate
    p.twin <- invlogit(m.par["t.a1.F.(Intercept)"] + m.par["t.a1.F.capWgt"]*mean.sz.F
                       + m.par["t.a1.F.poly(ageY, 2, raw = TRUE)1"]*mean.age.F
                       + m.par["t.a1.F.poly(ageY, 2, raw = TRUE)2"]*mean.age.F^2)
    prop.twin <- 2*p.twin/(1+p.twin)

    for (G in Gset) {
        for (S in Sset) {
            mu <- m.par["sz.off.(Intercept)"] +
                m.par["sz.off.capWgtMum"]*mean.sz.F +
                    m.par["sz.off.ageMum"]*mean.age.F +
                        m.par["sz.off.sexM"]*as.numeric(S=="M")

            z[[S, "0", G]] <- c(rnorm(n*prop.twin
                                      , mean = mu + m.par["sz.off.isTwnMat"]
                                      , sd = m.par["sz.off.sigma"]),
                                rnorm(n*(1-prop.twin)
                                      , mean = mu
                                      , sd = m.par["sz.off.sigma"]))

        }
    }
    return(z)
}

doIBM <- function (model.params, init.pop.size = 500, sim.length = 200) {

    ## check passed in params
    ## - model params: should contain all expected parameters
    ## - initial population size
    ## - simulation length

    ## Get initial population state
    z <- initState (model.params, init.pop.size)

    for (i in seq_len(sim.length)) {
        doIBMStep(z)
    }

}
