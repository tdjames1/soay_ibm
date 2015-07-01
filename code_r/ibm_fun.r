## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Functions for the Soay sheep coat colour IBM
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

library(arm) ## invlogit

## Population statistics from annual sheep data
sdNt <- 127.1302
midNt <- 431.7917
scaleN <- function(N) return((N-midNt)/sdNt)

## Mean age and capture weight of mothers of age 0 individuals
mean.ageMum <- 4.06147
mean.capWgtMum <- 3.082146

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

    ## Initialise state array
    z <- mk.flist(list(S=Sset, A=Aset, G=Gset))

    ## Following strategy in supplementary material to Rees et
    ## al. 2014, which is to get a random selection of new recruits
    ## from normal distribution based on mean parent size.  Offspring
    ## size model developed for the Soay coat colour IPM takes into
    ## account sex, twin status, maternal age, population density and
    ## observation year.  With centred standardized density and year
    ## this reduces to:

    ## mu = mPar["sz.off.(Intercept)"]+mPar["sz.off.capWgtMum"]*x+
    ## mPar["sz.off.ageMum"]*A+mPar["sz.off.sexM"]+mPar["sz.off.isTwnMat"]

    ## Assume equal numbers of offspring of each sex and genotype.
    n <- init.pop.size/6

    ## Get proportion of twins in the new population based on twinning rate
    p.twin <- invlogit(m.par["t.a1.F.(Intercept)"] + m.par["t.a1.F.capWgt"]*mean.capWgtMum
                       + m.par["t.a1.F.poly(ageY, 2, raw = TRUE)1"]*mean.ageMum
                       + m.par["t.a1.F.poly(ageY, 2, raw = TRUE)2"]*mean.ageMum^2)
    prop.twin <- 2*p.twin/(1+p.twin)

    for (G in Gset) {
        for (S in Sset) {
            mu <- m.par["sz.off.(Intercept)"] +
                m.par["sz.off.capWgtMum"]*mean.capWgtMum +
                    m.par["sz.off.ageMum"]*mean.ageMum +
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

doSim <- function (model.params, init.pop.size = 500, sim.length = 200, maxA = 15) {

    ## check passed in params
    ## - model params: should contain all expected parameters
    ## - initial population size
    ## - simulation length
    ## - misc params such as maxA etc

    ## this could be generalised to require that client passes in the
    ## functions required to set up the state array for the system and
    ## to iterate the array between time steps.

    ## Check the dimensions of model.params to determine whether
    ## to run simulation as averaged or variable environment.
    tEnv <- numeric(0)
    mPar <- model.params
    if (!is.null(dim(mPar))) {
        ## Sample environment parameterss from the set available
        tEnv <- sample(1:(dim(mPar)[2]), sim.length+1, replace=TRUE)
        mPar <- model.params[,tEnv[1]]
    }

    ## Get initial population state
    z <- initState (mPar, init.pop.size, maxA)

    ## Population structure - age, sex, genotype
    Aset <- dimnames(z)$A
    Gset <- dimnames(z)$G
    Sset <- dimnames(z)$S

    numAset <- as.numeric(Aset)
    names(numAset) <- Aset
    numGset <- seq.int(-1,1)
    names(numGset) <- Gset

    ## Storage for population size and mean body mass
    num <- num.F <- num.M <- num.GG <- num.GT <- num.TT <-
        zbar <- zbar.F <- zbar.M <- zbar.GG <- zbar.GT <- zbar.TT <- numeric(sim.length)

    for (i in seq_len(sim.length)) {
        if(length(tEnv)>0) {
            ## Select model parameters for this step
            mPar <- model.params[,tEnv[i+1]]
        }

        ## Get the current population density (standardized according
        ## to the mean and standard deviation observed in the St Kilda
        ## population)
        Nt <- scaleN(length(unlist(z)))

        ## Calculate male paternity probabilities for each genotype
        pGm <- calcMalePaternity(z, Nt, mPar)

        ## Apply life history events
        z1 <- mk.flist(list(S=Sset, A=Aset, G=Gset))
        z1.rec <- mk.flist(list(S=Sset, G=Gset))
        for (G in Gset) {
            numG <- switch(G, GG=-1, GT=0, TT=+1)
            for (A in Aset) {
                numA <- numAset[A]
                for (S in Sset) {
                    z.sub <- z[[S, A, G]]
                    if (!is.null(z.sub) & length(z.sub) > 0) {

                        ## survival
                        pS <- p.surv(z.sub, numG, S, numA, Nt, mPar)
                        surv <- rbinom(n=length(z.sub), prob=pS, size=1)
                        z.sub <- z.sub[which(surv == 1)]

                        if (length(z.sub) == 0) next

                        ## growth
                        if (numA < maxA) {
                            z1[[S, as.character(numA+1), G]] <- r.grow(z.sub, numG, S,
                                                                       numA, Nt, mPar)
                        }

                        ## reproduction
                        if (S == "F") {
                            pB <- p.repr(z.sub, numG, numA, Nt, mPar)
                            repr <- rbinom(n=length(z.sub), prob=pB, size=1)
                            z.repr <- z.sub[which(repr == 1)]

                            ## Twin or not?
                            pT <- p.twin(z.repr, numG, numA, Nt, mPar)
                            twin <- rbinom(n=length(z.repr), prob=pT, size=1)

                            for (t.off in seq.int(0,1)) {

                                z_ <- z.repr[which(twin == t.off)]
                                n.off <- ifelse(t.off == 1, 2*length(z_), length(z_))

                                ## Offspring sex
                                o.sex <- rbinom(n=n.off, prob=0.5, size=1)
                                o.sex <- paste(c("F","M")[match(o.sex, c(0,1))])

                                ## Offspring genotype
                                ## pGm is paternity probability for male genotypes
                                pOffgen <- p.offgenotype(Gset, G, pGm)
                                ## Assign (numeric) genotypes to
                                ## offspring; not taking account of
                                ## paternity of twins for now. (Pemberton
                                ## et al. 1999 recorded 26% of twins with
                                ## known paternity having same sire.)
                                o.gen <- rmultinom(n=n.off, prob=pOffgen, size=1)
                                o.gen <- t(o.gen)%*%numGset

                                ## Offspring survival
                                ## Need to do for each offspring sex (M/F) and twin (0/1) status
                                for (s.off in Sset) {
                                    ## Index into maternal set for this sex
                                    ## This should sort out the parents appropriately for twin offspring
                                    i.subset <- which(o.sex == s.off)
                                    i.par <- ifelse(t.off == 1, ceiling(i.subset/2), i.subset)
                                    pRec <- p.offsurv(z_[i.par], A=numA, Nt=Nt, mPar=mPar,
                                                      S.off=s.off, T.off=t.off)
                                    recr <- rep(NA, n.off)
                                    recr[i.subset] <- rbinom(n=length(i.subset), prob=pRec, size=1)

                                    ## Recruit size
                                    if (sum(recr,na.rm=TRUE)>0) {
                                        z1.rec <- rep(NA, n.off)
                                        i.recr <- which(recr == 1)
                                        i.par <- ifelse(t.off == 1, ceiling(i.recr/2), i.recr)
                                        z1.rec[i.recr] <- r.offsize(z_[i.par], numG, numA, Nt,
                                                                    mPar, S.o=s.off, T.off=t.off)
                                        ## Store result for each genotype from the survivors
                                        for (numG.off in unique(o.gen[which(recr==1)])) {
                                            G.off <- Gset[numG.off+2]
                                            z1[[s.off, "0", G.off]] <- c(z1[[s.off, "0", G.off]],
                                                                         z1.rec[which(recr == 1 &
                                                                                      o.gen == numG.off)])
                                        }
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }
        ## Store some population stats - number of individuals, mean body mass for different classes
        num[i] <- length(unlist(z[   ,,    ]))
        num.F[i] <- length(unlist(z["F",,    ]))
        num.M[i] <- length(unlist(z["M",,    ]))
        num.GG[i] <- length(unlist(z[   ,,"GG"]))
        num.GT[i] <- length(unlist(z[   ,,"GT"]))
        num.TT[i] <- length(unlist(z[   ,,"TT"]))
        zbar[i] <- mean(unlist(z[   ,,    ]))
        zbar.F[i] <-  mean(unlist(z["F",,    ]))
        zbar.M[i] <-  mean(unlist(z["M",,    ]))
        zbar.GG[i] <-  mean(unlist(z[   ,,"GG"]))
        zbar.GT[i] <-  mean(unlist(z[   ,,"GT"]))
        zbar.TT[i] <-  mean(unlist(z[   ,,"TT"]))
        z <- z1
    }
    return(list(z=z,
                pop.size=list(T=num,F=num.F,M=num.M,GG=num.GG,GT=num.GT,TT=num.TT),
                mean.z=list(T=zbar,F=zbar.F,M=zbar.M,GG=zbar.GG,GT=zbar.GT,TT=zbar.TT)))
}

calcMalePaternity <- function(z, Nt, model.params) {
    ## Calculate male paternity probabilities for each genotype.

    ## Args:
    ##   z: Population state vector
    ##   Nt: Population density
    ##   model.params: Vital rate model parameters

    ## Returns:
    ##   Vector indicating the probability of paternity for males of each genotype.
    Gset <- dimnames(z)$G
    Aset <- dimnames(z)$A
    num.off <- mk.flist(list(G=Gset))
    for (G in Gset) {
        num.off[[G]] <- 0
        numG <- switch(G, GG=-1, GT=0, TT=+1)
        for (A in Aset) {
            z.sub <- z[["M", A, G]]
            if (!is.null(z.sub) & length(z.sub) > 0) {
                numA <- as.numeric(A)
                ## n.mrepro gives mean num offspring for each individual
                lambda <- n.mrepro(z.sub, numG, numA, Nt, model.params)
                ## rpois gets random pick from Poisson distribution with the given mean
                num.off[[G]] <- sum(sapply(lambda, function(x) rpois(1,x))) + num.off[[G]]
            }
        }
    }
    ## obtain mating probs
    nTot <- sum(unlist(num.off))
    pGm <- sapply(num.off, function(x) x/nTot)
    return(pGm)
}

## Adapted from IPM code
r.grow.list <- mk.flist(list(S=c("F","M"), A=c(0,+1)))

r.grow.list[["F","0"]] <- function(x, Nt, mPar) {
    mu <- mPar["g.a0.F.(Intercept)"]+mPar["g.a0.F.capWgt"]*x+mPar["g.a0.F.Nt"]*Nt
    sg <- mPar["g.a0.F.sigma"]
    return(rnorm(length(x),mu,sg))
}
r.grow.list[["F","1"]] <- function(x, G, A, Nt, mPar) {
    mu <- mPar["g.a1.F.(Intercept)"] +
        mPar["g.a1.F.capWgt"]*x +
            mPar["g.a1.F.poly(ageY, 2, raw = TRUE)1"]*A +
                mPar["g.a1.F.Nt"]*Nt +
                    mPar["g.a1.F.ageY:Nt"]*Nt*A +
                        mPar["g.a1.F.poly(ageY, 2, raw = TRUE)2"]*A^2
    sg <- mPar["g.a1.F.sigma"]
    return(rnorm(length(x),mu,sg))
}
r.grow.list[["M","0"]] <- function(x, Nt, mPar) {
    mu <- mPar["g.a0.M.(Intercept)"]+
        mPar["g.a0.M.capWgt"]*x+
            mPar["g.a0.M.Nt"]*Nt+
                mPar["g.a0.M.obsY"]*mPar["obsY"]
    sg <- mPar["g.a0.M.sigma"]
    return(rnorm(length(x),mu,sg))
}
r.grow.list[["M","1"]] <- function(x, G, A, Nt, mPar) {
    mu <- mPar["g.a1.M.(Intercept)"] +
        mPar["g.a1.M.capWgt"]*x
    sg <- mPar["g.a1.M.sigma"]
    return(rnorm(length(x),mu,sg))
}

r.grow <- function(x, G, S, A, Nt, mPar)
{
    ## expect S and A to have length=1
    if (length(S) != 1 | length(A) != 1)
        stop("growth function not vectorised for S(ex) and (A)ge")
    ## assign the required function
    if       (S=="F") {
        if        (A ==  0) {
            f <- r.grow.list[["F", "0"]]
        } else if (A >=  1) {
            f <- r.grow.list[["F", "1"]]
        } else stop("invalid (A)ge")
    } else if (S=="M") {
        if        (A ==  0) {
            f <- r.grow.list[["M", "0"]]
        } else if (A >=  1) {
            f <- r.grow.list[["M", "1"]]
        } else stop("invalid (A)ge")
    } else stop("invalid (S)ex")
    ## get the required arguments as a list of atomic 'names'
    fargs <- names(formals(f))
    fargs <- sapply(fargs, as.name)
    ## call the function and return
    return(do.call("f", fargs))
}

## Adapted from IPM code
r.offsize.list <- mk.flist(list(S.off=c("F","M"), T.off=c("0","1")))

r.offsize.list[["F","0"]] <- function(x, A, Nt, mPar) {
    mu <- mPar["sz.off.(Intercept)"]+mPar["sz.off.capWgtMum"]*x+
          mPar["sz.off.ageMum"]*A+
          mPar["sz.off.Ntm1"]*Nt+
          mPar["sz.off.obsY"]*mPar["obsY"]+
          mPar["sz.off.ageMum:obsY"]*A*mPar["obsY"]
    sg <- mPar["sz.off.sigma"]
    return(rnorm(length(x),mu,sg))
}
r.offsize.list[["M","0"]] <- function(x, A, Nt, mPar) {
    mu <- mPar["sz.off.(Intercept)"]+mPar["sz.off.capWgtMum"]*x+
          mPar["sz.off.ageMum"]*A+
          mPar["sz.off.Ntm1"]*Nt+mPar["sz.off.sexM"]+
          mPar["sz.off.obsY"]*mPar["obsY"]+
          mPar["sz.off.ageMum:obsY"]*A*mPar["obsY"]
    sg <- mPar["sz.off.sigma"]
    return(rnorm(length(x),mu,sg))
}
r.offsize.list[["F","1"]] <- function(x, A, Nt, mPar) {
    mu <- mPar["sz.off.(Intercept)"]+mPar["sz.off.capWgtMum"]*x+
          mPar["sz.off.ageMum"]*A+
          mPar["sz.off.Ntm1"]*Nt+mPar["sz.off.isTwnMat"]+
          mPar["sz.off.obsY"]*mPar["obsY"]+
          mPar["sz.off.ageMum:obsY"]*A*mPar["obsY"]
    sg <- mPar["sz.off.sigma"]
    return(rnorm(length(x),mu,sg))
}
r.offsize.list[["M","1"]] <- function(x, A, Nt, mPar) {
    mu <- mPar["sz.off.(Intercept)"]+mPar["sz.off.capWgtMum"]*x+
          mPar["sz.off.ageMum"]*A+
          mPar["sz.off.Ntm1"]*Nt+mPar["sz.off.sexM"]+mPar["sz.off.isTwnMat"]+
          mPar["sz.off.obsY"]*mPar["obsY"]+
          mPar["sz.off.ageMum:obsY"]*A*mPar["obsY"]
    sg <- mPar["sz.off.sigma"]
    return(rnorm(length(x),mu,sg))
}

r.offsize <- function(x, G, A, Nt, mPar, S.off, T.off)
{
    ## expect S and A to have length=1
    if (length(S.off) != 1 | length(T.off) != 1)
        stop("offspring size function not vectorised for offspring S(ex), female (A)ge and (T)win status")
    ## assign the required function
    if (all(A >= 0)) {
        if        (S.off == "F") {
            if        (T.off == 0) {
                f <- r.offsize.list[["F","0"]]
            } else if (T.off == 1) {
                f <- r.offsize.list[["F","1"]]
            } else stop("invalid offspring (T)win status")
        } else if (S.off == "M") {
            if        (T.off == 0) {
                f <- r.offsize.list[["M","0"]]
            } else if (T.off == 1) {
                f <- r.offsize.list[["M","1"]]
            } else stop("invalid offspring (T)win status")
        } else stop("invalid offspring (S)ex")
    } else stop("invalid (A)ge")
    ## get the required arguments as a list of atomic 'names'
    fargs <- names(formals(f))
    fargs <- sapply(fargs, as.name)
    ## call the function and return
    return(do.call("f", fargs))
}
