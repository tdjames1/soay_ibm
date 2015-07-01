## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Demographic component functions
##
## x - size this autumn
## y - size next autumn
## S - sex
## A - age
## G - genotype
## Nt - standardised population density
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

mk.flist <- function(sets) {
    mydims <- sapply(sets,length)
    mylist <- vector("list", prod(mydims))
    dim(mylist) <- mydims
    dimnames(mylist) <- sets
    return(mylist)
}

## ~~~~~~~~~ SURVIVAL ~~~~~~~~~

p.surv.list <- mk.flist(list(S=c("F","M"), A=c(0,+1)))

p.surv.list[["F", "0"]] <- function(x, Nt, mPar) {
    nu <- mPar["s.a0.F.(Intercept)"]+mPar["s.a0.F.capWgt"]*x+
          mPar["s.a0.F.Nt"]*Nt+mPar["s.a0.F.obsY"]*mPar["obsY"]
    return(invlogit(nu))
}
p.surv.list[["F", "1"]] <- function(x, G, A, Nt, mPar, gSwitch) {
    prefix <- "s.a1.F"
    if (prefix %in% names(gSwitch) & !gSwitch[prefix]) { G = 0 }
    nu <- mPar["s.a1.F.(Intercept)"]+mPar["s.a1.F.capWgt"]*x+
          mPar["s.a1.F.APS071add"]*G+
          mPar["s.a1.F.poly(ageY, 2, raw = TRUE)1"]*A+mPar["s.a1.F.poly(ageY, 2, raw = TRUE)2"]*A^2+
          mPar["s.a1.F.Nt"]*Nt+mPar["s.a1.F.Nt:APS071add"]*G*Nt+mPar["s.a1.F.obsY"]*mPar["obsY"]
    return(invlogit(nu))
}
p.surv.list[["M", "0"]] <- function(x, G, Nt, mPar, gSwitch) {
    prefix <- "s.a0.M"
    if (prefix %in% names(gSwitch) & !gSwitch[prefix]) { G = 0 }
    nu <- mPar["s.a0.M.(Intercept)"]+mPar["s.a0.M.capWgt"]*x+
          mPar["s.a0.M.Nt"]*Nt+mPar["s.a0.M.capWgt:Nt"]*x*Nt+mPar["s.a0.M.obsY"]*mPar["obsY"]+
          mPar["s.a0.M.APS071add"]*G
    return(invlogit(nu))
}
p.surv.list[["M", "1"]] <- function(x, G, A, Nt, mPar, gSwitch) {
    prefix <- "s.a1.M"
    if (prefix %in% names(gSwitch) & !gSwitch[prefix]) { G = 0 }
    nu <- mPar["s.a1.M.(Intercept)"]+mPar["s.a1.M.capWgt"]*x+
          mPar["s.a1.M.APS071add"]*G+mPar["s.a1.M.ageY"]*A+
          mPar["s.a1.M.Nt"]*Nt+mPar["s.a1.M.obsY:APS071add"]*G*mPar["obsY"]+
          mPar["s.a1.M.obsY"]*mPar["obsY"]
    return(invlogit(nu))
}

p.surv <- function(x, G, S, A, Nt, mPar, gSwitch=numeric(0))
{
    ## expect S and A to have length=1
    if (length(S) != 1 | length(A) != 1)
        stop("survival function not vectorised for S(ex) and (A)ge")
    ## assign the required function
    if       (S=="F") {
        if        (A ==  0) {
            f <- p.surv.list[["F", "0"]]
        } else if (A >=  1) {
            f <- p.surv.list[["F", "1"]]
        } else stop("invalid (A)ge")
    } else if (S=="M") {
        if        (A ==  0) {
            f <- p.surv.list[["M", "0"]]
        } else if (A >=  1) {
            f <- p.surv.list[["M", "1"]]
        } else stop("invalid (A)ge")
    } else stop("invalid (S)ex")
    ## get the required arguments as a list of atomic 'names'
    fargs <- names(formals(f))
    fargs <- sapply(fargs, as.name)
    ## call the function and return
    return(do.call("f", fargs))
}

## ~~~~~~~~~ GROWTH ~~~~~~~~~

d.grow.list <- mk.flist(list(S=c("F","M"), A=c(0,+1)))

d.grow.list[["F","0"]] <- function(y, x, Nt, mPar) {
    mu <- mPar["g.a0.F.(Intercept)"]+mPar["g.a0.F.capWgt"]*x+mPar["g.a0.F.Nt"]*Nt
    sg <- mPar["g.a0.F.sigma"]
    return(dnorm(y,mu,sg))
}
d.grow.list[["F","1"]] <- function(y, x, G, A, Nt, mPar) {
    mu <- mPar["g.a1.F.(Intercept)"] + mPar["g.a1.F.capWgt"]*x + mPar["g.a1.F.poly(ageY, 2, raw = TRUE)1"]*A +
          mPar["g.a1.F.Nt"]*Nt + mPar["g.a1.F.ageY:Nt"]*Nt*A + mPar["g.a1.F.poly(ageY, 2, raw = TRUE)2"]*A^2
    sg <- mPar["g.a1.F.sigma"]
    return(dnorm(y,mu,sg))
}
d.grow.list[["M","0"]] <- function(y, x, Nt, mPar) {
    mu <- mPar["g.a0.M.(Intercept)"]+mPar["g.a0.M.capWgt"]*x+mPar["g.a0.M.Nt"]*Nt+mPar["g.a0.M.obsY"]*mPar["obsY"]
    sg <- mPar["g.a0.M.sigma"]
    return(dnorm(y,mu,sg))
}
d.grow.list[["M","1"]] <- function(y, x, G, A, Nt, mPar) {
    mu <- mPar["g.a1.M.(Intercept)"] + mPar["g.a1.M.capWgt"]*x
    sg <- mPar["g.a1.M.sigma"]
    return(dnorm(y,mu,sg))
}

d.grow <- function(y, x, G, S, A, Nt, mPar)
{
    ## expect S and A to have length=1
    if (length(S) != 1 | length(A) != 1)
        stop("growth function not vectorised for S(ex) and (A)ge")
    ## assign the required function
    if       (S=="F") {
        if        (A ==  0) {
            f <- d.grow.list[["F", "0"]]
        } else if (A >=  1) {
            f <- d.grow.list[["F", "1"]]
        } else stop("invalid (A)ge")
    } else if (S=="M") {
        if        (A ==  0) {
            f <- d.grow.list[["M", "0"]]
        } else if (A >=  1) {
            f <- d.grow.list[["M", "1"]]
        } else stop("invalid (A)ge")
    } else stop("invalid (S)ex")
    ## get the required arguments as a list of atomic 'names'
    fargs <- names(formals(f))
    fargs <- sapply(fargs, as.name)
    ## call the function and return
    return(do.call("f", fargs))
}

## ~~~~~~~~~ FEMALE PROBABILITY OF REPRODUCTION ~~~~~~~~~

p.repr.list <- mk.flist(list(A=c(0,+1)))

p.repr.list[["0"]] <- function(x, Nt, mPar) {
    nu <- mPar["r.a0.F.(Intercept)"]+mPar["r.a0.F.capWgt"]*x+
          mPar["r.a0.F.Nt"]*Nt+mPar["r.a0.F.obsY"]*mPar["obsY"]+
          mPar["r.a0.F.capWgt:obsY"]*x*mPar["obsY"]
    return(invlogit(nu))
}
p.repr.list[["1"]] <- function(A, Nt, mPar) {
    nu <- mPar["r.a1.F.(Intercept)"]+
          mPar["r.a1.F.poly(ageY, 2, raw = TRUE)1"]*A+mPar["r.a1.F.poly(ageY, 2, raw = TRUE)2"]*A^2+
          mPar["r.a1.F.obsY"]*mPar["obsY"]+mPar["r.a1.F.ageY:obsY"]*A*mPar["obsY"]
    return(invlogit(nu))
}

p.repr <- function(x, G, A, Nt, mPar)
{
    ## expect S and A to have length=1
    if (length(A) != 1)
        stop("female P(reproduction) function not vectorised for (A)ge")
    ## assign the required function
    if        (A ==  0) {
        f <- p.repr.list[["0"]]
    } else if (A >=  1) {
        f <- p.repr.list[["1"]]
    } else stop("invalid (A)ge")
    ## get the required arguments as a list of atomic 'names'
    fargs <- names(formals(f))
    fargs <- sapply(fargs, as.name)
    ## call the function and return
    return(do.call("f", fargs))
}

## ~~~~~~~~~ FEMALE PROBABILITY OF TWINNING ~~~~~~~~~

p.twin.list <- mk.flist(list(A=c(0,+1)))

p.twin.list[["0"]] <- function(x, Nt, mPar) {
    return(0)
}
p.twin.list[["1"]] <- function(x, A, Nt, mPar) {
    nu <- mPar["t.a1.F.(Intercept)"]+mPar["t.a1.F.capWgt"]*x+
          mPar["t.a1.F.poly(ageY, 2, raw = TRUE)1"]*A+mPar["t.a1.F.poly(ageY, 2, raw = TRUE)2"]*A^2+
          mPar["t.a1.F.Nt"]*Nt+mPar["t.a1.F.obsY"]*mPar["obsY"]
    return(invlogit(nu))
}

p.twin <- function(x, G, A, Nt, mPar)
{
    ## expect S and A to have length=1
    if (length(A) != 1)
        stop("female P(twinning) function not vectorised for (A)ge")
    ## assign the required function
    if        (A ==  0) {
        f <- p.twin.list[["0"]]
    } else if (A >=  1) {
        f <- p.twin.list[["1"]]
    } else stop("invalid (A)ge")
    ## get the required arguments as a list of atomic 'names'
    fargs <- names(formals(f))
    fargs <- sapply(fargs, as.name)
    ## call the function and return
    return(do.call("f", fargs))
}

## ~~~~~~~~~ OFFSPRING SPRING SURVIVAL FUNCTION ~~~~~~~~~

p.offsurv.list <- mk.flist(list(A=c(0,+1), S.off=c("F","M"), T.off=c("0","1")))

p.offsurv.list[["0","F","0"]] <- function(x, Nt, mPar) {
    nu <- mPar["s.off.a0.(Intercept)"]+mPar["s.off.a0.capWgtMum"]*x+
          mPar["s.off.a0.Ntm1"]*Nt+mPar["s.off.a0.obsY"]*mPar["obsY"]
    return(invlogit(nu))
}
p.offsurv.list[["0","M","0"]] <- function(x, Nt, mPar) {
    nu <- mPar["s.off.a0.(Intercept)"]+mPar["s.off.a0.capWgtMum"]*x+
          mPar["s.off.a0.Ntm1"]*Nt+mPar["s.off.a0.obsY"]*mPar["obsY"]+mPar["s.off.a0.sexM"]
    return(invlogit(nu))
}
p.offsurv.list[["1","F","0"]] <- function(x, A, Nt, mPar, isTwn) {
    nu <- mPar["s.off.a1.(Intercept)"]+mPar["s.off.a1.capWgtMum"]*x+
          mPar["s.off.a1.Ntm1"]*Nt+mPar["s.off.a1.obsY"]*mPar["obsY"]
    return(invlogit(nu))
}
p.offsurv.list[["1","M","0"]] <- function(x, A, Nt, mPar, isTwn) {
    nu <- mPar["s.off.a1.(Intercept)"]+mPar["s.off.a1.capWgtMum"]*x+
          mPar["s.off.a1.Ntm1"]*Nt+mPar["s.off.a1.obsY"]*mPar["obsY"]+mPar["s.off.a1.sexM"]
    return(invlogit(nu))
}
p.offsurv.list[["1","F","1"]] <- function(x, A, Nt, mPar, isTwn) {
    nu <- mPar["s.off.a1.(Intercept)"]+mPar["s.off.a1.capWgtMum"]*x+
          mPar["s.off.a1.Ntm1"]*Nt+mPar["s.off.a1.obsY"]*mPar["obsY"]+mPar["s.off.a1.isTwnMat"]
    return(invlogit(nu))
}
p.offsurv.list[["1","M","1"]] <- function(x, A, Nt, mPar, isTwn) {
    nu <- mPar["s.off.a1.(Intercept)"]+mPar["s.off.a1.capWgtMum"]*x+
          mPar["s.off.a1.Ntm1"]*Nt+mPar["s.off.a1.obsY"]*mPar["obsY"]+
          mPar["s.off.a1.sexM"]+mPar["s.off.a1.isTwnMat"]
    return(invlogit(nu))
}

p.offsurv <- function(x, G, A, Nt, mPar, S.off, T.off)
{
    ## expect S and A to have length=1
    if (length(A) != 1)
        stop("offspring survival function not vectorised for (A)ge")
    ## assign the required function
    if        (A == 0) {
        if        (S.off == "F") {
            if        (T.off == 0) {
                f <- p.offsurv.list[["0","F","0"]]
            } else if (T.off == 1) {
                return(0) # age 0 females cannot twin
            } else stop("invalid offspring (T)win status")
        } else if (S.off == "M") {
            if        (T.off == 0) {
                f <- p.offsurv.list[["0","M","0"]]
            } else if (T.off == 1) {
                return(0) # age 0 females cannot twin
            } else stop("invalid offspring (T)win status")
        } else stop("invalid offspring (S)ex")
    } else if (A >= 1) {
        if        (S.off == "F") {
            if        (T.off == 0) {
                f <- p.offsurv.list[["1","F","0"]]
            } else if (T.off == 1) {
                f <- p.offsurv.list[["1","F","1"]]
            } else stop("invalid offspring (T)win status")
        } else if (S.off == "M") {
            if        (T.off == 0) {
                f <- p.offsurv.list[["1","M","0"]]
            } else if (T.off == 1) {
                f <- p.offsurv.list[["1","M","1"]]
            } else stop("invalid offspring (T)win status")
        } else stop("invalid offspring (S)ex")
    } else stop("invalid (A)ge")
    ## get the required arguments as a list of atomic 'names'
    fargs <- names(formals(f))
    fargs <- sapply(fargs, as.name)
    ## call the function and return
    return(do.call("f", fargs))
}

## ~~~~~~~~~ OFFSPRING SIZE DISTRIBUTION ~~~~~~~~~

d.offsize.list <- mk.flist(list(S.off=c("F","M"), T.off=c("0","1")))

d.offsize.list[["F","0"]] <- function(y, x, A, Nt, mPar) {
    mu <- mPar["sz.off.(Intercept)"]+mPar["sz.off.capWgtMum"]*x+
          mPar["sz.off.ageMum"]*A+
          mPar["sz.off.Ntm1"]*Nt+
          mPar["sz.off.obsY"]*mPar["obsY"]+
          mPar["sz.off.ageMum:obsY"]*A*mPar["obsY"]
    sg <- mPar["sz.off.sigma"]
    return(dnorm(y,mu,sg))
}
d.offsize.list[["M","0"]] <- function(y, x, A, Nt, mPar) {
    mu <- mPar["sz.off.(Intercept)"]+mPar["sz.off.capWgtMum"]*x+
          mPar["sz.off.ageMum"]*A+
          mPar["sz.off.Ntm1"]*Nt+mPar["sz.off.sexM"]+
          mPar["sz.off.obsY"]*mPar["obsY"]+
          mPar["sz.off.ageMum:obsY"]*A*mPar["obsY"]

    sg <- mPar["sz.off.sigma"]
    return(dnorm(y,mu,sg))
}
d.offsize.list[["F","1"]] <- function(y, x, A, Nt, mPar) {
    mu <- mPar["sz.off.(Intercept)"]+mPar["sz.off.capWgtMum"]*x+
          mPar["sz.off.ageMum"]*A+
          mPar["sz.off.Ntm1"]*Nt+mPar["sz.off.isTwnMat"]+
          mPar["sz.off.obsY"]*mPar["obsY"]+
          mPar["sz.off.ageMum:obsY"]*A*mPar["obsY"]
    sg <- mPar["sz.off.sigma"]
    return(dnorm(y,mu,sg))
}
d.offsize.list[["M","1"]] <- function(y, x, A, Nt, mPar) {
    mu <- mPar["sz.off.(Intercept)"]+mPar["sz.off.capWgtMum"]*x+
          mPar["sz.off.ageMum"]*A+
          mPar["sz.off.Ntm1"]*Nt+mPar["sz.off.sexM"]+mPar["sz.off.isTwnMat"]+
          mPar["sz.off.obsY"]*mPar["obsY"]+
          mPar["sz.off.ageMum:obsY"]*A*mPar["obsY"]
    sg <- mPar["sz.off.sigma"]
    return(dnorm(y,mu,sg))
}

d.offsize <- function(y, x, G, A, Nt, mPar, S.off, T.off)
{
    ## expect S and A to have length=1
    if (length(S.off) != 1 | length(T.off) != 1)
        stop("offspring size function not vectorised for offspring S(ex), female (A)ge and (T)win status")
    ## assign the required function
    if (all(A >= 0)) {
        if        (S.off == "F") {
            if        (T.off == 0) {
                f <- d.offsize.list[["F","0"]]
            } else if (T.off == 1) {
                f <- d.offsize.list[["F","1"]]
            } else stop("invalid offspring (T)win status")
        } else if (S.off == "M") {
            if        (T.off == 0) {
                f <- d.offsize.list[["M","0"]]
            } else if (T.off == 1) {
                f <- d.offsize.list[["M","1"]]
            } else stop("invalid offspring (T)win status")
        } else stop("invalid offspring (S)ex")
    } else stop("invalid (A)ge")
    ## get the required arguments as a list of atomic 'names'
    fargs <- names(formals(f))
    fargs <- sapply(fargs, as.name)
    ## call the function and return
    return(do.call("f", fargs))
}

## ~~~~~~~~~ MALE MATING SUCCESS ~~~~~~~~~

n.mrepro.list <- mk.flist(list(A=c("0","1")))

n.mrepro.list[["0"]] <- function(x, G, Nt, mPar, gSwitch) {
    prefix <- "n.off.a0.M"
    if (prefix %in% names(gSwitch) & !gSwitch[prefix]) { G = 0 }
    nu <- mPar["n.off.a0.M.(Intercept)"] + mPar["n.off.a0.M.capWgt"]*x +
          mPar["n.off.a0.M.APS071add"]*G + mPar["n.off.a0.M.Nt"]*Nt + mPar["n.off.a0.M.obsY"]*mPar["obsY"]
    return(exp(nu))
}
n.mrepro.list[["1"]] <- function(x, G, Nt, mPar, gSwitch) {
    prefix <- "n.off.a1.M"
    if (prefix %in% names(gSwitch) & !gSwitch[prefix]) { G = 0 }
    nu <- mPar["n.off.a1.M.(Intercept)"] + mPar["n.off.a1.M.capWgt"]*x +
          mPar["n.off.a1.M.APS071add"]*G + mPar["n.off.a1.M.APS071add:Nt"]*G*Nt +
          mPar["n.off.a1.M.Nt"]*Nt + mPar["n.off.a1.M.obsY"]*mPar["obsY"]
    return(exp(nu))
}

n.mrepro <- function(x, G, A, Nt, mPar, gSwitch=numeric(0))
{
    ## expect S and A to have length=1
    if (length(A) != 1)
        stop("male reproduction function not vectorised for (A)ge")
    ## assign the required function
    if        (A ==  0) {
        f <- n.mrepro.list[["0"]]
    } else if (A >=  1) {
        f <- n.mrepro.list[["1"]]
    } else stop("invalid (A)ge")
    ## get the required arguments as a list of atomic 'names'
    fargs <- names(formals(f))
    fargs <- sapply(fargs, as.name)
    ## call the function and return
    return(do.call("f", fargs))
}
