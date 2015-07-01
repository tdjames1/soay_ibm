## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## population genetics
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Gset <- c("GG","GT","TT")
G_table <- expand.grid(Gf=Gset,Gm=Gset,GG=NA,GT=NA,TT=NA)
G_table[1,Gset] <- c(1  ,   0,   0)
G_table[2,Gset] <- c(1/2, 1/2,   0)
G_table[3,Gset] <- c(0  ,   1,   0)
G_table[4,Gset] <- c(1/2, 1/2,   0)
G_table[5,Gset] <- c(1/4, 1/2, 1/4)
G_table[6,Gset] <- c(0  , 1/2, 1/2)
G_table[7,Gset] <- c(0  ,   1,   0)
G_table[8,Gset] <- c(0  , 1/2, 1/2)
G_table[9,Gset] <- c(0  ,   0,   1)

pGenoOffCond <- vector("list",3)
names(pGenoOffCond) <- Gset
for (i in Gset) pGenoOffCond[[i]] <- t(subset(G_table, Gf==i)[,Gset])

p.offgenotype <- function(G_, Gf, pGm) {
    ## Look up the probability of an offspring genotype, given the maternal genotype and
    ## proportions of offspring sired by each paternal genotype.

    ## Args:
    ##   G_: Offspring genotype.
    ##   Gf: Maternal genotype.
    ##   pGm: Relative fecundities of different male genotypes.

    ## Returns:
    ##   The probability of offspring genotype G_.

    ## expect Gf to have length=1
    if (length(Gf) != 1)
        stop("offspring genotype function not vectorised for maternal genotype")
    ## return the required probability
    return((pGenoOffCond[[Gf]] %*% pGm)[G_,1])
}
