require(dplyr)

DD <- read.csv("DD.csv")

DD <- DD %>% filter(!is.na(EADDC))

# need numbers
n <- function(x) length(unique(x))
nidx <- function(x) as.numeric(factor(x,levels=unique(x)))


SPECIES = nidx(DD$Species)

## matching taxonomy for taxonomy indices
# match unique species to genus

GE <- nidx(DD$Genus)

# match unique species to family 
FAM <- nidx(DD$Family)

# match unique family to order
ORD <- nidx(DD$Order)


lpred <- 50
pred = sample.int(n = length(DD$EADDC),lpred)
DDdata_rfx_h <- list(lDD = log10(DD$EADDC),
                       pred =pred ,
                       NOBS = nrow(DD),
                       NPRED = lpred,
                       NSPECIES = n(DD$Species),
                       NGENUS = n(DD$Genus),
                       NFAMILIES = n(DD$Family),
                       NORDERS = n(DD$Order),
                       species = SPECIES,
                       genus = GE,
                       family = FAM,
                       order = ORD,
                     species_pred = SPECIES[pred],
                       genus_pred = GE[pred],
                       family_pred = FAM[pred],
                       order_pred = ORD[pred]
)

require(R2jags)

DD_rfx_h <- jags.parallel(model.file = 'DD-model_rfx_h.R',
                          n.iter = 30000,
                          n.burnin = 10000,
                          DIC = T,
                          n.thin = 10,
                          data=DDdata_rfx_h,
                          n.chains = 3,
                          parameters.to.save = c('sd.species',
                                                 'sd.genus',
                                                 'sd.family',
                                                 'sd.order',
                                                 'genus.scale',
                                                 'family.scale',
                                                 'order.scale',
                                                 'grandmu',
                                                 'grand.xi',
                                                 'mu_pred_g',
                                                 'mu_pred_f',
                                                 'mu_pred_o',
                                                 'mu_pred_s',
                                                 'sd.sigma.species',
                                                 'sd.sigma.genus',
                                                 'sd.sigma.family',
                                                 'sd.pop'))



#DD_rfx_h
#traceplot(DD_rfx_h)

DD_rfx <- jags.parallel(model.file = 'DD-model_rfx.R',
                        n.iter = 12500,
                        n.burnin = 2500,
                        DIC = T,
                        n.thin = 5,
                        data=DDdata_rfx_h,
                        n.chains = 3,
                        parameters.to.save = c('sd.species',
                                               'sd.genus',
                                               'sd.family',
                                               'sd.order',
                                               'genus.scale',
                                               'family.scale',
                                               'order.scale',
                                               'grandmu',
                                               'grand.xi',
                                               'mu_pred_g',
                                               'mu_pred_f',
                                               'mu_pred_o',
                                               'mu_pred_s',
                                               'sigma.species',
                                               'sigma.genus',
                                               'sigma.family',
                                               'sigma.order',
                                               'sd.pop'))



#DD_rfx
#traceplot(DD_rfx)
