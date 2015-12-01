require(dplyr)

sizes <- read.table("MOMv3.3.txt",sep='\t')
colnames(sizes) <- c("continent",'status','order','family','genus','species','lmass','cmass','ref')
sizes <- sizes %>% 
  filter(lmass>-100) %>% 
  unique() %>%
  group_by(order,family,genus,species) %>%
  summarise(lmass = mean(mean(lmass)))

# need numbers
n <- function(x) length(unique(x))
nidx <- function(x) as.numeric(factor(x,levels=unique(x)))


SPECIES = nidx(paste(sizes$genus,sizes$species,sep=' '))

## matching taxonomy for taxonomy indices
# match unique species to genus

GE <- nidx(sizes$genus)

# match unique species to family 
FAM <- nidx(sizes$family)

# match unique family to order
ORD <- nidx(sizes$order)

lpred <- 50
pred = sample.int(n = length(sizes$lmass),lpred)
jagsdata <- list(pred =pred ,
                 NPRED = lpred,
                 lmass = sizes$lmass,
                 CONTINENTS = n(sizes$continent),
                 NSPECIES = n(paste(sizes$genus,sizes$species,sep=' ')),
                 NGENUS = n(sizes$genus),
                 NFAMILIES = n(sizes$family),
                 NORDERS = n(sizes$order),
                 species = SPECIES,
                 genus = GE[match(unique(SPECIES),SPECIES)],
                 family = FAM[match(unique(GE),GE)],
                 order = ORD[match(unique(FAM),FAM)],
                 genus_pred = GE[pred],
                 family_pred = FAM[pred],
                 order_pred = ORD[pred],
                 tau=1/sizes$lmass*0.05
                 
)

require(R2jags)

JM <- jags.parallel(model.file = 'size-model.R',
                    n.iter = 12500,
                    n.burnin = 2500,
                    DIC = T,
                    n.thin = 5,
                    data=jagsdata,
                    n.chains = 3,
                    parameters.to.save = c('sd.species',
                                           'sd.genus',
                                           'sd.family',
                                           'sd.order',
                                           'mu_pred_g',
                                           'mu_pred_f',
                                           'mu_pred_o',
                                           'mu_pred_s',
                                           'species.scale',
                                           'genus.scale',
                                           'family.scale',
                                           'order.scale',
                                           'grandmu',
                                           'hc_scale'))



#JM
# traceplot(JM)

jagsdata_rfx_h <- list(lmass = sizes$lmass,
                       pred =pred ,
                       NPRED = lpred,
                       NSPECIES = n(paste(sizes$genus,sizes$species,sep=' ')),
                       NGENUS = n(sizes$genus),
                       NFAMILIES = n(sizes$family),
                       NORDERS = n(sizes$order),
                       species = SPECIES,
                       genus = GE,
                       family = FAM,
                       order = ORD,
                       genus_pred = GE[pred],
                       family_pred = FAM[pred],
                       order_pred = ORD[pred],
                       tau=1/sizes$lmass*0.05
                       
)

require(R2jags)

JM_rfx_h <- jags.parallel(model.file = 'size-model_rfx_h.R',
                          n.iter = 12500,
                          n.burnin = 2500,
                          DIC = T,
                          n.thin = 5,
                          data=jagsdata_rfx_h,
                          n.chains = 1,
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
                                                 'hc_scale'))



# JM_rfx_h
# traceplot(JM_rfx_h)

JM_rfx <- jags.parallel(model.file = 'size-model_rfx.R',
                        n.iter = 12500,
                        n.burnin = 2500,
                        DIC = T,
                        n.thin = 5,
                        data=jagsdata_rfx_h,
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
                                               'hc_scale'))



# JM_rfx
# traceplot(JM_rfx)
