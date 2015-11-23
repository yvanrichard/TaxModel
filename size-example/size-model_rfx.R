model{
  for(s in 1:NSPECIES){
    lmass[s] ~ dnorm(mu[s],tau[s])
    mu[s] <- grandmu + speciesmu[s] + genusmu[genus[s]] + familymu[family[s]] + ordermu[order[s]]
    
    speciesmu[s] <- genus.xi*species.eta[s]
    
      species.eta[s] ~ dnorm(0,genus.prec)
    
  }
  
  # genus fx from order distributions
  for(g in 1:NGENUS){
    # genus mean drawn from family dist
    genusmu[g] <- family.xi*genus.eta[g]
    # genus tau drawn from hierarchical distr over all families
    genus.eta[g] ~ dnorm(0,family.prec)
  }
  
  # family fx from order distribution
  for(f in 1:NFAMILIES){
    # family mean drawn from order dist
     familymu[f] <- order.xi*family.eta[f]
      family.eta[f] ~ dnorm(0,order.prec)
  }
  
  # order fx from class distribution
  for(o in 1:NORDERS){
    # order mean drawn from order dist
    ordermu[o] <- grand.xi*order.eta[o]
    order.eta[o] ~ dnorm(0,grand.prec)
  }
  
  # priors tax hierachy
  species.xi ~ dnorm(0,0.001)
  species.prec ~ dgamma(0.5,0.5)
  
  family.xi ~ dnorm(0,0.001)
  family.prec ~ dgamma(0.5,0.5)
  
  genus.xi ~ dnorm(0,0.001)
  genus.prec ~ dgamma(0.5,0.5)
  
  
  order.xi ~ dnorm(0,0.001)
  order.prec ~ dgamma(0.5,0.5)
  
  grand.xi ~ dnorm(0,0.0001)
  grand.prec ~ dgamma(0.5,0.5)
  #sigma.grand <- abs(grand.xi)/sqrt(grand.prec) 
  
  grandmu ~ dnorm(0,1e-6)
  
  
  #   # other rfx
  #   
  #   for (h in 1:CONTINENTS){
  #     continent_eff[h] <- cont.xi*cont.eta[h] 
  #     cont.eta[h] ~ dnorm(0,cont.prec)
  #   }
  #   
  #   cont.xi ~ dnorm(0,0.1)
  #   cont.prec ~ dgamma(0.5,0.5)
  
  # finite population sds
  sd.genus <- sd(genusmu)
  sd.species <- sd(speciesmu)
  sd.family  <- sd(familymu)
  sd.order   <- sd(ordermu)
  # sd.continent_eff <- sd(continent_eff)
  
  #tau ~ dgamma(0.001,0.001)
  
}