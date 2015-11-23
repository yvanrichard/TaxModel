model{
  for(s in 1:NSPECIES){
    lmass[s] ~ dnorm(speciesmu[s],tau[s])
    # species mean drawn from family dist 
    speciesdev[s] <- speciesmu[s] - genusmu[genus[s]]
    speciesmu[s] <- genusmu[genus[s]] + genus.xi[genus[s]]*species.eta[s]
    
    #l_obs_sp[s] ~ dnorm(speciesmu[s]+species.xi[s]*sp_pred.eta[s],tau)
    # species tau drawn from hierarchical distr over all families
    #sp_pred.eta[s] ~ dnorm(0,species.prec)
    species.eta[s] ~ dnorm(0,genus.prec)
    species.xi[s] ~ dnorm(0,species.scale)
  }
  
  # genus fx from order distributions
  for(g in 1:NGENUS){
    # genus mean drawn from family dist
    genusdev[g] <- genusmu[g] - familymu[family[g]]
    genusmu[g] <- familymu[family[g]]+family.xi[family[g]]*genus.eta[g]
    #l_obs_gen[g] ~ dnorm(genusmu[g]+genus.xi[g]*gen_pred.eta[g],tau)
    # genus tau drawn from hierarchical distr over all families
    #gen_pred.eta[g] ~ dnorm(0,genus.prec)
    genus.eta[g] ~ dnorm(0,family.prec)
    genus.xi[g] ~ dnorm(0,genus.scale)
  }
  
  # family fx from order distribution
  for(f in 1:NFAMILIES){
    # family mean drawn from order dist
    familydev[f] <- familymu[f] - ordermu[order[f]]
    familymu[f] <- ordermu[order[f]]+order.xi[order[f]]*family.eta[f]
    #l_obs_fam[f] ~ dnorm(familymu[f]+family.xi[f]*fam_pred.eta[f],tau)
    # family tau drawn from hierarchical distr over all families
    #fam_pred.eta[f] ~ dnorm(0,family.prec)
    family.eta[f] ~ dnorm(0,order.prec)
    family.xi[f] ~ dnorm(0,family.scale)
  }
  
  # order fx from class distribution
  for(o in 1:NORDERS){
    # order mean drawn from order dist
    ordermu[o] <- grandmu+grand.xi*order.eta[o]
    #l_obs_ord[o] ~ dnorm(ordermu[o]+order.xi[o]*ord_pred.eta[o],tau)
    # order tau drawn from hierarchical distr over all orders
    #ord_pred.eta[o] ~ dnorm(0,order.prec)
    order.eta[o] ~ dnorm(0,grand.prec)
    order.xi[o] ~ dnorm(0,order.scale)
  }
  
  # priors tax hierachy
  species.scale ~ dgamma(0.001,0.001)
  species.prec ~ dgamma(0.5,0.5)
  
  family.scale  ~ dgamma(0.001,0.001)
  family.prec ~ dgamma(0.5,0.5)
  
  genus.scale  ~ dgamma(0.001,0.001)
  genus.prec ~ dgamma(0.5,0.5)
  
  order.scale   ~ dgamma(0.001,0.001)
  order.prec ~ dgamma(0.5,0.5)
  
  grand.xi ~ dnorm(0,0.0001)
  grand.prec ~ dgamma(0.5,0.5)
  #sigma.grand <- abs(grand.xi)/sqrt(grand.prec) 
  
  #grandtau ~ dgamma(0.01,0.01)
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
  sd.genus <- sd(genusdev)
  sd.species <- sd(speciesdev)
  sd.family  <- sd(familydev)
  sd.order   <- sd(ordermu)
  # sd.continent_eff <- sd(continent_eff)
  
  #tau ~ dgamma(0.01,0.01)
}