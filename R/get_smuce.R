get_smuce = function(help, segs, alpha,ll,quant=0.35,rescale,constant,demography=F) { # removed: list.quantile.regs,
  gam = 0.25; eps = 0
  if(!demography) {pr1 = predict(LDJump::mod.full,help)} else {pr1 = predict(LDJump::mod.full.demo,help)}
  pr1[is.na(pr1)] = -1/gam;
  ind.q = which(seq(0.1, 0.5, by = 0.05) == quant)
  pr.cor = predict(LDJump::list.quantile.regs[[ind.q]], data.frame(x = pr1)); pr.cor = ifelse(pr.cor < -4, -4, pr.cor)
  pr.cor.nat = (pr.cor*gam+1)^(1/gam)-eps
  #### added lines 8-10 ####
  ind = as.numeric(which(rowSums(help) == 0))
  pr.cor.nat[ind] = NA
  pr.cor.nat = impute(pr.cor.nat, ind, two = T, segs = segs)
  if(constant) {return(list(pr.cor.nat))}
  else {
    idx = c(1,2,4:9,11,13);
    #  temp.cor.back=smuceR(pr.cor.nat, alpha = alpha, family = "gauss", confband = T)
    if(length(alpha) > 1) {
      temp.cor.back = list()
      for(i in 1:length(alpha)) {
        # temp.cor.back[[i]] = stepR::stepFit(pr.cor.nat, alpha = alpha[i], family = "gauss", confband = T)
        temp = stepR::stepFit(pr.cor.nat, alpha = alpha[i], family = "gauss", confband = T)
        if(rescale) {
          for(idxs in idx) {
            temp[[idxs]] = temp[[idxs]]/segs*ll
          }
        }
        temp.cor.back[[i]] = temp
      }
    } else {
      temp.cor.back = stepR::stepFit(pr.cor.nat, alpha = alpha, family = "gauss", confband = T)
      if(rescale) {
        for(idxs in idx) {
          temp.cor.back[[idxs]] = temp.cor.back[[idxs]]/segs*ll
        }  ## correction for 1-100% of sequence length
      }
    }
    return(list(temp.cor.back, pr.cor.nat, ind)) ### added ind
  }
}
