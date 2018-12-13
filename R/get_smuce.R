get_smuce = function(help, segs, alpha,ll,quant=0.35,rescale,constant,demography,regMod) { # removed: list.quantile.regs,
  gam = 0.5; eps = 0
  help$MaxChi[is.infinite(help$MaxChi)] = NA
  if(length(regMod) == 1 && !demography) {pr1 = predict(object = LDJump::mod.full,newdata = help)} ## mod.full
  if(length(regMod) == 1 && demography)  {pr1 = predict(object = LDJump::mod.full.demo,newdata = help)}
  if(length(regMod) > 1) {pr1 = predict(object = regMod, newdata = help)}
  pr1[is.na(pr1)] = -1/gam;
  ind.na = is.na(pr1)
  ind.q = which(round(seq(0.1, 0.5, by = 0.05),2 == quant))
  pr.cor = predict(object = LDJump::list.quantile.regs[[ind.q]], newdata = data.frame(x = pr1)); pr.cor = ifelse(pr.cor < -1/gam, -1/gam, pr.cor)
  pr.cor.nat = (pr.cor*gam+1)^(1/gam)-eps
  help = help[,1:8]
  #### added lines 8-10 ####
  # ind = as.numeric(which(rowSums(help) == 0))
  if(length(regMod) == 0) {
    ind = as.numeric(which(is.na(rowSums(help))))
  } else {
    ind = as.numeric(which(is.na(rowSums(help[,-c(ncol(help)-1, ncol(help))]))))
  }
  ind = c(ind, as.numeric(which(is.infinite(rowSums(help)))),as.numeric(which(rowSums(help) == 0)))
  ind = sort(unique(ind))
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
