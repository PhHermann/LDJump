impute = function(data, index, two, segs) {
  while(!(length(index) == 0) || (length(index) == segs)) {
    data.to.impute = sapply(index, get_impute_data, data, two = two, segs = segs)
    help.ind = which(colSums(apply(data.to.impute, 2, is.na))==0)
    if(length(help.ind) > 0) {
      data[index[help.ind]] = apply(matrix(data.to.impute[,help.ind],nrow=2), 2, mean)
      index = as.numeric(which(is.na(data)))
      # return(impute(data, index, two = T, segs = segs))
    } else {
      data.to.impute = sapply(index, get_impute_data, data, two = F, segs = segs)
      ct = 1; help.ind=c()
      while(length(help.ind) == 0) {
        help.ind = which(colSums(apply(data.to.impute,2,is.na))==ct)
        ct = ct+1
      }
      if(length(help.ind) == 1) {imp.ind = help.ind} else {
        imp.ind = help.ind[colSums(apply(data.to.impute[1:2,help.ind],2,is.na))<2][1]
        if(is.na(imp.ind)) {imp.ind = help.ind[1]}
      }
      data[index[imp.ind]] = weighted.mean(data.to.impute[,imp.ind], c(1/3,1/3,1/6,1/6),na.rm = T)
      index = as.numeric(which(is.na(data)))
      # return(impute(data, index, two = T, segs = segs))
    }
  }
  return(data)
}
