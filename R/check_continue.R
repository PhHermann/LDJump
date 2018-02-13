check_continue = function(seqName, segs, accept) {
  fasta = ape::read.FASTA(file = seqName)
  lower.thres = seq(1,length(fasta[[1]])-1, by = length(fasta[[1]])/segs)
  upper.thres = c(lower.thres[-1], length(fasta[[1]]))
  thres = cbind(lower.thres, upper.thres)
  seg.test = ape::seg.sites(fasta)
  snps.per.seg = sapply(1:nrow(thres), function(i) {return(sum(seg.test < thres[i,2] & seg.test >= thres[i,1]))})
  if(sum(snps.per.seg <= 1)>1 & !accept) {
    input = readline(paste("There are ", sum(snps.per.seg <= 1), " segments with less than 2 SNPs. Do you want to continue and impute the estimated rate? We recommend to use larger segment lengths. Input here y/n: ", sep = ""))
    return(ifelse(input=="y", T, F))
  } else {return(T)}
}
