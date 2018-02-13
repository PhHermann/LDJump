fgt_rrate_dpr = function(x,y,data1,data2) {
  temp = genetics::diseq(genetics::genotype(data1[,x],data1[,y],allow.partial.missing = T))
  t1 = temp$R2.overall
  t2 = temp$Dprime.overall
  temp = paste(data2[,x], data2[,y], sep = "")
  return(c(dim(table(temp))>=4,t1,t2))
}
