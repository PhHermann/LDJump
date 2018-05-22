summary_statistics = function(x,s,segLength,segs,seqName,nn,thth,cor = 1,pathLDhat, status, polyThres) {
  # ix = round(1+(x-1)/segs*ll); ex = round(x/segs*ll)
  ix = 1+(x-1)*segLength; ex = x*segLength
  # sub = Biostrings::subseq(s, start = round(1+(x-1)/segs*ll), end=round(x/segs*ll))
  sub = Biostrings::subseq(s, start = 1+(x-1)*segLength, end=x*segLength)
  seqNamePart    = paste(Biostrings::substring(seqName,1,nchar(seqName)-3),"_part_",x,".fa",sep="")
  # seqNamePartHel = paste(Biostrings::substring(seqName,1,nchar(seqName)-3),"Hel_part.fa",sep="")
  Biostrings::writeXStringSet(sub, seqNamePart, format = "fasta")
  # system(paste("cp ", seqNamePart, " ", seqNamePartHel, sep = ""))
  if(length(ape::seg.sites(ape::read.FASTA(seqNamePart)))>1){
    system(paste(paste(find.package("LDJump"),"/exec/Sums_LDHat_pack.sh", sep=""), seqNamePart,nn,ex-ix+1,pathLDhat,formatC(thth, format = "fg"),cor,x,sep=" "))
    samp = ape::read.dna(seqNamePart, as.matrix = T, format = "fasta") # seqNamePartHel
    tajd = pegas::tajima.test(samp)$D
    temp = try(adegenet::DNAbin2genind(samp, polyThres = polyThres)); hahe = adegenet::Hs(temp);
    g2dftemp = adegenet::genind2df(temp)
    indices = t(combn(1:ncol(g2dftemp),2))
    helper = mapply(fgt_rrate_dpr, x = indices[,1], y = indices[,2], MoreArgs = list(data1 = g2dftemp, data2 = samp))
    if(ncol(helper) > 1) {stats = c(sum(helper[1,]),apply(cbind(sapply(1:(max(indices)-2),function(jj) rowMeans(helper[2:3,indices[,1]==jj])),helper[2:3,indices[,1]==max(indices-1)]),1,mean))}
    else {stats = helper}
  } else {
    stats = rep(0,3); hahe = 0; tajd = 0
    cat("Segregating sites=0\nAverage PWD=0\nWatterson theta=0\nTajima D statistic=0\nFu and Li D* statistic=0\nVariance PWD=0\n",file = paste0("Sums_part_main_", x, ".txt"),append=T)
    cat("Maximum at 4Ner(region) = 0.000 : Lk = NA\n", file = paste0("resLDHats_pairwise_main_", x, ".txt"),append=T)
  }
  if(status) {
    # nrow(read.table(file = "LDJump_Status.txt"))
    if(!file.exists("LDJump_Status.txt")) {
      cat(paste0(round(1/segs*100, 2), "% of Segments calculated\n"), file = "LDJump_Status.txt", append = T)
      print(paste0(round(1/segs*100, 2), "% of Segments calculated"))
    } else {
      cat(paste0(round((nrow(read.table("LDJump_Status.txt"))+1)/segs*100, 2), "% of Segments calculated\n"), file = "LDJump_Status.txt", append = T)
      print(paste0(round((nrow(read.table("LDJump_Status.txt"))+1)/segs*100, 2), "% of Segments calculated"))
    }

  }
  return(c(stats, hahe, tajd))
}
