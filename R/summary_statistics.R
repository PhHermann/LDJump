summary_statistics = function(x,s,segLength,segs,seqName,nn,pathLDhat, pathPhi, status, polyThres, out) {
  ix = 1+(x-1)*segLength; ex = x*segLength
  sub = Biostrings::subseq(s, start = 1+(x-1)*segLength, end=x*segLength)
  seqNamePart    = paste(Biostrings::substring(seqName,1,nchar(seqName)-3),"_part_",x,out,".fa",sep="")
  Biostrings::writeXStringSet(sub, seqNamePart, format = "fasta")
  if(length(ape::seg.sites(ape::read.FASTA(seqNamePart)))>1){
    if(pathLDhat != "") {
      system(paste(paste0(find.package("LDJump"),"/exec/Sums_LDHat_pack.sh"),seqNamePart,nn,ex-ix+1,pathLDhat,x,out,sep=" "))
      apwd = vapw = NA
      }
    samp = ape::read.dna(seqNamePart, as.matrix = T, format = "fasta")
    tajd = pegas::tajima.test(samp)$D
    temp = try(adegenet::DNAbin2genind(samp, polyThres = polyThres)); hahe = adegenet::Hs(temp);
    if(pathLDhat != "") {system(paste0("sed '1d' ", seqNamePart, " > tmpfile",x,out,".fa"))}
    phis = getPhi(seqName = paste0("tmpfile",x,out,".fa"), pathPhi = pathPhi, out = paste0(x,out))
    haps = length(print(pegas::haplotype(samp)))
    if(pathLDhat == "") {
      s.Dist = Biostrings::stringDist(sub)
      apwd = mean(s.Dist); vapw = stats::var(s.Dist)/length(s.Dist)*(length(s.Dist)-1)
    }
    wath = pegas::theta.s(samp)
  } else {
    haps = hahe = tajd = apwd = vapw = wath = 0
    if(pathLDhat != "") {cat("Segregating sites=0\nAverage PWD=0\nWatterson theta=0\nTajima D statistic=0\nFu and Li D* statistic=0\nVariance PWD=0\n",file = paste0("Sums_part_main_",out,x,".txt"),append=T)}
    phis = rep(0, 4)
  }
  if(status) {
    # nrow(read.table(file = "LDJump_Status.txt"))
    if(!file.exists(paste0("LDJump_Status",out,".txt"))) {
      cat(paste0(round(1/segs*100, 2), "% of Segments calculated\n"), file = paste0("LDJump_Status",out,".txt"), append = T)
      print(paste0(round(1/segs*100, 2), "% of Segments calculated"))
    } else {
      cat(paste0(round((nrow(read.table(paste0("LDJump_Status",out,".txt")))+1)/segs*100, 2), "% of Segments calculated\n"), file = paste0("LDJump_Status",out,".txt"), append = T)
      print(paste0(round((nrow(read.table(paste0("LDJump_Status",out,".txt")))+1)/segs*100, 2), "% of Segments calculated"))
    }
  }
  system(paste0("rm ", seqNamePart))
  return(c(hahe, tajd, haps, apwd, vapw, wath, phis))
}
