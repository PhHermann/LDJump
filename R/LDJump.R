LDJump = function(seqName = "", alpha = 0.05, quant = 0.35, segLength = 1000, pathLDhat = "", pathPhi = "", format = "fasta", refName = NULL, start = NULL, constant = F, rescale = F, status = T, polyThres = 0, cores = 1, accept = F, demography = F, regMod = "", out = "", lengthofseq = NULL, chr = NULL, startofseq = NULL, endofseq = NULL) {
  if(pathPhi == "") {stop("Please provide the path of PhiPack. Beware that this package requires PathPhi to be installed for usage.")}
  if(seqName == "") {stop("Please provide the path of the sequence files in fasta/vcf format.")}
  if(format == "vcf") {
    
    dir.create("temp")  
    return(vcf_statistics(seqName, alpha, quant, segLength, pathLDhat, pathPhi, format, refName, start, constant, rescale, status, polyThres, cores, accept, demography, regMod, out, lengthofseq, chr, startofseq, endofseq))
    }
  if(format == "fasta") {
    if(length(ape::seg.sites(ape::read.dna(file = seqName, format = "fasta", as.matrix = T))) <= 1) {return("Data does not contain SNPs. Recombination rate cannot be estimated")}
    s = Biostrings::readDNAStringSet(seqName)
  }
  if(format == "DNABin") {
    s = DNAStringSet(get(seqName))
  }
  nn = length(s)
  if(constant) {
    ll = segLength = Biostrings::width(s)[1]
  }
  else {
    ll = floor(Biostrings::width(s)[1]/segLength)*segLength
  }
  segs = ll/segLength
  if(!check_continue(seqName = seqName, segs = segs, accept = accept, format = format)) {return()}
  if(pathLDhat != "") {system(paste("dos2unix -q ", paste(find.package("LDJump"),"/exec/Sums_LDHat_pack.sh", sep=""), " --quiet", sep = ""))}
  if(cores == 1) {
    helper = t(sapply(1:segs,summary_statistics,s=s,segLength=segLength,segs = segs, seqName=seqName,nn=nn,pathLDhat = pathLDhat, pathPhi = pathPhi, status = status, polyThres = polyThres, out = out))
  } else {
    cl <- makeCluster(cores, type = "SOCK")
    helper = t(parSapply(cl, 1:segs,summary_statistics,s=s,segLength=segLength,segs = segs, seqName=seqName,nn=nn,pathLDhat = pathLDhat, pathPhi = pathPhi, status = status, polyThres = polyThres, out = out))
    stopCluster(cl)
  }

  hahe = helper[,1]
  tajd = helper[,2]
  haps = helper[,3]/segLength/nn
  if(pathLDhat != "") {
    l.files = list.files(pattern=paste0("Sums_part_main_", out))
    l.files = l.files[order(as.numeric(regmatches(l.files, regexpr("[0-9]+", l.files))))]
    named.list <- lapply(l.files, read.table, sep = "=")
    sums = data.table::rbindlist(named.list)
    system(paste0("for f in Sums_part_main_", out, "*; do rm -f $f; done"))
    indices = seq(1,nrow(sums),by=6)
    apwd = sums[indices+1,2]/segLength
    vapw = sums[indices+5,2]/segLength
    colnames(apwd) = "apwd"; colnames(vapw) = "vapw"
  } else {
    apwd = helper[,4]/segLength
    vapw = helper[,5]/segLength
  }
  system(paste0("for f in *_part_",out,"*; do rm -f $f; done"))
  system(paste0("rm -f LDJump_Status", out, ".txt"))

  wath = helper[,6]/segLength
  MaxChi=helper[,7]
  NSS =  helper[,8]
  phi.mean = helper[,9]
  phi.var = helper[,10]
  helper = data.frame(cbind(hahe, tajd, haps, apwd, vapw, wath, MaxChi, NSS, phi.mean, phi.var), row.names = 1:nrow(helper))
  full.list = get_smuce(helper, segs, alpha, quant = quant, ll,rescale = rescale, constant=constant, demography = demography, regMod = regMod)
  if(!constant) {
    seq.full.cor = full.list[[1]]; pr.full.cor = full.list[[2]]; ind = full.list[[3]]
    if(length(ind) > 0) {warning(paste("Recombination rates were imputed for the following segments:", toString(ind), sep = ""))}
    return(list("Estimated recombination map" = seq.full.cor, "Constant estimates:" = pr.full.cor, "Summary Statistics" = helper, "alpha" = alpha, "Sample Size" = nn, "Sequence Length" = ll, "Segment Length" = segLength, "Imputed Segments" = ind))
  } else {
    pr.full.cor = full.list[[1]]
    return(list("Constant estimates" = pr.full.cor, "Summary Statistics" = helper, "alpha" = alpha, "Sample Size" = nn, "Sequence Length" = ll, "Segment Length" = segLength))
  }
}
