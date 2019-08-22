vcf_statistics = function(seqName = "", alpha = 0.05, quant = 0.35, segLength = 1000, pathLDhat = "", pathPhi = "", format = NULL, refName = NULL, start = NULL, constant = F, rescale = F, status = T, polyThres = 0, cores = 1, accept = F, demography = F, regMod = "", out = "", lengthofseq=NA, chr = NA, startofseq = NA, endofseq = NA) {
  
  segs <- ceiling(lengthofseq/segLength)
  y <- 0
  ref = seqinr::read.fasta(file = refName)
  attr_name = attributes(ref)$names
  fa_start = 1
  fa_end = 0
  helper_new = c()
  ll = lengthofseq
  
  
  for (x in 1:segs){
    ix =  startofseq + (x-1) * segLength + y
    ex = ix + segLength - y
    
    y <- 1
    
    fa_start= (x-1) * segLength + y
    fa_end = fa_start + segLength -y
    
    system(paste0("vcftools --gzvcf ", seqName, " --chr ", as.integer(chr), " --from-bp ", as.integer(ix), " --to-bp ", as.integer(ex), "_sel.id --recode --recode-INFO-all --out ", "temp/sel_", as.integer(ix), "_", as.integer(ex)))
    vcfR_to_fasta(paste0("temp/sel_", as.integer(ix), "_", as.integer(ex), ".recode.vcf"), ref = ref, start = as.integer(ix), fa_start = fa_start, fa_end = fa_end, attr_name = attr_name)
    fastaName <- paste0("temp/sel_", as.integer(ix), "_", as.integer(ex), ".recode.vcf.fasta")
    if(x == 1){nn = Biostrings::readDNAStringSet(paste0("temp/sel_", as.integer(ix), "_", as.integer(ex), ".recode.vcf.fasta")); nn = length(nn)}
    }
    
    if(pathLDhat != "") {system(paste("dos2unix -q ", paste(find.package("LDJump"),"/exec/Sums_LDHat_pack.sh", sep=""), " --quiet", sep = ""))}
    if(cores == 1) {
      helper = t(sapply(1:segs,summary_statistics,segLength=segLength,segs = segs, seqName=seqName,nn=nn,pathLDhat = pathLDhat, pathPhi = pathPhi, status = status, polyThres = polyThres, out = out, format = format, startofseq = startofseq))
    } else {
       cl <- makeCluster(cores, type = "SOCK")
       helper = t(parSapply(cl, 1:segs,summary_statistics,segLength=segLength,segs = segs, seqName=seqName,nn=nn,pathLDhat = pathLDhat, pathPhi = pathPhi, status = status, polyThres = polyThres, out = out, format = format, startofseq = startofseq))
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
    helper_new = rbind(helper_new, helper)
  
    full.list = get_smuce(helper_new, segs, alpha, quant = quant, ll,rescale = rescale, constant=constant, demography = demography, regMod = regMod)
    if(!constant) {
      seq.full.cor = full.list[[1]]; pr.full.cor = full.list[[2]]; ind = full.list[[3]]
      if(length(ind) > 0) {warning(paste("Recombination rates were imputed for the following segments:", toString(ind), sep = ""))}
      return(list("Estimated recombination map" = seq.full.cor, "Constant estimates:" = pr.full.cor, "Summary Statistics" = helper_new, "alpha" = alpha, "Sample Size" = nn, "Sequence Length" = ll, "Segment Length" = segLength, "Imputed Segments" = ind))
    } else {
      pr.full.cor = full.list[[1]]
      return(list("Constant estimates" = pr.full.cor, "Summary Statistics" = helper_new, "alpha" = alpha, "Sample Size" = nn, "Sequence Length" = ll, "Segment Length" = segLength))
    }
  }
}
}
