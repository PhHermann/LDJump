LDJump = function(seqName = "", alpha = 0.05, segLength = 1000, pathLDhat = "", format = "fasta", refName = NULL, start = NULL, thth = 0.01, constant = F, rescale = F, status = T, polyThres = 0, cores = 1, accept = F, demography = F) {
  if(pathLDhat == "") {stop("Please provide the path of LDhat. Beware that this package requires LDhat to be installed for usage.")}
  if(seqName == "") {stop("Please provide the path of the sequence files in fasta/vcf format.")}
  if(format == "vcf") {
    vcfR_to_fasta(seqName, refName, ext.ind = T, cons = F, ext.haps = T, start = start)
  }
    s = Biostrings::readDNAStringSet(seqName);
    nn = length(s)
    if(constant) {
      ll = segLength = width(s)[1]
    }
    ll = floor(width(s)[1]/segLength)*segLength
    segs = ll/segLength
    if(!check_continue(seqName, segs = segs, accept)) {return()}
    system(paste("dos2unix -q ", paste(find.package("LDJump"),"/exec/Sums_LDHat_pack.sh", sep=""), " --quiet", sep = ""))
    # help = t(sapply(1:segs,summary_statistics,s=s,segLength=segLength,segs = segs, thth=thth,seqName=seqName,nn=nn,pathLDhat = pathLDhat, status = status, polyThres = polyThres))
    # haps = read.table(paste("Confs_part_main.txt",sep=""),sep = ":")[,2]/ll*segs/nn
    if(cores == 1) {
      help = t(sapply(1:segs,summary_statistics,s=s,segLength=segLength,segs = segs, thth=thth,seqName=seqName,nn=nn,pathLDhat = pathLDhat, status = status, polyThres = polyThres))
    } else {
      cl <- makeCluster(cores, type = "SOCK")
      help = t(parSapply(cl, 1:segs,summary_statistics,s=s,segLength=segLength,segs = segs, thth=thth,seqName=seqName,nn=nn,pathLDhat = pathLDhat, status = status, polyThres = polyThres))
      stopCluster(cl)
    }
    l.files = list.files(pattern="Sums_part_main_")
    l.files = l.files[order(as.numeric(regmatches(l.files, regexpr("[0-9]+", l.files))))]
    named.list <- lapply(l.files, read.table, sep = "=")
    sums = data.table::rbindlist(named.list)
    l.files = list.files(pattern="resLDHats_pairwise_main_")
    l.files = l.files[order(as.numeric(regmatches(l.files, regexpr("[0-9]+", l.files))))]
    named.list <- lapply(l.files, read.table)
    hats = data.table::rbindlist(named.list)
    hats = hats$V5/ll*segs
    system("rm -f resLDHats_pairwise* Sums_part_main* *_part_*.fa")
    indices = seq(1,nrow(sums),by=6)
    apwd = sums[indices+1,2]/ll*segs
    wath = sums[indices+2,2]/ll*segs
    vapw = sums[indices+5,2]/ll*segs
    help[,1] = help[,1]/ll*segs
    colnames(help) = c("fgts", "rsqu", "ldpr", "hahe")
    colnames(apwd) = "apwd"; colnames(wath) = "wath"; colnames(vapw) = "vapw"
    # hats = read.table("resLDHats_pairwise_main.txt")[,5]/ll*segs
    # help = data.frame(cbind(help,apwd,wath,vapw,hats,haps), row.names = 1:nrow(help))
    help = data.frame(cbind(help,apwd,wath,vapw,hats), row.names = 1:nrow(help))
    full.list = get_smuce(help, segs, alpha,ll,rescale = rescale, constant=constant, demography = demography)
    if(!constant) {
      seq.full.cor = full.list[[1]]; pr.full.cor = full.list[[2]]; ind = full.list[[3]]
      if(length(ind) > 0) {warning(paste("Recombination rates were imputed for the following segments:", toString(ind), sep = ""))}
      return(list("Estimated recombination map" = seq.full.cor, "Constant estimates:" = pr.full.cor, "Summary Statistics" = help, "alpha" = alpha, "Sample Size" = nn, "Sequence Length" = ll, "Segment Length" = segs, "Imputed Segments" = ind))
    } else {
      pr.full.cor = full.list[[1]]
      return(list("Constant estimates" = pr.full.cor, "Summary Statistics" = help, "alpha" = alpha, "Sample Size" = nn, "Sequence Length" = ll, "Segment Length" = segs))
    }
}
