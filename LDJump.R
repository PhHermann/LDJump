LDJump = function(seqName = "", alpha = 0.05, segLength = 1000, pathLDhat = "", format = "fasta", refName = NULL, start = NULL, thth = 0.005, constant = F) {
  if(pathLDhat == "") {stop("Please provide the path of LDhat. Beware that this package requires LDhat to be installed for usage.")}
  if(seqName == "") {stop("Please provide the path of the sequence files in fasta/vcf format.")}
  if(format == "vcf") {
    vcfR_to_fasta(seqName, refName, ext.ind = T, cons = F, ext.haps = T, start = start)
  }
    s = readDNAStringSet(seqName); s = s[1:100]
    nn = length(s)
    if(constant) {
      segLength = nn
    }
    ll = floor(width(s)[1]/segLength)*segLength
    segs = ll/segLength
    system(paste("dos2unix ", paste(find.package("LDJump"),"/exec/Sums_LDHat_pack.sh", sep=""), sep = ""))
    help = t(sapply(1:segs,summary_statistics,s=s,segs=segs,thth=thth,seqName=seqName,nn=nn,ll = ll,pathLDhat = pathLDhat))
    # haps = read.table(paste("Confs_part_main.txt",sep=""),sep = ":")[,2]/ll*segs/nn
    sums = read.table("Sums_part_main.txt",sep = "=")
    indices = seq(1,nrow(sums),by=6)
    apwd = sums[indices+1,2]/ll*segs
    wath = sums[indices+2,2]/ll*segs
    vapw = sums[indices+5,2]/ll*segs
    help[,1] = help[,1]/ll*segs
    colnames(help) = c("fgts", "rsqu", "ldpr", "hahe")
    hats = read.table("resLDHats_pairwise_main.txt")[,5]/ll*segs
    # help = data.frame(cbind(help,apwd,wath,vapw,hats,haps), row.names = 1:nrow(help))
    help = data.frame(cbind(help,apwd,wath,vapw,hats), row.names = 1:nrow(help))
    full.list = get_smuce(help, segs, alpha,ll,constant=constant)
    if(!constant) {
      seq.full.cor = full.list[[1]]; pr.full.cor = full.list[[2]]
      return(list(seq.full.cor, pr.full.cor, help, alpha, nn, ll, segs))
    } else {
      pr.full.cor = full.list[[1]]
      return(list(pr.full.cor, help, alpha, nn, ll, segs))
    }
}
