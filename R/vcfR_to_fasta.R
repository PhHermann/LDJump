vcfR_to_fasta = function(seqName, refName = NULL, ext.ind = T, cons = F, ext.haps = T, start = NULL) {
  pop = vcfR::read.vcfR(file = paste(seqName, ".vcf", sep=""),
                   limit = 1e+07, nrows = -1, skip = 0, cols = NULL, convertNA = F, verbose = T)
  ref = ape::read.dna(file = paste(refName, ".fasta", sep = ""), format = "fasta")
  pop.dnabin = vcfR::vcfR2DNAbin(pop, extract.indels = ext.ind, consensus = cons,
                           extract.haps = ext.haps, ref.seq = ref,
                           start.pos = start, verbose = TRUE)
  ### with consensus = T we have 1 sequence per individual ####
  ape::write.dna(x = pop.dnabin, file = paste(seqName, ".fa", sep = ""), format = "fasta", colsep = "")
  print(paste(seqName, "converted to fasta file: ", seqName, ".fa", sep = ""))
}
