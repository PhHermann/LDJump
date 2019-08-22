vcfR_to_fasta = function(seqName, refName = NULL, ext.ind = T, cons = F, ext.haps = T, start = NULL, , ref = NULL, fa_start= NULL, fa_end = NULL, attr_name = NULL) {

  pop = vcfR::read.vcfR(file = seqName, limit = 1e+07, nrows = -1, skip = 0, cols = NULL, convertNA = F, verbose = T)
  
  seqinr::write.fasta(sequences = ref[[attr_name]][fa_start:fa_end], names = names(ref), nbchar = 80, file.out = paste0("temp/", fa_start, "_", fa_end, "_out.fa"))
  
  pop.dnabin = vcfR::vcfR2DNAbin(pop, extract.indels = T, consensus = F,
                                 extract.haps = T, ref.seq = read.dna(paste0("temp/", fa_start, "_", fa_end, "_out.fa"), format = "fasta"),
                                 start.pos = start, verbose = TRUE)

  ape::write.dna(x = pop.dnabin, file = paste0(seqName, ".fasta", sep = ""), format = "fasta", colsep = "")
  print(paste(seqName, " converted to fasta file: ", seqName, ".fasta", sep = ""))
}
