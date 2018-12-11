getPhi = function(seqName, pathPhi, out, rm = T) {
  phiName = paste0("Phi",out,".log")
  system(paste0(pathPhi, " -f ", seqName, " -o -v >", phiName))
  MaxChi = as.numeric(strsplit(system(paste0("grep 'Value of maximum breakpoint is:' ",phiName), intern = T), ": ")[[1]][2])
  if(is.infinite(MaxChi)) {MaxChi = NA}
  # MaxChi = ifelse(is.infinite(MaxChi), 0, MaxChi)
  NSS = as.numeric(strsplit(system(paste("grep 'The Neighbour Similarity score is ' ",phiName), intern = T), " ")[[1]][6])
  # NSS = min(NSS, 1, na.rm = T)
  phi.mean = unlist(strsplit(system(paste0("grep 'Mean: ' ", phiName), intern = T), " "))
  phi.mean = as.numeric(phi.mean[phi.mean != ""][2])
  # phi.mean = max(as.numeric(phi.mean[phi.mean != ""][2]), 0, na.rm = T)
  phi.var = unlist(strsplit(system(paste0("grep 'Variance: ' ", phiName), intern = T), " "))
  phi.var = as.numeric(phi.var[phi.var != ""][2])
  # phi.var = max(as.numeric(phi.var[phi.var != ""][2]), 0, na.rm = T)
  if(rm) {system((paste0("rm ", phiName, " ", seqName)))}
  return(c(MaxChi, NSS, phi.mean, phi.var))
}
