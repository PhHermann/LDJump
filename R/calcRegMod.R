calcRegMod = function(n = c(10,16,20), len = c(500,1000,2000,3000,5000), thth = 0.01, nsim = 100, fr = c(), pathToScrm, scenario, pathToMs2dna, status = T, pathLDhat, pathPhi) {
data.all = c()
total = nsim*5*length(n)*length(len)
for(nn in n) {
 for(ll in len) {
  if(length(fr) == 0) {fr = round(c(runif(nsim, 0,0.01),runif(nsim,0.01,0.02),runif(nsim,0.02,0.05),runif(nsim,0.05,0.1),runif(nsim,0.1,0.2)),4)}
  help = c()
  for(j in 1:length(fr)){
   frfr = fr[j];  add = frfr*100000
   frS = frfr*(ll-1)
   thS = thth*(ll-1)
   seqName = paste0(getwd(),"/SimulatedPopulationN",nn,"Len",as.integer(ll),"th",thth,"rr",add,"nsim",j,".fa")
   system(paste0(pathToScrm,"/scrm ",nn," 1 -r ",frS," ",as.integer(ll)," -t ",thS,scenario,"-T | ", pathToMs2dna, "/ms2dna >", seqName))
   while(length(seg.sites(read.dna(seqName, format = "fasta"))) <= 1) {
     system(paste0(pathToScrm,"/scrm ",nn," 1 -r ",frS," ",as.integer(ll)," -t ",thS,scenario,"-T | ", pathToMs2dna, "/ms2dna >", seqName))
   } # end while; loop in order to avoid populations with too little snps
   help = rbind(help, c(unlist(LDJump(seqName, pathLDhat = pathLDhat, pathPhi = pathPhi, format = "fasta", constant = T, status = F, accept = T, regMod = "")[[2]])))
   if(status) {print(paste0(round((max(nrow(data.all)+j, j))/total, 4)*100, "% of in total ", total, " Simulations performed"))}
  } # end fr
 data.all = rbind(data.all, cbind(help, fr))
 }  # end ll
} # end n
# colnames(data.all)[10] = "tajd"
colNam = colnames(data.all)
data.all = data.frame(matrix(unlist(data.all), ncol = 10, byrow = F))
colnames(data.all) = colNam
data.all$fr = (data.all$fr^0.25-1)/0.25
regMod = gam(fr~haps+I(haps^2)+s(vapw)+s(apwd)+s(hahe)+s(wath)+s(MaxChi)+s(NSS)+s(phi.mean)+s(phi.var)+s(tajd),data=data.all,family = gaussian(),select = T)
return(list(regMod = regMod, data = data.all))
}
