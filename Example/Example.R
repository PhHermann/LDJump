require(LDJump)
results = LDJump("/pathToSample/HatLandscapeN16Len1000000Nrhs15_th0.01_540_1.fa", alpha = 0.05, segLength = 1000, 
                 pathLDhat = "/pathToLDhat", format = "fasta", refName = NULL, thth = 0.01)
postscript("Results.eps", horiz = F)
plot(results[[1]], xlab = "Segments", ylab = "Estimated Recombination Rate", main = "Estimated recombination map with LDJump")
dev.off()
