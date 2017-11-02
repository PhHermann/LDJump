require(LDJump)
results = LDJump("/pathToSample/HatLandscapeN16Len1000000Nrhs15_th0.01_540_1.fa", alpha = 0.05, segLength = 1000, 
                 pathLDhat = "/pathToLDhat", format = "fasta", refName = NULL, thth = 0.01)
plot(results)
