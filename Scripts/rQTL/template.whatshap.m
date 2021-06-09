
SetDirectory["/scratch/aob2x/daphnia_hwe_sims/RABBIT/RABBIT_Packages/"]
Needs["MagicReconstruct`"]
Needs["MagicMap`"]

SetDirectory["/scratch/aob2x/daphnia_hwe_sims/Rabbit_phase_10cm/STEM"]
popScheme="STEM.ped"
ngroup = 1
model = "jointModel"
estfun = "origViterbiDecoding"
inputfile = "STEM.all.in"
resultFile = "STEM.all.out"
imputeTarget = "Offspring"

magicImpute[inputfile,model,popScheme,isFounderInbred -> False, outputFileID -> resultFile, isPrintTimeElapsed -> True, imputingTarget -> imputeTarget]


magicReconstruct[inputfile, model, popScheme, isFounderInbred -> False, outputFileID -> resultFile, reconstructAlgorithm -> estfun, isPrintTimeElapsed -> True]
summaryFile = StringDrop[resultFile, -4] <> ".csv"
saveAsSummaryMR[resultFile<>"_magicReconstruct.txt", summaryFile]

estfun = "origPosteriorDecoding"
resultFile = "STEM.all.out.post"
magicReconstruct[imputed, model, popScheme, isFounderInbred -> False, outputFileID -> resultFile, reconstructAlgorithm -> estfun, isPrintTimeElapsed -> True]
summaryFile = StringDrop[resultFile, -5] <> ".post.csv"
saveAsSummaryMR[resultFile<>"_magicReconstruct.txt", summaryFile]

Exit
