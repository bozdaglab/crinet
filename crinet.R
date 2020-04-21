# Initialization
{
  # Set workingspace
  workingSpace = getwd()
  saveDir = "supplementary"
}

# Load functions
source("crinet.functions.R")
# load data
load(paste(workingSpace, saveDir, paste("CRINET_data.rda", sep = ""), sep = "/"))
load(paste(workingSpace, saveDir, paste("CRINET_interactions.rda", sep = ""), sep = "/"))

core.no = 68
case.no = 1

# hyperparameters
singleErExp.thr = 0; comMir.thr = 1; hgeo.thr = 0.05; ceRNApartcor.thr = 0.3; ceRNApartcorPvalue.thr = 0.05; sumErExp.thr = -0.01

# generate Effective Regulation
ERs = generate_ERs(RNA = RNA,
                   miRNA = miRNA,
                   interactions = miRNA.target.interactions,
                   singleErExp.thr = singleErExp.thr,
                   core.no = core.no)

# Generate candidate pairs after eliminating wrt Effective Regulation
candidate.pairs = generate_candidate_pairs(ERs)
# Eliminate candidate pairs with respect to common number of miRNAs
candidate.pairs$common.miRs = eliminate_wrt_commonRegulator(pairs = candidate.pairs,
                                                            interactions = miRNA.target.interactions,
                                                            core.no = core.no)
candidate.pairs = candidate.pairs[candidate.pairs$common.miRs>=comMir.thr,]
# Eliminate candidate pairs with respect to hypergeometric test based on common number of miRNAs
candidate.pairs$hypergeo = eliminate_wrt_hypergeo(pairs = candidate.pairs,
                                                  interactions = miRNA.target.interactions,
                                                  core.no = core.no)
candidate.pairs = candidate.pairs[candidate.pairs$hypergeo<hgeo.thr,]
# Eliminate candidate pairs with respect to significant correlation among their RNA expression excluding copy number aberration effect
candidate.pairs$ceRNA.part.cor = eliminate_wrt_partialCor(RNA = RNA,
                                                          cna = cna,
                                                          pairs = candidate.pairs,
                                                          core.no = core.no)
candidate.pairs = candidate.pairs[candidate.pairs$ceRNA.part.cor>ceRNApartcor.thr,]
candidate.pairs$ceRNA.part.corP = eliminate_wrt_partialCorP(RNA = RNA,
                                                            cna = cna,
                                                            pairs = candidate.pairs,
                                                            core.no = core.no)
candidate.pairs = candidate.pairs[candidate.pairs$ceRNA.part.corP<ceRNApartcorPvalue.thr,]
# Eliminate candidate pairs with respect to significant regulation collectively
candidate.pairs$sumErExp = eliminate_wrt_collectiveRegulation(RNA = RNA,
                                                              ERs = ERs,
                                                              pairs = candidate.pairs,
                                                              core.no = core.no)
candidate.pairs = candidate.pairs[candidate.pairs$sumErExp<sumErExp.thr,]
# 1.bootstrapping
candidate.pairs$ceRNA.part.cor.bootstrapP = eliminate_wrt_partialCor_bootstrapPvalue(RNA = RNA,
                                                                                     cna = cna,
                                                                                     pairs = candidate.pairs,
                                                                                     thr = ceRNApartcor.thr,
                                                                                     core.no = core.no)
candidate.pairs = candidate.pairs[candidate.pairs$ceRNA.part.cor.bootstrapP<0.05,]
# 2.bootstrapping
candidate.pairs$sumErExp.bootstrapP = eliminate_wrt_collectiveRegulation_bootstrapPvalue(RNA = RNA,
                                                                                         ERs = ERs,
                                                                                         pairs = candidate.pairs,
                                                                                         thr = sumErExp.thr,
                                                                                         core.no = core.no)

candidate.pairs = candidate.pairs[candidate.pairs$sumErExp.bootstrapP<0.05,]
# Network Deconvolution
df = generate_deconv_input(pairs = candidate.pairs, scores = candidate.pairs$ceRNA.part.cor, needNormalization = T, deconv_in = paste0(def,"Deconv_input.mat"))

# Run network deconvolution using Matlab [generated input is named as: "Deconv_input.mat". Output should be put into workingSpace named as "Deconv_output.mat"]

inferred.pairs = eliminate_wrt_networkDeconvolution(pairs = candidate.pairs, deconv_out = paste0(def,"Deconv_output.mat"), df = df)
