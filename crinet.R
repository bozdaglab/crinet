# Load functions
source("crinet.functions.R")

# load sample data
load("CRINET_sample.rda")

# Specify the core no to use in parallel processing
core.no = 10

# Select hyperparameters
# singleErExp.thr = -0.01; comMir.thr = 1; hgeo.thr = 0.01; ceRNApartcor.thr = 0.55; ceRNApartcorPvalue.thr = 0.01; sumErExp.thr = -0.01
# Hyperparameters for sample data
singleErExp.thr = -0.01; comMir.thr = 1; hgeo.thr = 0.1; ceRNApartcor.thr = 0.1; ceRNApartcorPvalue.thr = 0.1; sumErExp.thr = 0.1

# Generate Effective Regulation
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
candidate.pairs = candidate.pairs[candidate.pairs$common.miRs >= comMir.thr,]

# Eliminate candidate pairs with respect to hypergeometric test based on common number of miRNAs
candidate.pairs$hypergeo = eliminate_wrt_hypergeo(pairs = candidate.pairs,
                                                  interactions = miRNA.target.interactions,
                                                  core.no = core.no)
candidate.pairs = candidate.pairs[candidate.pairs$hypergeo < hgeo.thr,]

# Eliminate candidate pairs with respect to significant correlation among their RNA expression excluding copy number aberration effect
candidate.pairs$ceRNA.part.cor = eliminate_wrt_partialCor(RNA = RNA,
                                                          cna = cna,
                                                          pairs = candidate.pairs,
                                                          core.no = core.no)
candidate.pairs = candidate.pairs[candidate.pairs$ceRNA.part.cor > ceRNApartcor.thr,]

candidate.pairs$ceRNA.part.corP = eliminate_wrt_partialCorP(RNA = RNA,
                                                            cna = cna,
                                                            pairs = candidate.pairs,
                                                            core.no = core.no)

candidate.pairs = candidate.pairs[candidate.pairs$ceRNA.part.corP < ceRNApartcorPvalue.thr,]

# Eliminate candidate pairs with respect to significant regulation collectively
candidate.pairs$sumErExp = eliminate_wrt_collectiveRegulation(RNA = RNA,
                                                              ERs = ERs,
                                                              pairs = candidate.pairs,
                                                              core.no = core.no)
candidate.pairs = candidate.pairs[candidate.pairs$sumErExp < sumErExp.thr,]

# 1.bootstrapping
candidate.pairs$ceRNA.part.cor.bootstrapP = eliminate_wrt_partialCor_bootstrapP(RNA = RNA,
                                                                                cna = cna,
                                                                                pairs = candidate.pairs,
                                                                                core.no = core.no)
candidate.pairs = candidate.pairs[candidate.pairs$ceRNA.part.cor.bootstrapP > ceRNApartcor.thr,]

# 2.bootstrapping
candidate.pairs$sumErExp.bootstrapP = eliminate_wrt_collectiveRegulation_bootstrapP(RNA = RNA,
                                                                                    ERs = ERs,
                                                                                    pairs = candidate.pairs,
                                                                                    core.no = core.no)

candidate.pairs = candidate.pairs[candidate.pairs$sumErExp.bootstrapP < sumErExp.thr,]

# Network Deconvolution
df = generate_deconv_input(pairs = candidate.pairs, scores = candidate.pairs$ceRNA.part.cor, needNormalization = T, deconv_in = paste0("Deconv_input.mat"))

# Run network deconvolution algorithm (ND.m) using Matlab [Input for Matlab is generated and named as: "Deconv_input.mat". This file will be used as input for Matlab]
# [Output of Matlab should be put into workingSpace named as "Deconv_output.mat"]

# Generated ceRNA interaction pairs
inferred.pairs = eliminate_wrt_networkDeconvolutionWithPercentage(pairs = candidate.pairs, deconv_out = paste0("Deconv_output.mat"), df = df, perc = .33)
