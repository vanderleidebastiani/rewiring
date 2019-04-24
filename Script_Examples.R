### Examples

### Requerid pacakges
require(bipartite)
require(plotrix) # for ploting

### Load the functions

source("./Functions/extinction.mod.R")
source("./Functions/one.second.extinct.mod.R")
source("./Functions/calc.mean.one.second.extinct.mod.R")
source("./Functions/matrix.p1.R")
source("./Functions/IC.R") # calc of 95 percent confidence interval

### Load the data

# Network - mutualistic plant-hummingbird networks
network <- as.matrix(read.csv("./DataSetExamples/SantaVirginia/SantaVirginia_dataset_network.csv", row.names = 1))
network

# Raw abundances of hummingbirds
h_abund <- as.matrix(read.csv("./DataSetExamples/SantaVirginia/SantaVirginia_dataset_h_abund.csv", row.names = 1))
h_abund

# Raw abundances of plants
pl_abund <- as.matrix(read.csv("./DataSetExamples/SantaVirginia/SantaVirginia_dataset_pl_abund.csv", row.names = 1))
pl_abund

# Trait (bill lenght) of hummingbirds
h_morph <- as.matrix(read.csv("./DataSetExamples/SantaVirginia/SantaVirginia_dataset_h_morph.csv", row.names = 1))
h_morph

# Trait (Cololla lenght) of plants
pl_morph <- as.matrix(read.csv("./DataSetExamples/SantaVirginia/SantaVirginia_dataset_pl_morph.csv", row.names = 1))
pl_morph

# Phenological distribution of hummingbirds
h_phen <- as.matrix(read.csv("./DataSetExamples/SantaVirginia/SantaVirginia_dataset_h_phen.csv", row.names = 1))
h_phen

# Phenological distribution of plants
pl_phen <- as.matrix(read.csv("./DataSetExamples/SantaVirginia/SantaVirginia_dataset_pl_phen.csv", row.names = 1))
pl_phen

### Define rewiring probabilities

## Abundances

# Relative abundances of pairwise bird and plant species
abundance <- (pl_abund/sum(pl_abund))%*%t(h_abund/sum(h_abund))
abundance

# Relative abundances of plants
abundance_pl <- sweep(abundance, 2, colSums(abundance), "/")
abundance_pl

# Relative abundances of hummingbirds
abundance_h <- sweep(abundance, 1, rowSums(abundance), "/") 
abundance_h

## Morphological matching

# Morphological matching
morphological <- 1-(as.matrix(vegdist(rbind(h_morph,pl_morph), method = "gower"))[(length(h_morph)+1):(length(h_morph)+length(pl_morph)),(1:length(h_morph))])
morphological

## Temporal coexistence

# Relative temporal coexistence (phenological overlap)
temporal <- pl_phen%*%t(h_phen)
temporal <- temporal/max(temporal)
temporal

## Fuzzy sets

# Distance between plant species
pl_morph_dist <- as.matrix(vegdist(pl_morph, method = "euclidean", na.rm = TRUE)) # Euclidean distance between species
pl_morph_dist

# Distance between hummingbirds species
h_morph_dist <- as.matrix(vegdist(h_morph, method = "euclidean", na.rm = TRUE)) # Euclidean distance between species
h_morph_dist

# Probability fuzzy trait similarity to plant
Tpl <- t(matrix.p1(t(network), pl_morph_dist)$matrix.P)
Tpl

# Probability fuzzy trait similarity to hummingbirds
Th <- matrix.p1(network, h_morph_dist)$matrix.P
Th

## All equal to one
one <- network
one[] <- 1
one

## Combined probabilities

# Morphological matching and phenological overlap
MP <- morphological*temporal # Hadamard product
MP

# Probability fuzzy trait similarity to hummingbirds and phenological overlap
ThP <- Th*temporal

# Probability fuzzy trait similarity to hummingbirds and phenological overlap
TplP <- Tpl*temporal

### Simule secondary extinctions in networks 

## Define the number of replications
nrep <- 1000

## Define the participant to be extinct
participant <- "lower"
participant

## Define method to extinct
method <- "random"
method

## Define method to rewiring
method.rewiring <- "one.try.single.partner"
method.rewiring

## Define the matrices of probability of rewiring

# To plants loss ("lower" participants)
probabilities.rewiring1 <- abundance_pl # relative abundances of plants
probabilities.rewiring1 # Probability of choice of a potential partner 

## Run secondary extinctions with specified parameters (results are list)
# The argument probabilities.rewiring2 to specify the probabilities of rewiring 
RES.without.rewiring <- replicate(nrep, one.second.extinct.mod(network, participant = participant, method = method, rewiring = FALSE), simplify = FALSE)
RES.with.rewiring.M <- replicate(nrep, one.second.extinct.mod(network, participant = participant, method = method, rewiring = TRUE, probabilities.rewiring1 = probabilities.rewiring1, probabilities.rewiring2 = morphological, method.rewiring = method.rewiring), simplify = FALSE)
RES.with.rewiring.P <- replicate(nrep, one.second.extinct.mod(network, participant = participant, method = method, rewiring = TRUE, probabilities.rewiring1 = probabilities.rewiring1, probabilities.rewiring2 = temporal, method.rewiring = method.rewiring), simplify = FALSE)
RES.with.rewiring.T <- replicate(nrep, one.second.extinct.mod(network, participant = participant, method = method, rewiring = TRUE, probabilities.rewiring1 = probabilities.rewiring1, probabilities.rewiring2 = Tpl, method.rewiring = method.rewiring), simplify = FALSE)
RES.with.rewiring.MP <- replicate(nrep, one.second.extinct.mod(network, participant = participant, method = method, rewiring = TRUE, probabilities.rewiring1 = probabilities.rewiring1, probabilities.rewiring2 = MP, method.rewiring = method.rewiring), simplify = FALSE)
RES.with.rewiring.TP <- replicate(nrep, one.second.extinct.mod(network, participant = participant, method = method, rewiring = TRUE, probabilities.rewiring1 = probabilities.rewiring1, probabilities.rewiring2 = TplP, method.rewiring = method.rewiring), simplify = FALSE)

## Calculates robustness to species extinctions
RES.robustness.without.rewiring <- sapply(RES.without.rewiring, robustness)
RES.robustness.with.rewiring.M <- sapply(RES.with.rewiring.M, robustness)
RES.robustness.with.rewiring.P <- sapply(RES.with.rewiring.P, robustness)
RES.robustness.with.rewiring.T <- sapply(RES.with.rewiring.T, robustness)
RES.robustness.with.rewiring.MP <- sapply(RES.with.rewiring.MP, robustness)
RES.robustness.with.rewiring.TP <- sapply(RES.with.rewiring.TP, robustness)

# Organize the results
res.robustness <- data.frame(robustness.without.rewiring = RES.robustness.without.rewiring, 
                             robustness.with.rewiring.M = RES.robustness.with.rewiring.M,
                             robustness.with.rewiring.P = RES.robustness.with.rewiring.P,
                             robustness.with.rewiring.T = RES.robustness.with.rewiring.T,
                             robustness.with.rewiring.MP = RES.robustness.with.rewiring.MP,
                             robustness.with.rewiring.TP = RES.robustness.with.rewiring.TP)
res.robustness

# Compute 95% confidence intervals
res.robustness.summary <- as.data.frame(t(sapply(res.robustness, IC)))
res.robustness.summary$rewiring <- c("0", "M", "P", "T", "MP", "TP") # add code for rewiring
res.robustness.summary

# Plot the results
plot(NA, xlim = c((1-0.5),(nrow(res.robustness.summary)+0.5)), ylim = c(0.8,1), 
     type = "n", las = 1, xaxt ="n", 
     ylab = "Robustness", xlab ="Mode of rewiring",
     main = paste("Participant =", participant, ";", "Method =", method))
axis(1, at = 1:nrow(res.robustness.summary), labels = res.robustness.summary$rewiring)
plotCI(1:nrow(res.robustness.summary), res.robustness.summary$mean, 
       ui = res.robustness.summary$upper, 
       li = res.robustness.summary$lower, 
       add = TRUE, pch = 19) # mean and confidence intervals
abline(h = c(res.robustness.summary$lower[1], res.robustness.summary$upper[1]), lty = 3) # lines of confidence interval to scenario without rewiring
