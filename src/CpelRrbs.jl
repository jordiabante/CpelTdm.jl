module CpelTdm

## Dependences
using XAM                       # For BAM/SAM files
using FASTX                     # For FASTA/FASTQ files
using BED                       # For BED files
using Calculus                  # For gradient computation
using Combinatorics             # For permutation test
using Dates                     # For printing with date and hour
using Distributed               # For pmap
using Distributions             # For simulations
using Optim                     # For estimating parameters
using Random                    # For simulations

## Includes
include("Structures.jl")
include("Bioinformatics.jl")
include("Inference.jl")
include("HypothesisTesting.jl")
include("Simulations.jl")

## Exports
export cpeltdm

end # module
