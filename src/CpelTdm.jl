module CpelTdm

## Dependencies
using XAM                       # For BAM/SAM files
using FASTX                     # For FASTA/FASTQ files
using BED                       # For BED files
using GenomicFeatures           # Find overlaps between coordinates
using Calculus                  # For gradient computation
using Combinatorics             # For permutation test
using Dates                     # For printing with date and hour
using Distributed               # For pmap
using Distributions             # For simulations
using Optim                     # For estimating parameters
using Random                    # For simulations
using LinearAlgebra             # For matrix operations

## Includes
include("Structures.jl")
include("Bioinformatics.jl")
include("Inference.jl")
include("HypothesisTesting.jl")
include("Simulations.jl")

## Exports
export cpel_tdm

end # module
