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
using DelimitedFiles            # For reading and writing delimited files
using MultipleTesting           # For BH correction

## CONSTANTS
const BG_BUFFER = 50000                     # Bytes of bedGraph records until write
const GFF_BUFFER = 500000                   # Bytes of GFF records until write
const THRESH_MAPQ = 39                      # MAPQ threshold (-10*log(p)) only true uni-reads
const FLAGS_ALLOWED = [0,16,83,99,147,163]  # Flags allowed in BAM recs
const ETA_MAX_ABS=5.0                       # Parameters space
const LOG2 = log(2)                         # Ln(2)

## Includes
include("Structures.jl")
include("IO.jl")
include("Bioinformatics.jl")
include("Inference.jl")
include("HypothesisTesting.jl")
include("Simulations.jl")

## Exports
export cpel_tdm

end # module
