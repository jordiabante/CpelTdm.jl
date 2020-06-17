![CI](https://github.com/jordiabante/CpelTdm.jl/workflows/CI/badge.svg)
![Docs](https://github.com/jordiabante/CpelTdm.jl/workflows/Docs/badge.svg)
[![MIT license](https://img.shields.io/badge/license-MIT-green.svg)](https://github.com/jordiabante/CpelTdm.jl/blob/master/LICENSE.md)

# Main commands

## Unmatched group comparison

In order to perform an unmatched group comparison, we first need
to load the package. Then, the right paths need to be specified,
including paths to the reference FASTA file, given by `fa` below,
the path to the BED file containing the regions of interest, given
by `bed` below, and the path to each bam file for each groups, given
by the `bams1` & `bams2` vectors below. In addition, the output 
directory `outdir` and `prefix` need to be defined as well.
Finally, we call function `cpel_tdm` setting the optional argument
`matched` equal to `false`, since in this case we wish to perform
an unmatched group comparison. In addition, here we assume that
we are using paired end sequencing data, hence we set `pe=true`.
However, if the sequencing data were single end, we would set
`pe=false` instead.

```julia

# Load dependencies
using Distributed
@everywhere using CpelTdm

# Set paths
dirname = "/path/to/data/"
fa = "$(dirname)/fasta/reference.fa"
bed = "$(dirname)/bed/regions_of_interest.bed"
bams1 = "$(dirname)/bam/".*["g1_s1.bam","g1_s2.bam","g1_s3.bam","g1_s4.bam"]
bams2 = "$(dirname)/bam/".*["g2_s2.bam","g2_s2.bam","g2_s3.bam","g2_s4.bam","g2_s5.bam"]

# Output paths
prefix = "cpeltdm"
outdir = "$(dirname)/cpeltdm-unmatched/"

# Run matched analysis
cpel_tdm(bams1,bams2,bed,fa,outdir,prefix;pe=true,matched=false)

```

If you wish to run CpelTdm on a cluster, we recommend you store
the previous code into a file, say `cpeltdm_unmatched_call.jl`. 
Then, in the submission script you can execute the content of 
this file using the command `julia -p 10 cpeltdm_unmatched_call.jl`, 
which will execute the code providing CpelTdm with 10 CPUs to 
perform the analysis.

## Matched group comparison

In order to perform an matched group comparison we need to follow
the same steps as in the unmatched group comparison as far as 
the paths go. In this case, however, CpelTdm assumes that the
i-th sample in `bams1` matches that of `bams2`. Finally, we call 
function `cpel_tdm` setting the optional argument `matched` equal 
to `true`, since in this case we wish to perform a matched group 
comparison. In addition, here we assume that we are using paired 
end sequencing data, hence we set `pe=true`. 

```julia

# Load dependencies
using Distributed
@everywhere using CpelTdm

# Set paths
dirname = "/path/to/data/"
fa = "$(dirname)/fasta/reference.fa"
bed = "$(dirname)/bed/regions_of_interest.bed"
bams1 = "$(dirname)/bam/".*["g1_s1.bam","g1_s2.bam","g1_s3.bam","g1_s4.bam","g1_s5.bam"]
bams2 = "$(dirname)/bam/".*["g2_s2.bam","g2_s2.bam","g2_s3.bam","g2_s4.bam","g2_s5.bam"]

# Output paths
prefix = "cpeltdm"
outdir = "$(dirname)/cpeltdm-matched/"

# Run matched analysis
cpel_tdm(bams1,bams2,bed,fa,outdir,prefix;pe=true,matched=true)

```

Similar to the unmatched case, if you wish to run CpelTdm on a 
cluster, we recommend you store the previous code into a file, 
say `cpeltdm_matched_call.jl`. Then, in the submission script you 
can execute the content of this file using the command 
`julia -p 10 cpeltdm_matched_call.jl`, which will execute the 
code providing CpelTdm with 10 CPUs to perform the analysis.
