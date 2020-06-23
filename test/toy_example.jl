# Deps
using CpelTdm

# Set dir of test data
pack_dir = dirname(dirname(pathof(CpelTdm)))
data_dir = "$(pack_dir)/test/" 

# Set paths
fa = "$(data_dir)/fasta/ref.fa"
bed = "$(data_dir)/bed/roi.bed"
bams1 = "$(data_dir)/bam/g1_s" .* string.(collect(1:8)) .* ".bam"
bams2 = "$(data_dir)/bam/g2_s" .* string.(collect(1:8)) .* ".bam"

# Output
outprefix = "cpeltdm"
outdir = "$(data_dir)/cpeltdm/"

# Run analysis
rm(outdir;force=true,recursive=true)
cpel_tdm(bams1,bams2,bed,fa,outdir,outprefix;matched=false)
