###################################################################################################
# STRUCTS
###################################################################################################

# Alignment
struct AlignTemp
    strand::String      # Methylation call strand
    R1::BAM.Record      # First record from left to right
    R2::BAM.Record      # Second record from left to right
end
mutable struct AllAlignTemps
    paired::Bool
    templates::Array{AlignTemp,1}
end

# CpelTdm Configuration
mutable struct CpeltdmConfig
    pe::Bool
    min_cov::Int64
    matched::Bool
    bound_check::Bool
    trim::NTuple{4,Int64}
    max_size_subreg::Int64
    max_size_anal_reg::Int64
    # Initializing method
    CpeltdmConfig(pe,min_cov,matched,bound_check,trim,max_size_subreg,max_size_anal_reg) = 
        new(pe,min_cov,matched,bound_check,trim,max_size_subreg,max_size_anal_reg)
end

# ROI data
mutable struct RoiData
    # Chromosome properties
    chr::String
    chrst::Int64
    chrend::Int64
    n_vec::Vector{Int64}
    cpg_pos::Vector{Int64}
    # Analysis flags
    diff_anal_done::Bool
    analyzed1::Vector{Bool}
    analyzed2::Vector{Bool}
    # Data
    mmls1::Vector{Float64}
    mmls2::Vector{Float64}
    nmes1::Vector{Float64}
    nmes2::Vector{Float64}
    θ1s::Vector{Vector{Float64}}
    θ2s::Vector{Vector{Float64}}
    mml_test::NTuple{2,Float64}
    nme_test::NTuple{2,Float64}
    cmd_test::NTuple{2,Float64}
    # Initializing method
    RoiData(s1,s2) = new("",0,0,[],[],false,fill(false,s1),fill(false,s2),fill(NaN,s1),fill(NaN,s2),
        fill(NaN,s1),fill(NaN,s2),fill([],s1),fill([],s2),(NaN,NaN),(NaN,NaN),(NaN,NaN))    
end