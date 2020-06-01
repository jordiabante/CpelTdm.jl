###################################################################################################
# CONSTANTS
###################################################################################################
const BG_BUFFER = 50000                     # Bytes of bedGraph records until write
const GFF_BUFFER = 500000                   # Bytes of GFF records until write
const THRESH_MAPQ = 39                      # MAPQ threshold (-10*log(p)) only true uni-reads
const FLAGS_ALLOWED = [0,16,83,99,147,163]  # Flags allowed in BAM recs
###################################################################################################
# FUNCTIONS
###################################################################################################
"""
    `print_log(MESSAGE)`

    Function that prints MESSAGE to stderr.

    # Examples
    ```julia-repl
    julia> CpelTdm.print_log("Hello")
    [2020-03-30 16:24:18]: Hello
    ```
"""
function print_log(mess::String)

    # Right format for date
    date = Dates.format(now(), "yyyy-mm-dd HH:MM:SS")
    println(stderr,"[$(date)]: "  * mess)
    flush(stderr)

    # Return
    return nothing

end # end print_log
"""
    `get_align_strand(PAIRED_END,FLAG1,FLAG2)`

    Function that returns the strand of methylation call based on Bismark's logic. In single-end mode,
    OT and CTOT reads will both receive a FLAG of 0 while OB and CTOB both receive a FLAG of 16. In
    paired end mode:

                Read 1       Read 2
    OT:         99          147
    OB:         83          163
    CTOT:      147           99
    CTOB:      163           83

    This table was extracted from

        https://github.com/FelixKrueger/Bismark/issues/151

    # Examples
    ```julia-repl
    julia> CpelTdm.get_align_strand(true,UInt16(99),UInt16(147))
    "OT"
    ```
"""
function get_align_strand(pe::Bool,flag1::UInt16,flag2::UInt16)::String

    # Treat SE and PE separately
    if !pe
        if flag1==0
            s="OT"
        elseif flag1==16
            s="OB"
        else
            print_log("SE: was expecting FLAG 0 or and encountered $(flag1) instead.")
            print_log("Exiting julia ...")
            exit(1)
        end
    else
        if flag1==99 && flag2==147
            s="OT"
        elseif flag1==83 && flag2==163
            s="OB"
        elseif flag1==147 && flag2==99
            s="CTOT"
        elseif flag1==163 && flag2==83
            s="CTOB"
        else
            print_log("PE: unexpected flag combination. Expected 99/147, 147/99, 83/163, 163/83.")
            print_log("Exiting julia ...")
            exit(1)
        end
    end

    # Return strand
    return s

end # end get_align_strand
"""
    `order_bams(PAIRED_END,RECORDS)`

    Function that returns an AlignTemp object with R1 as the first record in the forward strand and
    with the methylation call strand taken from `get_align_strand`.

    # Examples
    ```julia-repl
    julia> CpelTdm.order_bams(true,RECORDS)
    ```
"""
function order_bams(pe::Bool,records::Vector{BAM.Record})::AlignTemp

    # Check which record appears first in forward strand
    if BAM.position(records[1])<=BAM.position(records[2])
        s = get_align_strand(pe, BAM.flag(records[1]), BAM.flag(records[2]))
        return AlignTemp(s,records[1],records[2])
    else
        s = get_align_strand(pe, BAM.flag(records[2]), BAM.flag(records[1]))
        return AlignTemp(s,records[2],records[1])
    end

end
"""
    `clean_records(PAIRED_END,RECORDS)`

    Function that takes a set of records and returns an AllAlignTemps object that contains all the
    properly aligned templates as an array of AlignTemp, which contains information about the
    methylation call strand as well as the relevant BAM records. In the PE case, R1 corresponds to
    the BAM record that appears before in the forward strand.

    # Examples
    ```julia-repl
    julia> CpelTdm.clean_records(true,RECORDS)
    ```
"""
function clean_records(pe::Bool,records::Vector{BAM.Record})::AllAlignTemps

    # Initialize struct
    out = AllAlignTemps(pe,[])

    # Consider PE vs. SE
    if pe
        # Handle PE case, first retrieve unique template names
        temp_names = unique([BAM.tempname(x) for x in records])
        for temp_name in temp_names
            # Get records with same template name
            temp_recs = filter(x->BAM.tempname(x)==temp_name,records)

            # There should only be two records with the same template name
            length(temp_recs)==2 && push!(out.templates,order_bams(pe,temp_recs))
        end
    else
        # Handle SE case
        out.templates = [AlignTemp(get_align_strand(pe,BAM.flag(x),UInt16(0)),x,BAM.Record()) for
            x in records]
    end

    # Return struct
    return out

end # end clean_records
"""
    `try_olaps(READER,CHR,WINDOW)`

    Function that tries to find BAM records overlaping with `CHR` at positions `WINDOW`.

    # Examples
    ```julia-repl
    julia> CpelTdm.try_olaps(reader,chr,win)
    ```
"""
function try_olaps(reader::BAM.Reader,chr::String,win::Vector{Int64})::Vector{BAM.Record}

    # NOTE: known bug that has to do w/ BioAlignments when close to end?
    records_olap =
    try
        collect(eachoverlap(reader, chr, win[1]:win[2]))
    catch
        Array{BAM.Record,1}()
    end

    return records_olap

end # try_olaps
"""
    `read_bam(BAM_PATH,CHR,FEAT_ST,FEAT_END,CPG_POS,CHR_SIZE,PE,TRIM)`

    Function that reads in BAM file in `BAM_PATH` and returns methylation vectors for those records
    that overlap with (1-based) genomic coordinates `chr:chrSt-chrEnd` at `cpg_pos`. A trimming given
    by TRIM=(Rf_5p,Rf_3p,Rr_5p,Rr_3p) is applied to the reads. The information was taken from:

        XM: meth calls (https://github.com/FelixKrueger/Bismark/tree/master/Docs)
        XS: uniqueness (http://bowtie-bio.sourceforge.net/bowtie2/manual.shtml)

    For info on OT, CTOT, OB, CTOB nomenclature see

        http://software.broadinstitute.org/software/igv/book/export/html/37.

    # Examples
    ```julia-repl
    julia> CpelTdm.read_bam(BAM_PATH,"chr1",30,80,[40,60],1000,false,(0,0,0,0))
    ```
"""
function read_bam(bam::String,chr::String,roi_st::Int64,roi_end::Int64,cpg_pos::Vector{Int64},
                  chr_size::Int64,pe::Bool,trim::NTuple{4,Int64})::Array{Vector{Int64},1}

    # Number of CpG sites is determined by that in the region
    N = length(cpg_pos)

    # Expand window in PE case to include pair even if it's outside the window
    exp_win = pe ? [max(1,roi_st-75),min(chr_size,roi_end+75)] : [roi_st,roi_end]

    # Get records overlapping window.
    reader = open(BAM.Reader,bam,index=bam*".bai")
    records_olap = try_olaps(reader,chr,exp_win)
    length(records_olap)>0 || return Vector{Int64}[]
    close(reader)

    # Relevant flags in BAM file (both Bismark and Arioc)
    filter!(x-> (BAM.ismapped(x)) && (haskey(x,"XM")) && (!haskey(x,"XS")) && (BAM.flag(x) in
           FLAGS_ALLOWED) && (BAM.mappingquality(x)>THRESH_MAPQ),records_olap)

    # Clean records & organize them
    records_org = clean_records(pe,records_olap)

    # Loop over organized records
    xobs = Vector{Int64}[]
    for record in records_org.templates
        # Get methylation call and offset (depends on strand where call is made)
        if !(records_org.paired)
            # Obtain methylation call from single end ()
            meth_call = record.R1[:"XM"]
            meth_call = SubString(meth_call,(1+trim[1]):(length(meth_call)-trim[2]))
            OFFSET = record.strand in ["OT","CTOT"] ? 1 : 2
            OFFSET -= BAM.position(record.R1)+trim[1]
        else
            # Obtain methylation call
            R1_call = record.R1[:"XM"]
            R2_call = record.R2[:"XM"]
            R1_call = SubString(R1_call,(1+trim[1]):(length(R1_call)-trim[2]))
            R2_call = SubString(R2_call,(1+trim[4]):(length(R2_call)-trim[3]))
            dist_st = abs(BAM.position(record.R2)+trim[4] - BAM.position(record.R1)-trim[1])
            meth_call = R1_call * "."^max(0,dist_st-length(R1_call))
            meth_call *= R2_call[max(1,(length(R1_call)+1-dist_st)):end]
            OFFSET = record.strand in ["OT","CTOT"] ? 1 : 2
            OFFSET -= BAM.position(record.R1)+trim[1]
        end

        # Cross positions of CpG sites if template contains CpG sites
        obs_cpgs = [x.offset for x in eachmatch(r"[zZ]",meth_call)] .- OFFSET
        length(obs_cpgs)>0 || continue

        # If overlapping CpG sites then store and add to xobs
        olap_cpgs = findall(x-> x in obs_cpgs, cpg_pos)
        length(olap_cpgs)>0 || continue
        x = zeros(Int64,N)
        obs_cpgs = obs_cpgs[findall(x-> x in cpg_pos,obs_cpgs)] .+ OFFSET
        x[olap_cpgs] = reduce(replace,["Z"=>1,"z"=>-1],init=split(meth_call[obs_cpgs],""))

        # Add to set of observations
        push!(xobs,x)

    end # end loop over templates sequenced

    # Return
    return xobs

end # end read_bam
"""
    `get_paths(BAMS1,BAMS2,FA,BED,OUTDIR,OUTPREFIX)`

    Function that returns output file names.

    # Examples
    ```julia-repl
    julia> CpelTdm.get_paths(bams1,bams2,fa,bed,outdir,outprefix)
    ```
"""
function get_paths(bams1::Vector{String},bams2::Vector{String},fa::String,bed::String,outdir::String,
                   outprefix::String)::NTuple{3,Vector{String}}

    # Check index files exist
    ind_miss = !isfile(fa*".fai")
    for bam in vcat(bams1,bams2)
        if !isfile(bam*".bai")
            ind_miss |= true
            print_log("Index files missing for")
            print_log("$(bam)")
        end
    end
    ind_miss && return ([],[],[])

    # Individual bedGraph output files
    prefixes = [String(split(basename(bam),".")[1]) for bam in bams1]
    append!(prefixes,[String(split(basename(bam),".")[1]) for bam in bams2])
    out_mml_paths = ["$(outdir)/$(prefix)_mml.bedGraph" for prefix in prefixes]
    out_nme_paths = ["$(outdir)/$(prefix)_nme.bedGraph" for prefix in prefixes]
    
    # Check for existance of at least an output files
    if all(isfile.(vcat(out_mml_paths,out_nme_paths)))
        print_log("At least an MML or NME file already exists...")
        return ([],[],[])
    end

    # Differentiial analysis bedGraph output files
    out_diff_paths = "$(outdir)/$(outprefix)_" * ["Tmml","Tnme","Tcmd"] .* "_diff_analysis.bedGraph"
    
    # Check for existance of at least an output files
    if all(isfile.(out_diff_paths))
        print_log("At least a differential output file already exists...")
        return ([],[],[])
    end

    # Check BED file exists
    if !(isfile(bed))
        print_log("BED file does not exist...")
        return ([],[],[])
    end

    # Return paths
    return out_mml_paths,out_nme_paths,out_diff_paths
    
end
"""
    `split_bed_record(bed_rec,max_size_anal_reg)`

    Function that splits regions of interest in BED into analysis regions.

    # Examples
    ```julia-repl
    julia> str_j_rec = "chr1\t1000\t2001\t."
    julia> using BED; bed_rec = BED.Record(str_j_rec);
    julia> CpelTdm.split_bed_record(bed_rec,500)
    2-element Array{BED.Record,1}:
    BED.Record:
        chromosome: chr1
                start: 1000
                end: 1500
                name: .   
    BED.Record:
        chromosome: chr1
                start: 1500
                end: 2001
                name: .
    ```
"""
function split_bed_record(bed_rec::BED.Record,max_size_anal_reg::Int64)::Vector{BED.Record}

    ## BED files are:
    ## - 0-based coordinates
    ## - chromstart is inclusive
    ## - chromend is non-inclusive
    ## When read in by BED.jl the start becomes one-based.

    # Get size of record
    size_rec = BED.chromend(bed_rec) - BED.chromstart(bed_rec)

    # Initialize
    recs = Vector{BED.Record}()

    # Divide into analysis regions
    chr = BED.chrom(bed_rec)
    n_anal_reg = Int(ceil(size_rec/max_size_anal_reg))
    size_anal_reg = Int(floor(size_rec/n_anal_reg))

    # Create n_anal_reg analysis regions
    @inbounds for j=1:n_anal_reg    
        
        # Get coordinates
        chr_st = BED.chromstart(bed_rec) + (j-1) * size_anal_reg - 2
        chr_end = j==n_anal_reg ? BED.chromend(bed_rec) : chr_st + size_anal_reg + 1

        # Create string with BED record
        str_j_rec = "$(chr)\t$(chr_st)\t$(chr_end)\t."

        # Push record
        push!(recs,BED.Record(str_j_rec))

    end
    
    # Return vector
    return recs

end
"""
    `get_bed_out_filename(bed)`

    Function that returns name of output BED file.

    # Examples
    ```julia-repl
    julia> CpelTdm.get_bed_out_filename(bed)
    ```
"""
function get_bed_out_filename(bed::String)::String

    # Return name
    return dirname(bed) * "/" * split(basename(bed),".")[1] * "_cpel.bed" 
    
end
"""
    `gen_anal_reg_file(bed,max_size_anal_reg)`

    Function that dvides regions of interest in BED into analysis regions.

    # Examples
    ```julia-repl
    julia> max_size_anal_reg = 500;
    julia> bed_out = CpelTdm.get_bed_out_filename(bed);
    julia> CpelTdm.gen_anal_reg_file(bed,bed_out,max_size_anal_reg)
    ```
"""
function gen_anal_reg_file(bed::String,bed_out::String,max_size_anal_reg::Int64)::Nothing

    # Check if BED file already exists
    isfile(bed_out) && (print_log("Found CpelTdm BED file ..."); return nothing)

    # Read new BED file
    reader = open(BED.Reader,bed)
    record = BED.Record()
    anal_regs = Vector{BED.Record}()
    while !eof(reader)
        # Get record
        read!(reader, record)
        # Append new analysis regions
        append!(anal_regs,split_bed_record(record,max_size_anal_reg))
    end
    close(reader)

    # Write new BED file
    output = open(bed_out, "a")
    writer = BED.Writer(output)
    while length(anal_regs)>0
        write(writer,popfirst!(anal_regs))
    end
    close(writer)

    # Return nothing
    return nothing
    
end
"""
    `get_anal_reg_chr(bed,chr)`

    Function that returns analyses regions in chr.

    # Examples
    ```julia-repl
    julia> CpelTdm.get_anal_reg_chr(bed,"chr1")
    ```
"""
function get_anal_reg_chr(bed::String,chr::String)::Vector{BED.Record}

    # Read new BED file
    record = BED.Record()
    reader = open(BED.Reader,bed)
    anal_regs = Vector{BED.Record}()
    while !eof(reader)
        # Get record
        read!(reader,record)
        # Keep if in chromosome
        BED.chrom(record)==chr && CpelTdm.print_log(BED.chrom(record))
        BED.chrom(record)==chr && push!(anal_regs,record)
    end
    close(reader)

    # Return nothing
    return anal_regs
    
end
"""
    `write_bedGraph(CHR,STATS,PATH)`

    Function that writes bedGraph files.

    # Examples
    ```julia-repl
    julia> CpelTdm.write_bedGraph(chr,stats,path)
    ```
"""
function write_bedGraph(chr::String,stats::Vector{Tuple{Int64,Int64,Float64}},path::String)::Nothing
    
    # Write
    open(path, "a") do f
        for i=1:length(stats)
            write(f,"$(chr)\t"*join(string.(collect(stats[i])),"\t"),"\n")
        end
    end

    # Return nothing
    return nothing

end
"""
    `write_output(CHR,OUTPMAP,MML_PATHS,NME_PATHS,DIFF_PATHS)`

    Function that write output of pmap into respective files.

    # Examples
    ```julia-repl
    julia> CpelTdm.write_output(chr,outpmap,mml_paths,nme_paths,diff_paths)
    ```
"""
function write_output(chr::String,outpmap::Vector{Vector{Float64}},mml_paths::Vector{String},
                      nme_paths::Vector{String},diff_paths::Vector{String})::Nothing
    
    # Write MMLs
    for i=1:length(mml_paths)
        write_bedGraph([x[i] for x in outpmap],chr,mml_paths[i])
    end

    # Write NMEs
    offset = length(mml_paths)
    for i=1:length(nme_paths)
        write_bedGraph([x[offset+i] for x in outpmap],chr,nme_paths[i])
    end

    # Write differential analysis
    offset += length(nme_paths)
    write_bedGraph([x[offset+1] for x in outpmap],chr,diff_paths[1])
    write_bedGraph([x[offset+2] for x in outpmap],chr,diff_paths[2])
    write_bedGraph([x[offset+3] for x in outpmap],chr,diff_paths[3])

    # Return nothing
    return nothing

end
"""
    `cpel_tdm(BAMS1,BAMS2,BED,FA,OUTDIR)`

    Function to call to CpelTdm analysis. 
        - BAMS1: BAM files associated to group 1.
        - BAMS2: BAM files associated to group 2.
        - BED: BED file containing regions of interest.
        - FA: reference genome.
        - OUTDIR: output directory to store results.

    # Examples
    ```julia-repl
    julia> CpelTdm.cpel_tdm(bams1,bams2,bed,fa,outdir)
    ```
"""
function cpel_tdm(bams1::Vector{String},bams2::Vector{String},bed::String,fa::String,outdir::String,prefix::String;
                   pe::Bool=true,max_size_anal_reg::Int64=1000,max_size_subreg::Int64=250,cov_ths::Int64=5,
                   matched::Bool=false,trim::NTuple{4,Int64}=(5,0,5,0),bound_check::Bool=false)::Nothing

    # Print initialization of juliASM
    print_log("Starting CpelTdm ...")

    ## IO
    print_log("Checking IO ...")
    
    # Get paths
    mml_paths,nme_paths,diff_paths = get_paths(bams1,bams2,fa,bed,outdir,prefix)

    # Exist if any output file exists
    if length.(mml_paths).==0
        print_log("Output files already exist. Delete them to run the analysis again.")
        exit(0)
    end

    # Create output folder if it doesn't exist
    isdir(outdir) || mkdir(outdir)

    ## Define analysis regions
    print_log("Defining analysis regions from BED file ...")

    # Divide BED file divided into analysis regions
    bed_out = get_bed_out_filename(bed)
    gen_anal_reg_file(bed,bed_out,max_size_anal_reg)

    ## Inference

    # Estimate θ, compute MML & NME, and do test for each ROI
    print_log("Running analysis ...")
    anal_bed_file(bams1,bams2,bed_out,fa,mml_paths,nme_paths,diff_paths;pe=pe,max_size_anal_reg=max_size_anal_reg,
                  max_size_subreg=max_size_subreg,cov_ths=cov_ths,matched=matched,trim=trim,bound_check=bound_check)

    # Print done
    print_log("Done.")

    # Return
    return nothing

end
"""
    `anal_bed_file(BAMS1,BAMS2,BED,FA,MML_PATHS,NME_PATHS,DIFF_PATHS)`

    Function that estimates θ for all samples, computes statistics, and does a matched/unmatched 
    permutation test at each ROI.

    # Examples
    ```julia-repl
    julia> CpelTdm.anal_bed_file(bams1,bams2,bed,fa,mml_paths,nme_paths,diff_paths)
    ```
"""
function anal_bed_file(bams1::Vector{String},bams2::Vector{String},bed::String,fa::String,mml_paths::String,
                       nme_paths::String,diff_paths::String;pe::Bool=true,max_size_anal_reg::Int64=1000,
                       max_size_subreg::Int64=250,cov_ths::Int64=5,matched::Bool=false,
                       trim::NTuple{4,Int64}=(5,0,5,0),bound_check::Bool=false)::Nothing

    # Find chromosomes
    reader_fa = open(FASTA.Reader,fa,index=fa*".fai")
    chr_sizes = reader_fa.index.lengths
    chr_names = reader_fa.index.names
    close(reader_fa)

    # Loop over chromosomes
    for chr in chr_names

        # Get windows pertaining to current chromosome
        print_log("Processing chr: $(chr) ...")
        reg_chr = get_anal_reg_chr(bed,chr)
        chr_size = chr_sizes[findfirst(x->x==chr,chr_names)]

        # Process regions in chromosome
        # out_pmap = pmap(x->pmap_dmr_finder(x,chr,chr_size,bam,gff,fa,pe,g_max,cov_ths,trim,bound_check),regions_chr)
        # length(out_pmap)>0 || continue

        # Add last to respective bedGraph file
        # write_output(chr,out_pmap,mml_paths,nme_paths,diff_paths)

    end

    # Sort bedGraph if there's output
    # if all(isfile.(out_paths))
    #     sort_bedgraphs(out_paths)
    # else
    #     print_log("No output was created. Check the chromosome names in FASTA and GFF match ...")
    # end

    # Return nothing
    return nothing

end # end dmr_finder
# """
#     `pmap_dmr_finder(FEAT,CHR,CHR_SIZE,BAM,GFF,FA,PE,G,COV_THS,TRIM)`

#     Function that estimates θ and computes MML/NME for all samples, and does permutation test for a given ROI.

#     # Examples
#     ```julia-repl
#     julia> pmap_dmr_finder(FEAT,CHR,CHR_SIZE,BAM,GFF,FA,PE,G,COV_THS,TRIM)
#     ```
# """
# function pmap_dmr_finder(feat::GFF3.Record,chr::String,chr_size::Int64,bam::String,gff::String,
#                          fa::String,pe::Bool,g_max::Int64,cov_ths::Int64,trim::NTuple{4,Int64},
#                          bound_check::Bool)#::Vector{Tuple{Int64,Int64,Float64,Int64,Int64}}

#     # Empty output
#     # nan_out = [(0,0,0.0,0,0),(0,0,0.0,0,0)]

#     # Get window of interest
#     # f_st = GFF3.seqstart(feat)
#     # f_end = GFF3.seqend(feat)

#     # Get CpG sites
#     # cpg_pos = get_cpg_pos(Dict(GFF3.attributes(feat)))
#     # length(cpg_pos[1])>0 || return nan_out

#     # Get vector of Ns
#     # nvec = get_ns(cpg_pos[1],g_max,f_st,f_end)

#     # Get vectors from BAM overlapping region
#     # xobs = read_bam(bam,chr,f_st,f_end,cpg_pos[1],chr_size,pe,trim)
#     # obs_per_cpg = get_obs_per_cpg(xobs)
#     # (cov_ths<=mean_cov(xobs)<=400) && (sum(obs_per_cpg.==0)<=1.0/5.0*sum(nvec)) || return nan_out

#     # Estimate parameters of CPEL models
#     # θhat = est_theta_sa(nvec,xobs)
#     # bound_check && check_boundary(θhat) && return nan_out

#     # Estimate intermediate quantities
#     # ∇1 = get_grad_logZ(nvec,θhat)

#     # Compute output
#     # mml = comp_mml_∇(nvec,∇1)
#     # nme = comp_nme_∇(nvec,θhat,∇1)

#     # Return output
#     return true

# end # end pmap_allele_agnostic_chr