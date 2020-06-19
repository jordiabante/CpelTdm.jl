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
function print_log(mess::String)::Nothing

    # Right format for date
    date = Dates.format(now(), "yyyy-mm-dd HH:MM:SS")
    println(stderr,"[$(date)]: "  * mess)
    flush(stderr)

    # Return
    return nothing

end
"""
    `sleep_exit(MESSAGE)`

    Function that prints MESSAGE, waits 5 seconds and exits without error.

    # Examples
    ```julia-repl
    julia> CpelTdm.sleep_exit("Bye")
    [2020-03-30 16:24:18]: Bye
    ```
"""
function sleep_exit(mess::String)::Nothing

    # Print message
    print_log(mess)
    sleep(10)
    exit(0)

    # Return
    return nothing

end
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
    `get_paths(BAMS1,BAMS2,FASTA,BED,OUTDIR,OUTPREFIX)`

    Function that returns output file names.

    # Examples
    ```julia-repl
    julia> out_mml_paths,out_nme_paths,out_diff_paths = CpelTdm.get_paths(bams1,bams2,fasta,bed,outdir,outprefix)
    ```
"""
function get_paths(bams1::Vector{String},bams2::Vector{String},fasta::String,bed::String,outdir::String,
                   outprefix::String)::NTuple{3,Vector{String}}

    # Check BED file exists
    if !isfile(bed)
        sleep_exit("BED file was not found ...")
    end

    # Check FASTA index file exists
    if !all([isfile(fasta),isfile(fasta*".fai")])
        sleep_exit("Fasta reference and/or index were not found ...")
    end
    
    # Check BAM index files exists
    ind_miss = false
    for bam in vcat(bams1,bams2)
        if !all([isfile(bam),isfile(bam*".bai")])
            ind_miss |= true
            print_log("BAM and/or its index file missing for")
            print_log("$(bam)")
        end
    end
    ind_miss && sleep_exit("Index files missing ...")

    # Individual bedGraph output files
    prefixes = [String(split(basename(bam),".")[1]) for bam in bams1]
    append!(prefixes,[String(split(basename(bam),".")[1]) for bam in bams2])
    out_mml_paths = ["$(outdir)/$(prefix)_mml.bedGraph" for prefix in prefixes]
    out_nme_paths = ["$(outdir)/$(prefix)_nme.bedGraph" for prefix in prefixes]
    
    # Check for existance of at least an output files
    if any(isfile.(vcat(out_mml_paths,out_nme_paths)))
        sleep_exit("At least an MML or NME file already exists...")
    end

    # Differentiial analysis bedGraph output files
    out_diff_paths = "$(outdir)/$(outprefix)_" .* ["tmml","tnme","tcmd"] .* "_diff_analysis.bedGraph"
    
    # Check for existance of at least an output files
    if any(isfile.(out_diff_paths))
        sleep_exit("At least a differential output file already exists...")
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
        chr_st = BED.chromstart(bed_rec) + (j-1) * size_anal_reg - 1
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
    isfile(bed_out) && (print_log("Found existing CpelTdm BED file ..."); return nothing)   

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
    reader = open(BED.Reader,bed)
    anal_regs = Vector{BED.Record}()
    while !eof(reader)
        # Get record
        # Note: when reading converts to a 1-based and start and end are inclusive
        record = BED.Record()
        read!(reader,record)
        # Keep if in chromosome
        BED.chrom(record)==chr && push!(anal_regs,record)
    end
    close(reader)

    # Return nothing
    return anal_regs
    
end
"""
    `get_cpg_pos!(ROI,FASTA)`

    Function that stores vector of CpG site positions in ROI.

    # Examples
    ```julia-repl
    julia> CpelTdm.get_cpg_pos!(roi,fasta)
    ```
"""
function get_cpg_pos!(roi_data::RoiData,fasta::String)::Nothing

    # Get sequence
    fasta_reader = open(FASTA.Reader,fasta,index=fasta*".fai")
    fa_record = fasta_reader[roi_data.chr]
    close(fasta_reader)

    # Get position CpG sites
    roi_seq = convert(String,FASTA.sequence(fa_record,roi_data.chrst:roi_data.chrend))
    cpg_pos = map(x->getfield(x,:offset),eachmatch(r"CG",roi_seq)) .+ roi_data.chrst .- 1
    
    # Set cpg_pos
    roi_data.cpg_pos = cpg_pos

    # Return nothing
    return nothing

end
"""
    `get_ns!(ROI,MAX_SUBREGION_SIZE)`

    Function that stores vector [N1,...,NK] given the position of the CpG sites & maximum subregion size.

    # Examples
    ```julia-repl
    julia> roi_data = CpelTdm.RoiData(10,10);
    julia> roi_data.chr = "chr22"; roi_data.chrst = 10; roi_data.chrend = 400;
    julia> roi_data.cpg_pos = [100,200,300,350];
    julia> max_size_subreg = 200;
    julia> CpelTdm.get_ns!(roi_data,max_size_subreg)
    julia> roi_data.n_vec
    2-element Array{Int64,1}:
    2
    2
    ```
"""
function get_ns!(roi_data::RoiData,max_size_subreg::Int64)::Nothing

    # Check if need to partition
    if (roi_data.chrend-roi_data.chrst)<=max_size_subreg
        roi_data.n_vec = [length(roi_data.cpg_pos)]
        return nothing
    end

    # Get K of model
    k = ceil((roi_data.chrend-roi_data.chrst+1)/max_size_subreg)

    # Get delimiters
    delimiters = vcat(collect(roi_data.chrst:(roi_data.chrend-roi_data.chrst+1)/k:roi_data.chrend),roi_data.chrend)

    # Partition CpG sites
    ids = [findfirst(y->y>=x,delimiters)-1 for x in roi_data.cpg_pos]

    # Set n vector
    roi_data.n_vec = [sum(ids.==i) for i in unique(ids)]

    # Return nothing
    return nothing

end
"""
    `get_obs_per_cpg(XOBS)`

    Function that returns the number of observations per CpG site as a vector.

    # Examples
    ```julia-repl
    julia> CpelTdm.get_obs_per_cpg(xobs)
    ```
"""
function get_obs_per_cpg(xobs::Array{Vector{Int64},1})::Vector{Int64}

    # Get observations per CpG site
    out = length(xobs)>0 ? vcat(sum(abs.(hcat(xobs...)),dims=2)...) : [0]

    # Return
    return out

end
"""
    `mean_cov(XOBS)`

    Function returns the average coverage per CpG given some observations `XOBS`.

    # Examples
    ```julia-repl
    julia> xobs=[[1,-1] for i=1:10]; append!(xobs,[[1,0] for i=1:10]);
    julia> CpelTdm.mean_cov(xobs)
    15.0
    ```
"""
function mean_cov(xobs::Array{Vector{Int64},1})::Float64

    # Return 0 if no observations
    return length(xobs)>0 ? norm(hcat(xobs...),1)/length(xobs[1]) : 0.0

end
"""
    `write_mml_out(OUT_PMAP,MML_PATHS)`

    Function that writes MML in corresponding bedGraph files.

    # Examples
    ```julia-repl
    julia> CpelTdm.write_mml_out(out_pmap,mml_paths)
    ```
"""
function write_mml_out(out_pmap::Vector{RoiData},mml_paths::Vector{String})::Nothing
    
    # Open streams
    ios = [open(path,"a") for path in mml_paths]
    
    # Loop over ROIs
    @inbounds for roi in out_pmap
        
        # Get data from ROI
        mmls = vcat(roi.mmls1,roi.mmls2)
        
        # Write
        @inbounds for i=1:length(ios)
            write(ios[i],"$(roi.chr)\t$(roi.chrst-1)\t$(roi.chrend)\t$(mmls[i])\n")
        end

    end

    # Close all streams
    close.(ios)

    # Return nothing
    return nothing

end
"""
    `write_nme_out(OUT_PMAP,NME_PATHS)`

    Function that writes NME in corresponding bedGraph files.

    # Examples
    ```julia-repl
    julia> CpelTdm.write_nme_out(out_pmap,nme_paths)
    ```
"""
function write_nme_out(out_pmap::Vector{RoiData},nme_paths::Vector{String})::Nothing
    
    # Open streams
    ios = [open(path,"a") for path in nme_paths]
    
    # Loop over ROIs
    @inbounds for roi in out_pmap

        # Get data from ROI
        nmes = vcat(roi.nmes1,roi.nmes2)
        
        # Write
        @inbounds for i=1:length(ios)
            write(ios[i],"$(roi.chr)\t$(roi.chrst-1)\t$(roi.chrend)\t$(nmes[i])\n")
        end

    end

    # Close all streams
    close.(ios)

    # Return nothing
    return nothing

end
"""
    `write_diff_out(OUT_PMAP,DIFF_PATHS)`

    Function that writes differential output in corresponding bedGraph files.

    # Examples
    ```julia-repl
    julia> CpelTdm.write_diff_out(out_pmap,diff_paths)
    ```
"""
function write_diff_out(out_pmap::Vector{RoiData},diff_paths::Vector{String})::Nothing
    
    # Open streams
    ios = [open(path,"a") for path in diff_paths]

    # Loop over ROIs
    @inbounds for roi in out_pmap
        
        # Get data from ROI
        tmml,pmml = roi.mml_test
        tnme,pnme = roi.nme_test
        tcmd,pcmd = roi.cmd_test

        # Write
        write(ios[1],"$(roi.chr)\t$(roi.chrst-1)\t$(roi.chrend)\t$(tmml)\t$(pmml)\n")
        write(ios[2],"$(roi.chr)\t$(roi.chrst-1)\t$(roi.chrend)\t$(tnme)\t$(pnme)\n")
        write(ios[3],"$(roi.chr)\t$(roi.chrst-1)\t$(roi.chrend)\t$(tcmd)\t$(pcmd)\n")

    end

    # Close all streams
    close.(ios)

    # Return nothing
    return nothing

end
"""
    `write_output(OUT_PMAP,OUT_PATHS)`

    Function that write output of pmap into respective files.

    # Examples
    ```julia-repl
    julia> CpelTdm.write_output(out_pmap,out_paths)
    ```
"""
function write_output(out_pmap::Vector{RoiData},out_paths::NTuple{3,Vector{String}})::Nothing

    # Assign paths 
    mml_paths = out_paths[1]
    nme_paths = out_paths[2]
    diff_paths = out_paths[3]

    # Write MML
    write_mml_out(out_pmap,mml_paths)

    # Write NMEs
    write_nme_out(out_pmap,nme_paths)

    # Write differential analysis
    write_diff_out(out_pmap,diff_paths)

    # Return nothing
    return nothing

end
"""
    `add_qvalues_to_path(PATH)`

    Function that takes in the differential output and adds BH q-value.

    # Examples
    ```julia-repl
    julia> CpelAsm.add_qvalues_to_path(path)
    ```
"""
function add_qvalues_to_path(path::String)::Nothing

    # Leave if no data
    filesize(path)>0 || return nothing
    
    # Get data
    all_data = readdlm(path,'\t',Any)
    qvals = fill(NaN,size(all_data)[1])
    
    # Multiple hypothesis testing correction
    ind = .!isnan.(all_data[:,5])
    if sum(ind)>0 
        qvals[ind] = MultipleTesting.adjust(convert(Vector{Float64},all_data[ind,5]),BenjaminiHochberg())
    end
    
    # Append to output matrix
    all_data = hcat(all_data,qvals)
    
    # Write to temp output
    temp_path = path * ".tmp"
    open(temp_path,"w") do io
        writedlm(io,all_data,'\t')
    end

    # Move to original file
    mv(temp_path,path,force=true)

    # Return
    return nothing

end
"""
    `mult_hyp_corr(DIFF_PATHS)`

    Function that takes in all the differential output and adds BH q-value in each one.

    # Examples
    ```julia-repl
    julia> CpelAsm.mult_hyp_corr(diff_paths)
    ```
"""
function mult_hyp_corr(diff_paths::Vector{String})

    # Add q-values
    add_qvalues_to_path(diff_paths[1])
    add_qvalues_to_path(diff_paths[2])
    add_qvalues_to_path(diff_paths[3])

    # Return
    return nothing

end
"""
    `cpel_tdm(BAMS1,BAMS2,BED,FASTA,OUTDIR,PREFIX)`

    Function to call to CpelTdm analysis. 

        - BAMS1: BAM files associated to group 1 (require index file).
        - BAMS2: BAM files associated to group 2 (require index file).
        - BED: BED file containing regions of interest.
        - FASTA: reference genome (requires index file).
        - OUTDIR: output directory to store results.
        - PREFIX: prefix of output differential files.

    # Examples
    ```julia-repl
    julia> CpelTdm.cpel_tdm(bams1,bams2,bed,fasta,outdir,prefix)
    ```
"""
function cpel_tdm(bams1::Vector{String},bams2::Vector{String},bed::String,fasta::String,outdir::String,prefix::String;
                  pe::Bool=false,max_size_anal_reg::Int64=1000,max_size_subreg::Int64=250,min_cov::Int64=5,
                  matched::Bool=false,trim::NTuple{4,Int64}=(0,0,0,0),bound_check::Bool=false)::Nothing

    # Print initialization of juliASM
    print_log("Starting CpelTdm ...")

    ## IO
    print_log("Checking IO ...")
    
    # Get paths
    out_paths = get_paths(bams1,bams2,fasta,bed,outdir,prefix)

    # Create output folder if it doesn't exist
    isdir(outdir) || mkdir(outdir)

    ## Configure run
    print_log("Configuring CpelTdm ...")
    matched ? print_log("Matched comparison ...") : print_log("Unmatched comparison ...")
    config = CpeltdmConfig(pe,min_cov,matched,bound_check,trim,max_size_subreg,max_size_anal_reg)

    # Check same number of samles is matched
    matched && (length(bams1)!=length(bams2)) && sleep_exit("Unmatched samples found. Exiting ...")
    
    ## Define analysis regions
    print_log("Defining analysis regions from BED file ...")

    # Divide BED file divided into analysis regions
    bed_out = get_bed_out_filename(bed)
    gen_anal_reg_file(bed,bed_out,max_size_anal_reg)

    ## Inference

    # Estimate θ, compute MML & NME, and do test for each ROI
    print_log("Running differential analysis ...")
    anal_bed_file(bams1,bams2,bed_out,fasta,out_paths,config)
    
    ## Done

    # Print done
    print_log("Done.")

    # Return
    return nothing

end
"""
    `anal_bed_file(BAMS1,BAMS2,BED,FASTA,OUT_PATHS,CONFIG)`

    Function that estimates θ for all samples, computes statistics, and does a matched/unmatched 
    permutation test at each ROI.

    # Examples
    ```julia-repl
    julia> CpelTdm.anal_bed_file(bams1,bams2,bed,fasta,out_paths,config)
    ```
"""
function anal_bed_file(bams1::Vector{String},bams2::Vector{String},bed::String,fasta::String,
                       out_paths::NTuple{3,Vector{String}},config::CpeltdmConfig)::Nothing

    # Find chromosomes
    fasta_reader = open(FASTA.Reader,fasta,index=fasta*".fai")
    chr_sizes = fasta_reader.index.lengths
    chr_names = fasta_reader.index.names
    close(fasta_reader)

    # Loop over chromosomes
    for chr in chr_names

        # Print info
        print_log("Processing chr: $(chr) ...")

        # Get windows pertaining to current chromosome
        print_log("Obtaining corresponding ROIs ...")
        roi_chr = get_anal_reg_chr(bed,chr)
        chr_size = chr_sizes[findfirst(x->x==chr,chr_names)]
        length(roi_chr)>0 || continue

        # Process regions of interest in chromosome
        print_log("Performing differential analysis on ROIs ...")
        out_pmap = pmap(x->pmap_anal_roi(x,chr,chr_size,bams1,bams2,fasta,config),roi_chr)

        # Add last to respective bedGraph file
        length(out_pmap)>0 && write_output(out_pmap,out_paths)

    end

    # Correct for multiple hypothesis in differential analysis
    print_log("Performing multiple hypothesis correction ...")
    mult_hyp_corr(out_paths[3])

    # Return nothing
    return nothing

end
"""
    `pmap_anal_roi(ROI,CHR,CHR_SIZE,BAMS1,BAMS2,FASTA,CONFIG)`

    Function that estimates θ and computes MML/NME for all samples, and does permutation test for a given ROI.

    # Examples
    ```julia-repl
    julia> CpelTdm.pmap_anal_roi(roi,chr,chr_size,bams1,bams2,fasta,config)
    ```
"""
function pmap_anal_roi(roi::BED.Record,chr::String,chr_size::Int64,bams1::Vector{String},bams2::Vector{String},
                       fasta::String,config::CpeltdmConfig)::RoiData

    # Print ROI for testing
    # println(roi)

    # Get sample size
    s1 = length(bams1)
    s2 = length(bams2)

    # Empty output
    roi_data = RoiData(s1,s2)

    # Store chromosomal info
    roi_data.chr = BED.chrom(roi)
    roi_data.chrst = BED.chromstart(roi)
    roi_data.chrend = BED.chromend(roi)

    # Get CpG sites
    get_cpg_pos!(roi_data,fasta)
    length(roi_data.cpg_pos)>0 || return roi_data

    # Get vector of Ns
    get_ns!(roi_data,config.max_size_subreg)

    ## Estimate θs

    # Loop over group 1 samples
    # print_log("First group ...")
    @inbounds for i=1:s1

        # Get BAM file
        bam = bams1[i]

        # Get vectors from BAM overlapping region
        xobs = read_bam(bam,chr,roi_data.chrst,roi_data.chrend,roi_data.cpg_pos,chr_size,config.pe,config.trim)
        obs_per_cpg = get_obs_per_cpg(xobs)
        (config.min_cov<=mean_cov(xobs)<=400) && (sum(obs_per_cpg.==0)<=1.0/5.0*sum(roi_data.n_vec)) || continue

        # Estimate parameters of CPEL models
        θhat = est_theta_sa(roi_data.n_vec,xobs)
        config.bound_check && check_boundary(θhat) && continue
        roi_data.θ1s[i] = θhat

        # Estimate intermediate quantities
        ∇logZ = get_∇logZ(roi_data.n_vec,θhat)

        # Compute output
        roi_data.mmls1[i] = comp_mml(roi_data.n_vec,∇logZ)
        roi_data.nmes1[i] = comp_nme(roi_data.n_vec,θhat,∇logZ)

        # Set i-th sample done
        roi_data.analyzed1[i] = true

    end

    # Loop over group 2 samples
    # print_log("Second group ...")
    @inbounds for i=1:s2

        # Get BAM file
        bam = bams2[i]

        # Get vectors from BAM overlapping region
        xobs = read_bam(bam,chr,roi_data.chrst,roi_data.chrend,roi_data.cpg_pos,chr_size,config.pe,config.trim)
        obs_per_cpg = get_obs_per_cpg(xobs)
        (config.min_cov<=mean_cov(xobs)<=400) && (sum(obs_per_cpg.==0)<=1.0/5.0*sum(roi_data.n_vec)) || continue

        # Estimate parameters of CPEL models
        θhat = est_theta_sa(roi_data.n_vec,xobs)
        config.bound_check && check_boundary(θhat) && continue
        roi_data.θ2s[i] = θhat

        # Estimate intermediate quantities
        ∇logZ = get_∇logZ(roi_data.n_vec,θhat)

        # Compute output
        roi_data.mmls2[i] = comp_mml(roi_data.n_vec,∇logZ)
        roi_data.nmes2[i] = comp_nme(roi_data.n_vec,θhat,∇logZ)

        # Set i-th sample done
        roi_data.analyzed2[i] = true

    end

    ## Differential analysis
    # print_log("Differential analysis ...")
    if config.matched
        
        # Find pairs with data. TODO: issue with -1 power
        ind = roi_data.analyzed1 .& roi_data.analyzed2
        (0.5^(sum(ind)-1))<0.05 || return roi_data
        
        # Get sample pairs with data
        θ1s = roi_data.θ1s[ind]
        θ2s = roi_data.θ2s[ind]
        
        # Perform matched differential analysis
        tmml_test,tnme_test,tcmd_test = mat_tests(roi_data.n_vec,θ1s,θ2s)

    else
        
        # Check if sufficient data
        s1 = sum(roi_data.analyzed1)
        s2 = sum(roi_data.analyzed2)
        1/binomial(s1+s2,s1)<0.05 || return roi_data
        
        # Get samples with data
        θ1s = roi_data.θ1s[roi_data.analyzed1]
        θ2s = roi_data.θ2s[roi_data.analyzed2]
        
        # Perform unmatched differential analysis
        tmml_test,tnme_test,tcmd_test = unmat_tests(roi_data.n_vec,θ1s,θ2s)    
        
    end
    
    # Set differential analysis done
    roi_data.diff_anal_done = true

    # Store results
    roi_data.mml_test = tmml_test
    roi_data.nme_test = tnme_test
    roi_data.cmd_test = tcmd_test

    # Return
    return roi_data

end