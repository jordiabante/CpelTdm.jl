###################################################################################################
# IO FUNCTIONS
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
    out_diff_paths = "$(outdir)/$(outprefix)_" .* ["tmml","tnme","tpdm"] .* "_diff_analysis.bedGraph"
    
    # Check for existance of at least an output files
    if any(isfile.(out_diff_paths))
        sleep_exit("At least a differential output file already exists...")
    end

    # Return paths
    return out_mml_paths,out_nme_paths,out_diff_paths
    
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
        tpdm,ppdm = roi.pdm_test

        # Write
        write(ios[1],"$(roi.chr)\t$(roi.chrst-1)\t$(roi.chrend)\t$(tmml)\t$(pmml)\n")
        write(ios[2],"$(roi.chr)\t$(roi.chrst-1)\t$(roi.chrend)\t$(tnme)\t$(pnme)\n")
        write(ios[3],"$(roi.chr)\t$(roi.chrst-1)\t$(roi.chrend)\t$(tpdm)\t$(ppdm)\n")

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