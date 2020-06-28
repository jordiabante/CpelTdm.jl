var documenterSearchIndex = {"docs":
[{"location":"main_commands/#","page":"Main Commands","title":"Main Commands","text":"(Image: CI) (Image: Docs) (Image: MIT license)","category":"page"},{"location":"main_commands/#Main-commands-1","page":"Main Commands","title":"Main commands","text":"","category":"section"},{"location":"main_commands/#Unmatched-group-comparison-1","page":"Main Commands","title":"Unmatched group comparison","text":"","category":"section"},{"location":"main_commands/#","page":"Main Commands","title":"Main Commands","text":"In order to perform an unmatched group comparison, we first need to load the package. Then, the right paths need to be specified, including paths to the reference FASTA file, given by fa below, the path to the BED file containing the regions of interest, given by bed below, and the path to each bam file for each groups, given by the bams1 & bams2 vectors below. Note that the BED file is expected to have for columns, being the first three chr, st,  and end. In addition, the output directory outdir and prefix  need to be defined as well. Finally, we call function cpel_tdm  setting the optional argument matched equal to false, since in  this case we wish to perform an unmatched group comparison. In  addition, here we assume that we are using paired end sequencing  data, hence we set pe=true. However, if the sequencing data were  single end, we would set pe=false instead.","category":"page"},{"location":"main_commands/#","page":"Main Commands","title":"Main Commands","text":"\n# Load dependencies\nusing Distributed\n@everywhere using CpelTdm\n\n# Set paths\ndir_name = \"/path/to/data/\"\nfa = \"$(dir_name)/fasta/reference.fa\"\nbed = \"$(dir_name)/bed/regions_of_interest.bed\"\nbams1 = \"$(dir_name)/bam/\".*[\"g1_s1.bam\",\"g1_s2.bam\",\"g1_s3.bam\",\"g1_s4.bam\"]\nbams2 = \"$(dir_name)/bam/\".*[\"g2_s1.bam\",\"g2_s2.bam\",\"g2_s3.bam\",\"g2_s4.bam\",\"g2_s5.bam\"]\n\n# Output paths\nprefix = \"cpeltdm\"\noutdir = \"$(dir_name)/cpeltdm-unmatched/\"\n\n# Run matched analysis\ncpel_tdm(bams1,bams2,bed,fa,outdir,prefix;pe=true,matched=false)\n","category":"page"},{"location":"main_commands/#","page":"Main Commands","title":"Main Commands","text":"If you wish to run CpelTdm on a cluster, we recommend you store the previous code into a file, say cpeltdm_unmatched_call.jl.  Then, in the submission script you can execute the content of  this file using the command julia -p 10 cpeltdm_unmatched_call.jl,  which will execute the code providing CpelTdm with 10 CPUs to  perform the analysis.","category":"page"},{"location":"main_commands/#Matched-group-comparison-1","page":"Main Commands","title":"Matched group comparison","text":"","category":"section"},{"location":"main_commands/#","page":"Main Commands","title":"Main Commands","text":"In order to perform an matched group comparison we need to follow the same steps as in the unmatched group comparison as far as  the paths go. In this case, however, CpelTdm assumes that the i-th sample in bams1 matches that of bams2. Finally, we call  function cpel_tdm setting the optional argument matched equal  to true, since in this case we wish to perform a matched group  comparison. In addition, here we assume that we are using paired  end sequencing data, hence we set pe=true. ","category":"page"},{"location":"main_commands/#","page":"Main Commands","title":"Main Commands","text":"\n# Load dependencies\nusing Distributed\n@everywhere using CpelTdm\n\n# Set paths\ndir_name = \"/path/to/data/\"\nfa = \"$(dir_name)/fasta/reference.fa\"\nbed = \"$(dir_name)/bed/regions_of_interest.bed\"\nbams1 = \"$(dir_name)/bam/\".*[\"g1_s1.bam\",\"g1_s2.bam\",\"g1_s3.bam\",\"g1_s4.bam\",\"g1_s5.bam\"]\nbams2 = \"$(dir_name)/bam/\".*[\"g2_s1.bam\",\"g2_s2.bam\",\"g2_s3.bam\",\"g2_s4.bam\",\"g2_s5.bam\"]\n\n# Output paths\nprefix = \"cpeltdm\"\noutdir = \"$(dir_name)/cpeltdm-matched/\"\n\n# Run matched analysis\ncpel_tdm(bams1,bams2,bed,fa,outdir,prefix;pe=true,matched=true)\n","category":"page"},{"location":"main_commands/#","page":"Main Commands","title":"Main Commands","text":"Similar to the unmatched case, if you wish to run CpelTdm on a  cluster, we recommend you store the previous code into a file,  say cpeltdm_matched_call.jl. Then, in the submission script you  can execute the content of this file using the command  julia -p 10 cpeltdm_matched_call.jl, which will execute the  code providing CpelTdm with 10 CPUs to perform the analysis.","category":"page"},{"location":"toy_example/#","page":"Toy Example","title":"Toy Example","text":"(Image: CI) (Image: Docs) (Image: MIT license)","category":"page"},{"location":"toy_example/#Toy-Example-1","page":"Toy Example","title":"Toy Example","text":"","category":"section"},{"location":"toy_example/#","page":"Toy Example","title":"Toy Example","text":"The package includes a small toy example consisting of an unmatched group comparison, which should take approximately 20 seconds to run. In particular,  samples in group 1 were generated under an uncorrelated CPEL model with μ1=0.9, while samples in group 2 were generated under an uncorrelated CPEL model with  μ2=0.1. Thus, the example in particular consists of a region where the two groups  strongly differ in mean methylation level (MML) and probability distribution of  methylation (PDM), while they do not show any normalized methylation entropy  (NME) difference. To run the toy example, go to the package folder and simply run  julia test/toy_example.jl. CpelTdm will generate a total of 35 output files consisting  of:","category":"page"},{"location":"toy_example/#","page":"Toy Example","title":"Toy Example","text":"A file per sample containing MML μ (16 total).\nA file per sample containing NME h (16 total).\nA file per differential analysis type (3 total).","category":"page"},{"location":"toy_example/#","page":"Toy Example","title":"Toy Example","text":"in folder test/cpeltdm/.","category":"page"},{"location":"#","page":"Home","title":"Home","text":"(Image: CI) (Image: Docs) (Image: MIT license)","category":"page"},{"location":"#Description-1","page":"Home","title":"Description","text":"","category":"section"},{"location":"#","page":"Home","title":"Home","text":"CpelTdm is a julia package designed to perform DNA methylation differential analysis.  The package is based on a recently published method [1], and has several advantages when  compared to traditional DNA methylation differential analysis packages. CpelTdm not only  detects mean methylation level differences between groups, but it can also detect  significant differences in methylation entropy as well as in the probability distribution of methylation (see [1] for technical details).","category":"page"},{"location":"#Testing-1","page":"Home","title":"Testing","text":"","category":"section"},{"location":"#","page":"Home","title":"Home","text":"CpelTdm is tested against Julia 1.3.0 on the latest versions of Linux, macOS and Windows.","category":"page"},{"location":"#Getting-Started-1","page":"Home","title":"Getting Started","text":"","category":"section"},{"location":"#Prerequisites-1","page":"Home","title":"Prerequisites","text":"","category":"section"},{"location":"#","page":"Home","title":"Home","text":"julia v1.3.0\ngit.","category":"page"},{"location":"#Installing-1","page":"Home","title":"Installing","text":"","category":"section"},{"location":"#","page":"Home","title":"Home","text":"CpelTdm and dependencies can be installed via the following command in julia's REPL:","category":"page"},{"location":"#","page":"Home","title":"Home","text":"(v1.3) pkg> add https://github.com/jordiabante/CpelTdm.jl.git","category":"page"},{"location":"#Running-the-tests-1","page":"Home","title":"Running the tests","text":"","category":"section"},{"location":"#","page":"Home","title":"Home","text":"In order to test the package to ensure that it has been properly installed, run the following command in a julia session:","category":"page"},{"location":"#","page":"Home","title":"Home","text":"(v1.3) pkg> test CpelTdm","category":"page"},{"location":"#","page":"Home","title":"Home","text":"If the package has been properly installed, then all tests will be successful.","category":"page"},{"location":"#Authors-1","page":"Home","title":"Authors","text":"","category":"section"},{"location":"#","page":"Home","title":"Home","text":"Jordi Abante","category":"page"},{"location":"#License-1","page":"Home","title":"License","text":"","category":"section"},{"location":"#","page":"Home","title":"Home","text":"This project is licensed under the MIT License - see the LICENSE.md file for details.","category":"page"},{"location":"#References-1","page":"Home","title":"References","text":"","category":"section"},{"location":"#","page":"Home","title":"Home","text":"[1] Abante, J., Goutsias, J., CpelTdm.jl: a Julia package for targeted differential  DNA methylation analysis, Bioinformatics 2020 XYZ.","category":"page"}]
}
