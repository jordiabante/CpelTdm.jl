using Documenter,CpelTdm

makedocs(
    sitename="CpelTdm",
    pages = [
        "Home"           => "index.md",
        "Main Commands"  => "main_commands.md",
        "Optimal BED"    => "optimal_bed.md",
        "Toy Example"    => "toy_example.md"
    ],
    authors = "Jordi Abante"
)

deploydocs(
    repo = "github.com/jordiabante/CpelTdm.jl.git",
)