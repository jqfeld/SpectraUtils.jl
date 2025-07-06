using SpectraUtils
using Documenter

DocMeta.setdocmeta!(SpectraUtils, :DocTestSetup, :(using SpectraUtils); recursive=true)

makedocs(;
    modules=[SpectraUtils],
    authors="Jan Kuhfeld <jan@kuhfeld.net> and contributors",
    sitename="SpectraUtils.jl",
    format=Documenter.HTML(;
        canonical="https://jqfeld.github.io/SpectraUtils.jl",
        edit_link="main",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)

deploydocs(;
    repo="github.com/jqfeld/SpectraUtils.jl",
    devbranch="main",
)
