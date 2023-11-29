using GeomOpt
using Documenter

DocMeta.setdocmeta!(GeomOpt, :DocTestSetup, :(using GeomOpt); recursive=true)

makedocs(;
    modules=[GeomOpt],
    authors="CedricTravelletti <cedrictravelletti@gmail.com> and contributors",
    repo="https://github.com/CedricTravelletti/GeomOpt.jl/blob/{commit}{path}#{line}",
    sitename="GeomOpt.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://CedricTravelletti.github.io/GeomOpt.jl",
        edit_link="main",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)

deploydocs(;
    repo="github.com/CedricTravelletti/GeomOpt.jl",
    devbranch="main",
)
