using GaPLAC
using Documenter

makedocs(;
    modules=[GaPLAC],
    authors="Kevin Bonham, PhD; Jason Lloyde-Price, PhD",
    repo="https://github.com/biobakery/GaPLAC/blob/{commit}{path}#L{line}",
    sitename="GaPLAC",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://biobakery.github.io/GaPLAC",
        assets=String[]
    ),
    pages=[
        "Home" => "index.md"
    ],
)

deploydocs(;
    repo="github.com/biobakery/GaPLAC",
)
