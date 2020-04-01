using gptool
using Documenter

makedocs(;
    modules=[GPTool],
    authors="",
    repo="https://github.com/biobakery/gptool.jl/blob/{commit}{path}#L{line}",
    sitename="gptool.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://biobakery.github.io/gptool.jl",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)

deploydocs(;
    repo="github.com/biobakery/gptool.jl",
)