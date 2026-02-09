using Geodynamics
using Documenter

DocMeta.setdocmeta!(Geodynamics, :DocTestSetup, :(using Geodynamics); recursive = true)

makedocs(;
    modules = [Geodynamics],
    authors = "EarthSciML authors and contributors",
    repo = "https://github.com/EarthSciML/Geodynamics.jl/blob/{commit}{path}#{line}",
    sitename = "Geodynamics.jl",
    format = Documenter.HTML(;
        prettyurls = get(ENV, "CI", "false") == "true",
        canonical = "https://geodynamics.earthsci.dev",
        assets = String[],
        repolink = "https://github.com/EarthSciML/Geodynamics.jl",
    ),
    pages = [
        "Home" => "index.md",
        "Mogi & McTigue Models" => "mogi_mctigue.md",
        "API" => "api.md",
    ],
)

deploydocs(; repo = "github.com/EarthSciML/Geodynamics.jl.git")
