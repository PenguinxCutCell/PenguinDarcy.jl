using Documenter
using PenguinDarcy

makedocs(
    modules = [PenguinDarcy],
    authors = "PenguinxCutCell contributors",
    sitename = "PenguinDarcy.jl",
    format = Documenter.HTML(
        canonical = "https://PenguinxCutCell.github.io/PenguinDarcy.jl",
        repolink = "https://github.com/PenguinxCutCell/PenguinDarcy.jl",
        collapselevel = 2,
    ),
    pages = [
        "Home" => "index.md",
        "API" => "api.md",
        "Examples" => "examples.md",
        "Algorithms" => "algorithms.md",
        "Darcy Model" => "darcy.md",
    ],
    pagesonly = true,
    warnonly = false,
    remotes = nothing,
)

if get(ENV, "CI", "") == "true"
    deploydocs(
        repo = "github.com/PenguinxCutCell/PenguinDarcy.jl",
        push_preview = true,
    )
end
