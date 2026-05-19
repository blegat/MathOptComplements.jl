using Documenter
using DocumenterInterLinks
using MathOptComplements

links = InterLinks("MathOptInterface" => "https://jump.dev/MathOptInterface.jl/stable/")

makedocs(
    sitename = "MathOptComplements.jl",
    format = Documenter.HTML(
        assets = ["assets/favicon.ico"],
        prettyurls = Base.get(ENV, "CI", nothing) == "true",
        mathengine = Documenter.KaTeX(),
    ),
    modules = [MathOptComplements],
    repo = "https://github.com/jump-dev/MathOptComplements.jl/blob/{commit}{path}#{line}",
    checkdocs = :none,
    clean = true,
    pages = [
        "Home" => "index.md",
        "Quickstart" => "quickstart.md",
        "Tutorials" => ["Equilibrium problem" => "equilibrium.md"],
        "API reference" => [
            "Reformulations" => "api/reformulation.md",
            "Relaxations" => "api/relaxation.md",
            "Bridges" => "api/bridges.md",
        ],
    ],
    plugins = [links],
)

deploydocs(
    repo = "github.com/jump-dev/MathOptComplements.jl.git",
    target = "build",
    devbranch = "main",
    devurl = "dev",
    push_preview = true,
)
