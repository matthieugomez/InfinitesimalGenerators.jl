using Documenter
using InfinitesimalGenerators
import Plots

# Render GR plots without a display (needed on CI), and give every figure enough margin that
# the axis labels are never clipped.
ENV["GKSwstype"] = "100"
Plots.default(; left_margin = 5Plots.mm, bottom_margin = 5Plots.mm)

makedocs(;
    sitename = "InfinitesimalGenerators.jl",
    authors = "Matthieu Gomez",
    modules = [InfinitesimalGenerators],
    format = Documenter.HTML(;
        prettyurls = get(ENV, "CI", "false") == "true",
        edit_link = "main",
        assets = ["assets/custom.css"],
    ),
    checkdocs = :none,
    pagesonly = true,
    pages = [
        "Home" => "index.md",
        "Manual" => [
            "Univariate processes"   => "univariate.md",
            "Multivariate processes" => "multivariate.md",
            "Operators"              => "operators.md",
        ],
        "Tutorials" => [
            "Expected values (Feynman–Kac)" => "expectations.md",
            "Distribution dynamics"         => "distributions.md",
            "Solving HJB equations"         => "hjb.md",
            "Tail indices"                  => "tail_index.md",
        ],
        "API reference" => "api.md",
    ],
)

deploydocs(; repo = "github.com/matthieugomez/InfinitesimalGenerators.jl.git", devbranch = "main", push_preview = true)
