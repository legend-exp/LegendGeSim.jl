# Use
#
#     DOCUMENTER_DEBUG=true julia --color=yes make.jl local [nonstrict] [fixdoctests]
#
# for local builds.

using Documenter
using LegendGeSim

# Doctest setup
DocMeta.setdocmeta!(
    LegendGeSim,
    :DocTestSetup,
    :(using LegendGeSim);
    recursive=true,
)

makedocs(
    sitename = "LegendGeSim",
    modules = [LegendGeSim],
    format = Documenter.HTML(
        prettyurls = !("local" in ARGS),
        canonical = "https://legend-exp.github.io/LegendGeSim.jl/stable/"
    ),
    pages = [
        "Home" => "index.md",
        "API" => "api.md",
        "LICENSE" => "LICENSE.md",
    ],
    doctest = ("fixdoctests" in ARGS) ? :fix : true,
    linkcheck = !("nonstrict" in ARGS),
    warnonly = ("nonstrict" in ARGS),
)

deploydocs(
    repo = "github.com/legend-exp/LegendGeSim.jl.git",
    forcepush = true,
    push_preview = true,
)
