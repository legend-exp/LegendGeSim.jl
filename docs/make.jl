# Use
#
#     DOCUMENTER_DEBUG=true julia --color=yes make.jl local [nonstrict] [fixdoctests]
#
# for local builds.

using Documenter
using LegendGeSim
using Literate

# Doctest setup
DocMeta.setdocmeta!(
    LegendGeSim,
    :DocTestSetup,
    :(using LegendGeSim);
    recursive=true,
)

function fix_literate_output(content)
    content = replace(content, "EditURL = \"@__REPO_ROOT_URL__/\"" => "")
    return content
end

gen_content_dir = joinpath(@__DIR__, "src", "tutorials")
for tut_lit_fn in filter(fn -> endswith(fn, "_lit.jl"), readdir(gen_content_dir))
    lit_src_fn = joinpath(gen_content_dir, tut_lit_fn)
    tut_basename = tut_lit_fn[1:end-7] # remove "_lit.jl"
    Literate.markdown(lit_src_fn, gen_content_dir, name = tut_basename, documenter = true, credit = true, postprocess = fix_literate_output)
    Literate.notebook(lit_src_fn, gen_content_dir, execute = false, name = tut_basename, documenter = true, credit = true)
    Literate.script(lit_src_fn, gen_content_dir, keep_comments = false, name = tut_basename, documenter = true, credit = false)
end

makedocs(
    sitename = "LegendGeSim",
    modules = [LegendGeSim],
    format = Documenter.HTML(
        prettyurls = !("local" in ARGS),
        canonical = "https://legend-exp.github.io/LegendGeSim.jl/stable/"
    ),
    pages = [
        "Home" => "index.md",
        "Manual" => Any[
            "The LEGEND Germanium Simulation Chain" => "man/LegendGeSim.md",
            "Field Simulation" => "man/field_simulation.md",
            "Ideal Pulse Simulation" => "man/ideal_pulse_simulation.md",
            "Realistic Waveform Simulation" => "man/realistic_waveform_simulation.md"
        ],
        "Tutorials" => [
            "Detector Geometry" => "tutorials/detector_geometry.md",
            "Simulate Fields" => "tutorials/simulate_fields.md",
            "Simulate Ideal Pulses" => "tutorials/simulate_ideal_pulses.md",
            "Simulate Realistic Waveforms" => "tutorials/simulate_realistic_waveforms.md"
        ],
        "API" => "api.md",
        "LICENSE" => "LICENSE.md",
    ],
    doctest = ("fixdoctests" in ARGS) ? :fix : true,
    linkcheck = !("nonstrict" in ARGS),
    warnonly = ("nonstrict" in ARGS),
)

deploydocs(
    repo = "github.com/legend-exp/LegendGeSim.jl.git",
    devbranch = "dev",
    devurl = "dev",
    forcepush = true,
    push_preview = true,
)
