using Documenter
using Crystalline

# make sure we actually do `using Crystalline` before calling doctests (as described in
# https://juliadocs.github.io/Documenter.jl/stable/man/doctests/#Module-level-metadata)
DocMeta.setdocmeta!(Crystalline, :DocTestSetup, :(using Crystalline); recursive=true)

makedocs(
    modules = [Crystalline],
    sitename = "Crystalline",
    format = Documenter.HTML(prettyurls = get(ENV, "CI", nothing) == "true"),
    authors = "Thomas Christensen and contributors",
    pages   = [
        "Home" => "index.md",
        "API"  => "api.md"
    ]
)

# Documenter can also automatically deploy documentation to gh-pages.
# See "Hosting Documentation" and deploydocs() in the Documenter manual
# for more information.
deploydocs(
    repo   = "github.com/thchr/Crystalline.jl.git",
    target = "build",
    deps   = nothing,
    make   = nothing,
    push_preview = true
)
