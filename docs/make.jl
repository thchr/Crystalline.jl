using Documenter
using Crystalline

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
