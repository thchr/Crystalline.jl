using Documenter
using Crystalline

# make sure we actually do `using Crystalline` before calling doctests (as described in
# https://juliadocs.github.io/Documenter.jl/stable/man/doctests/#Module-level-metadata)
DocMeta.setdocmeta!(Crystalline, :DocTestSetup, :(using Crystalline); recursive=true)

makedocs(
    modules = [Crystalline, Bravais],
    sitename = "Crystalline.jl",
    authors = "Thomas Christensen <tchr@mit.edu> and contributors",
    repo = "https://github.com/thchr/Crystalline.jl/blob/{commit}{path}#L{line}",
    format=Documenter.HTML(;
        prettyurls = get(ENV, "CI", nothing) == "true",
        canonical = "https://thchr.github.io/Crystalline.jl",
        assets = ["assets/custom.css"], # increase logo size
    ),
    pages = [
        "Home"                  => "index.md",
        "Symmetry operations"   => "operations.md",
        "Groups"                => "groups.md",
        "Irreps"                => "irreps.md",
        "Bravais types & bases" => "bravais.md",
        "Band representations"  => "bandreps.md",
        "Lattices"              => "lattices.md",
        "API"                   => "api.md",
        "Internal API"          => "internal-api.md",
    ],
    warnonly = Documenter.except(
        :autodocs_block, :cross_references, :docs_block, :doctest, :eval_block, 
        :example_block, :footnote, :linkcheck_remotes, :linkcheck, :meta_block, 
        :parse_error, :setup_block,
        #:missing_docs # necessary due to docstrings from SmithNormalForm vendoring :(
    ),
    clean = true
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
