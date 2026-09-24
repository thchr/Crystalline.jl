# Package Crystalline's data sets into the tarballs that back its Julia artifacts, and
# print everything needed to publish them. `build/DATA-RELEASE.md` is the procedure this
# script belongs to; read it first.
#
#   julia --project=build build/data_release.jl <tag> [<name>...]
#
# e.g. `julia --project=build build/data_release.jl data-v0.0.2 isotropy`. With no names,
# every data set is packaged.

using Tar, SHA

const DATA_DIR  = joinpath(dirname(@__DIR__), "data")
const STAGE_DIR = joinpath(@__DIR__, "data-release")
const REPO      = "thchr/Crystalline.jl"

# Each data set becomes one artifact and one `<name>.tar.gz` release asset. Source paths
# are given relative to `data/`, and land at the root of the artifact under their own file
# name: an artifact is its own namespace, so there is nothing to be gained by rebuilding the
# directories a file happens to sit in here. A data set whose files do need a structure of
# their own — several dimensions or space groups, say — gives `source => path-in-artifact`
# pairs instead.
const DATASETS = Dict(
    "bilbao_spinless_irreps" => ["irreps/lgs/3d/irreps_data_spinless_bilbao.jld2"],
    "isotropy"               => ["misc/ISOTROPY/CIR_data.txt",
                                 "misc/ISOTROPY/PIR_data.txt"],
    )

"""
    paths(name) --> Vector{Pair{String,String}}

The files of data set `name`, as `source => destination` pairs: the source relative to
`data/`, the destination relative to the root of the artifact.
"""
paths(name::AbstractString) = [p isa Pair ? p : (p => basename(p)) for p in
    get(DATASETS, name) do
        error("unknown data set $(repr(name)); known: $(join(sort(collect(keys(DATASETS))), ", "))")
    end]

"""
    stage(name; from = DATA_DIR) --> String

Assemble the files of data set `name` in a fresh temporary directory, laid out as they will
be inside the artifact, and return that directory.

Files are taken from `from`, a `data/`-like tree — by default the local `data/` directory,
which is where a freshly built data set lands. A data set that no longer lives in the
repository (because it is published only as an artifact) must instead be staged from the
currently published artifact; `stage_from_artifact` does that.

Staging happens outside the repository on purpose. A git tree hash records the executable
bit, and a Windows drive mounted under WSL reports every file as executable and silently
ignores `chmod` — so staging there would give a different tree hash for identical data than
staging on a POSIX filesystem would.
"""
stage(name::AbstractString; from::AbstractString = DATA_DIR) =
    _stage(name, from, first)

"""
    stage_from_artifact(name) --> String

As [`stage`](@ref), but taking the files from the artifact that `Artifacts.toml` currently
points at, downloading it if necessary. Use this to amend a data set that is no longer kept
in the repository. Its files are already laid out as the artifact wants them, so they are
looked up by their destination names rather than their `data/` ones.
"""
function stage_from_artifact(name::AbstractString)
    @eval using Pkg.Artifacts: ensure_artifact_installed
    toml = joinpath(dirname(@__DIR__), "Artifacts.toml")
    from = Base.invokelatest(ensure_artifact_installed, name, toml)
    return _stage(name, from, last)
end

function _stage(name::AbstractString, from::AbstractString, srcof::Function)
    dir = mktempdir()
    for pair in paths(name)
        src = joinpath(from, srcof(pair))
        isfile(src) || error("$src does not exist; stage from the published artifact instead")
        dst = joinpath(dir, last(pair))
        mkpath(dirname(dst))
        cp(src, dst)
        chmod(dst, 0o644)
        filemode(dst) & 0o777 == 0o644 ||
            error("could not normalize permissions of $dst (got $(string(filemode(dst) & 0o777, base=8, pad=3))); \
                   stage on a POSIX filesystem, or the tree hash will not be reproducible")
    end
    return dir
end

"""
    package(name, dir) --> (tarball, tree_hash, sha256sum)

Compress the data set staged in `dir` into `build/data-release/<name>.tar.gz` and return
the two hashes that its `Artifacts.toml` entry needs: the git tree hash of the unpacked
tree, which Pkg verifies after unpacking, and the SHA-256 of the tarball, which it verifies
on download.

The tarball is reproducible: `Tar.create` normalizes timestamps, ownership and permissions,
and `gzip -n` omits its own timestamp, so identical input gives a byte-identical archive.
"""
function package(name::AbstractString, dir::AbstractString)
    mkpath(STAGE_DIR)
    tarball = joinpath(STAGE_DIR, name * ".tar.gz")
    mktemp() do tarpath, io
        close(io)
        Tar.create(dir, tarpath)
        # `run` waits for gzip to exit; piping into an `open(…, "w", io)` process and
        # closing it does not, and then the tarball may still be short when it is hashed
        run(pipeline(`gzip -9 -n -c $tarpath`, tarball))
        tree_hash = open(Tar.tree_hash, tarpath)
        return tarball, tree_hash, bytes2hex(open(sha256, tarball))
    end
end

"""
    artifacts_entry(name, tag, tree_hash, sha256sum) --> String

The `Artifacts.toml` stanza for data set `name` published under release `tag`.
"""
function artifacts_entry(name, tag, tree_hash, sha256sum)
    url = "https://github.com/$REPO/releases/download/$tag/$name.tar.gz"
    return """
    [$name]
    git-tree-sha1 = "$tree_hash"
    lazy = true

        [[$name.download]]
        url = "$url"
        sha256 = "$sha256sum"
    """
end

function main(tag, names = sort(collect(keys(DATASETS))))
    entries, tarballs = String[], String[]
    for name in names
        dir = isfile(joinpath(DATA_DIR, first(first(paths(name))))) ? stage(name) :
                                                                       stage_from_artifact(name)
        tarball, tree_hash, sha256sum = package(name, dir)
        push!(tarballs, tarball)
        push!(entries, artifacts_entry(name, tag, tree_hash, sha256sum))
        println("packaged $name: $(length(paths(name))) file(s), ",
                "$(round(filesize(tarball)/2^20; digits=2)) MiB → $(basename(tarball))")
    end
    println("\n--- Artifacts.toml entries ---\n")
    foreach(entry -> println(entry), entries)
    println("--- upload ---\n")
    println("gh release create $tag --title ... --notes ...   # if it does not exist yet")
    println("gh release upload $tag \\\n    ", join(tarballs, " \\\n    "), "\n")
    return nothing
end

if abspath(PROGRAM_FILE) == @__FILE__
    isempty(ARGS) && error("usage: julia --project=build build/data_release.jl <tag> [<name>...]")
    main(ARGS[1], length(ARGS) > 1 ? ARGS[2:end] : sort(collect(keys(DATASETS))))
end
