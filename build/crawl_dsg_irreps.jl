# Crawl the irreps of the *double* space groups (and double point groups) from the Bilbao
# Crystallographic Server, i.e. the data underlying REPRESENTATIONS DSG & DPG. This is the
# source for Crystalline's spinful (double-valued) little group irreps; ISOTROPY, which
# supplies the spinless irreps, has no spinful data.
#
# This script only fetches & caches raw HTML; parsing lives in the companion
# `parse_dsg_irreps.jl` (deliberately separate: parsing is re-run often, the crawl only when
# the cache is incomplete). The cache is gitignored; only parsed output is ever committed.
#
# ## Request flow (space groups)
#
# Bilbao's DSG program is a three-step form; step 1 is just a group picker, so we start at 2:
#   2. POST `representations_vec.pl`  {tipogrupo=dbg, super=<sgnum>}
#        → an HTML page of radio buttons, one per k-vector, with values "<klab>&(<kv>)";
#          this *is* the authoritative k-vector list for the group.
#   3. POST `representations_out.pl`  {tipogrupo=dbg, super=<sgnum>, texto=Double,
#                                      phys=, vecfinal=<klab>&(<kv>)}
#        → matrices of the little group irreps, single- *and* double-valued, plus the SU(2)
#          spin lift and the Seitz symbol of every operation.
# The single-valued irreps come along for free; `write_dsg_irreps.jl` stores them
# separately, and `validate_bilbao_vs_isotropy.jl` checks them against the ISOTROPY irreps
# of `lgirreps`, before we trust any of the double-valued data.
#
# ## Request flow (point groups)
#
# A single GET, identical in shape to `crawl_and_write_pgs.jl`'s `BILBAO_PG_URL_BASE` but
# with `tipogrupo=dbg` in place of `spg`; we reuse that file's parent-space-group lookup.
#
# ## Anti-bot gate
#
# BCS gates its CGI programs behind a Cloudflare Turnstile challenge; see the header of
# `../BandGraphs/build/crawl_mbandpaths.jl` for the full account. In short: without a
# clearance cookie the `*.pl` programs return a bare `500 Internal Server Error`, which is
# the gate and not downtime. The widget does not auto-complete headless, so a one-time
# interactive solve in a real browser is required.
#
# To obtain the cookie: open
#   https://cryst.ehu.es/cgi-bin/cryst/programs/representations.pl?tipogrupo=dbg
# in a browser, pass the challenge, then copy the cookie for cryst.ehu.es from DevTools →
# Application → Cookies, as "turnstile_passed=<value>". The value is the *issue epoch*, so
# its freshness can be checked arithmetically (`cookie_age`, below); lifetime is ≈ 1 h.
#
# ## Usage
#
#   julia --project=build build/crawl_dsg_irreps.jl probe  [--cookie "<c>"]
#   julia --project=build build/crawl_dsg_irreps.jl crawl  --cookie "<c>" \
#       [--sgnums 1:230] [--what vecs|irreps|pgs|all] [--delay 2]
#   julia --project=build build/crawl_dsg_irreps.jl status
#
# `crawl` is resumable: cached files are skipped, so re-running after the cookie expires
# picks up where it stopped. A full sweep is ≈ 230 + 4200 requests; at the default delay
# that is several hours, i.e. several cookies. Do not lower `--delay` below ~2 s — BCS has
# temporarily blocked an IP for hammering before.

using HTTP
using Dates
# NB: Crystalline itself is loaded *lazily*, only for the point-group crawl (`crawl_pgs`),
#     which needs its parent-point-group machinery. The space-group crawl, `probe` and
#     `status` deliberately need nothing from it: a crawler that cannot run because the
#     package it feeds fails to load would be a poor tool, and `using Crystalline` opens the
#     JLD2 irrep files at `__init__` for no benefit here.

const BCS_HOST     = "https://cryst.ehu.es"
const VEC_URL      = BCS_HOST * "/cgi-bin/cryst/programs/representations_vec.pl"
const OUT_URL      = BCS_HOST * "/cgi-bin/cryst/programs/representations_out.pl"
const ENTRY_URL    = BCS_HOST * "/cgi-bin/cryst/programs/representations.pl?tipogrupo=dbg"
const PG_URL_BASE  = BCS_HOST * "/cgi-bin/cryst/programs/representations_out.pl?" *
                     "tipogrupo=dbg&pointspace=point&"
const USER_AGENT   = "Mozilla/5.0 (Windows NT 10.0; Win64; x64) AppleWebKit/537.36 " *
                     "(KHTML, like Gecko) Chrome/126.0.0.0 Safari/537.36"

const CRAWL_DIR = joinpath(@__DIR__, "crawls", "dsg")
vec_path(sgnum)        = joinpath(CRAWL_DIR, "vec", "sg$(sgnum).html")
out_path(sgnum, klab)  = joinpath(CRAWL_DIR, "out", "sg$(sgnum)-$(klab).html")
pg_path(pgiuc)         = joinpath(CRAWL_DIR, "pg", replace(pgiuc, '/'=>"_") * ".html")
sel_path(sgnum, key)   = joinpath(CRAWL_DIR, "sel", "sg$(sgnum)-$(key).html")

# ---------------------------------------------------------------------------------------- #
# Fetch layer

struct GateBlocked <: Exception
    status :: Int
    body   :: String
end
Base.showerror(io::IO, e::GateBlocked) = print(io,
    "BCS anti-bot gate: HTTP $(e.status), body ≈ $(repr(first(e.body, 80))). ",
    "Pass a fresh clearance cookie (see the header of this script).")

"""
    cookie_age(cookie) --> Union{Int, Nothing}

Seconds since the `turnstile_passed` cookie was issued, or `nothing` if the cookie carries no
such field. The cookie's value is its issue epoch, so a stale cookie can be detected *before*
spending a request on it — worth doing, since a stale cookie fails identically to no cookie.
"""
function cookie_age(cookie::AbstractString)
    m = match(r"turnstile_passed=(\d+)", cookie)
    m === nothing && return nothing
    return Int(round(datetime2unix(now(UTC)))) - parse(Int, m.captures[1])
end

function default_headers(cookie::AbstractString, referer::AbstractString=ENTRY_URL)
    return ["User-Agent" => USER_AGENT,
            "Accept" => "text/html,application/xhtml+xml,application/xml;q=0.9,*/*;q=0.8",
            "Accept-Language" => "en-US,en;q=0.9",
            "Referer" => referer,
            "Cookie" => cookie]
end

"""
    fetch(method, url, cookie; body, referer, tries) --> String

One request, retrying *transient* failures with a linear backoff. A multi-hour sweep will
meet the occasional dropped connection, and dying on it wastes the rest of a cookie for no
reason.

`GateBlocked` is deliberately **not** retried: that is a cookie problem, not a network one, and
retrying it would hammer the server precisely when it is refusing us.
"""
function fetch(method::AbstractString, url::AbstractString, cookie::AbstractString;
               body=nothing, referer::AbstractString=ENTRY_URL, tries::Int=4)
    for attempt in 1:tries
        try
            r = HTTP.request(method, url, default_headers(cookie, referer), something(body, "");
                             status_exception=false, redirect=false)
            s = String(r.body)
            (r.status ≠ 200 || occursin("turnstile", lowercase(s))) && throw(GateBlocked(r.status, s))
            return s
        catch e
            e isa GateBlocked && rethrow()      # a cookie problem; retrying would only hammer
            attempt == tries && rethrow()
            @warn "transient fetch failure, retrying" attempt exception=e
            sleep(5*attempt)
        end
    end
    error("unreachable")
end

# `vecfinal` values contain a literal '&' ("X&(1/2,0,0)"), so every field must be escaped
formencode(kvs::Pair...) = join((HTTP.escapeuri(k) * "=" * HTTP.escapeuri(v) for (k,v) in kvs), '&')

# ---------------------------------------------------------------------------------------- #
# k-vector listing

"""
    kvec_specifiers(vec_html) --> Vector{String}

The `vecfinal` form values of a `representations_vec.pl` page, e.g. `"X&(1/2,0,0)"`. This is
Bilbao's own k-vector list for the group and is what we iterate over; we deliberately do not
reconstruct it from Crystalline, so that a coverage difference shows up as data rather than
being silently imposed.
"""
function kvec_specifiers(vec_html::AbstractString)
    return [m.captures[1] for m in eachmatch(r"name=\"vecfinal\"\s+value=\"([^\"]+)\"", vec_html)]
end
# ⚠ Bilbao sometimes *bundles* several k-manifolds behind a single form entry, e.g. sg 16's
#   "GP,W,V,N,M,L,K&(u,v,w)s(u,v,1/2)s…", or sg 1's "GP,GM,Z,Y,X,V,U,T,R&…" — seemingly where
#   the manifolds share a little-group structure. Submitting such an entry does **not** return
#   irreps: it returns a second selection page, whose radios carry the *bare* labels ("GP",
#   "W", …), and one must submit again to reach the irreps. So a bundle costs 1 + n requests
#   and yields n pages. `cachekey_of` names the intermediate page; `klabels_of` the manifolds
#   behind it. Every cached `out/` page is exactly one k-manifold.
cachekey_of(spec::AbstractString) = String(first(split(spec, '&')))
klabels_of(spec::AbstractString) = String.(split(cachekey_of(spec), ','))
isbundle(spec::AbstractString) = length(klabels_of(spec)) > 1

# ---------------------------------------------------------------------------------------- #
# Crawl steps

function crawl_vecs(sgnums, cookie; delay::Real=2)
    mkpath(joinpath(CRAWL_DIR, "vec"))
    for sgnum in sgnums
        path = vec_path(sgnum)
        isfile(path) && continue
        html = fetch("POST", VEC_URL, cookie;
                     body = formencode("tipogrupo"=>"dbg", "super"=>string(sgnum),
                                       "list"=>"Submit"))
        write(path, html)
        println("vec: sg $sgnum → $(length(kvec_specifiers(html))) k-vectors")
        sleep(delay)
    end
end

function fetch_irreps(sgnum, spec, cookie)
    return fetch("POST", OUT_URL, cookie;
                 body = formencode("tipogrupo"=>"dbg", "super"=>string(sgnum),
                                   "texto"=>"Double", "symbol"=>"", "phys"=>"",
                                   "vecfinal"=>spec, "list"=>"Submit"),
                 referer = VEC_URL)
end

function crawl_irreps(sgnums, cookie; delay::Real=2)
    mkpath(joinpath(CRAWL_DIR, "out"))
    mkpath(joinpath(CRAWL_DIR, "sel"))
    for sgnum in sgnums
        isfile(vec_path(sgnum)) || (@warn "no vec page for sg $sgnum; run `--what vecs` first"; continue)
        for spec in kvec_specifiers(read(vec_path(sgnum), String))
            klabs = klabels_of(spec)
            if !isbundle(spec)
                path = out_path(sgnum, only(klabs))
                isfile(path) && continue
                write(path, fetch_irreps(sgnum, spec, cookie))
                println("out: sg $sgnum, k = $(only(klabs))")
                sleep(delay)
            else
                # step 1: the intermediate selection page (cached so a resume can skip it)
                selpath = sel_path(sgnum, cachekey_of(spec))
                if !all(klab -> isfile(out_path(sgnum, klab)), klabs)
                    if !isfile(selpath)
                        write(selpath, fetch_irreps(sgnum, spec, cookie))
                        sleep(delay)
                    end
                end
                # step 2: one request per manifold, submitting the *bare* label
                for klab in klabs
                    path = out_path(sgnum, klab)
                    isfile(path) && continue
                    write(path, fetch_irreps(sgnum, klab, cookie))
                    println("out: sg $sgnum, k = $klab (from bundle $(cachekey_of(spec)))")
                    sleep(delay)
                end
            end
        end
    end
end

# NB: deliberately *not* `include`ing `crawl_and_write_pgs.jl` to reuse this — that file ends
#     in a top-level `__crawl_and_write_3d_pgirreps()` call, which would re-crawl and
#     overwrite `data/irreps/pgs/3d/irreps_data.jld2` as a side effect of the include.
#     The lookup is small enough to restate; keep it in sync with that file if it changes.
function findfirst_matching_parent_sgnum(pgiuc::String)
    for sgnum in 1:Crystalline.MAX_SGNUM[3]
        pg = Crystalline.find_parent_pointgroup(Crystalline.spacegroup(sgnum, Val(3)))
        Crystalline.label(pg) == pgiuc && return sgnum, Crystalline.num(pg)
    end
    throw(DomainError(pgiuc, "requested label cannot be found"))
end

# The Crystalline-dependent body is split out and reached through `invokelatest` below: the
# `import` happens while `crawl_pgs` is already running, so its methods are *newer* than the
# world age of the running frame and a direct call fails with "method too new to be called
# from this world context". `invokelatest` is what lets the freshly-imported methods be seen.
function crawl_pgs(cookie; delay::Real=2)
    mkpath(joinpath(CRAWL_DIR, "pg"))
    @eval import Crystalline # lazy; see the note by the imports
    return Base.invokelatest(_crawl_pgs, cookie, delay)
end

function _crawl_pgs(cookie, delay)
    for pgiuc in Crystalline.PG_IUCs[3]
        path = pg_path(pgiuc)
        isfile(path) && continue
        parent_sgnum, pgnum = findfirst_matching_parent_sgnum(pgiuc)
        url = PG_URL_BASE * "num=$(parent_sgnum)&super=$(pgnum)&symbol=$(HTTP.escapeuri(pgiuc))"
        write(path, fetch("GET", url, cookie))
        println("pg: $pgiuc")
        sleep(delay)
    end
end

# ---------------------------------------------------------------------------------------- #
# Status

function status(sgnums=1:230)
    nvec = count(sgnum -> isfile(vec_path(sgnum)), sgnums)
    println("vec pages: $nvec/$(length(sgnums))")
    total = have = 0
    missing_ks = Pair{Int, Vector{String}}[]
    for sgnum in sgnums
        isfile(vec_path(sgnum)) || continue
        specs = kvec_specifiers(read(vec_path(sgnum), String))
        klabs = reduce(vcat, klabels_of.(specs); init=String[])
        total += length(klabs)
        gone = filter(klab -> !isfile(out_path(sgnum, klab)), klabs)
        have += length(klabs) - length(gone)
        isempty(gone) || push!(missing_ks, sgnum => gone)
    end
    println("irrep pages: $have/$total")
    # 37, not 32: `PG_IUCs` distinguishes setting variants (312/321, 3m1/31m, -31m/-3m1, …)
    npg = isdir(joinpath(CRAWL_DIR, "pg")) ? length(readdir(joinpath(CRAWL_DIR, "pg"))) : 0
    println("point groups: $npg/37")
    if !isempty(missing_ks)
        println("incomplete space groups ($(length(missing_ks))):")
        for (sgnum, gone) in first(missing_ks, 20)
            println("  sg $sgnum: ", join(gone, ", "))
        end
        length(missing_ks) > 20 && println("  … and $(length(missing_ks)-20) more")
    end
    return nothing
end

# ---------------------------------------------------------------------------------------- #
# CLI

function getarg(args, flag, default=nothing)
    i = findfirst(==(flag), args)
    return i === nothing ? default : args[i+1]
end

function main(args=ARGS)
    cmd = isempty(args) ? "status" : args[1]
    cookie = something(getarg(args, "--cookie"), get(ENV, "BCS_COOKIE", ""))

    if cmd ≠ "status"
        isempty(cookie) && error("no cookie: pass --cookie or set BCS_COOKIE (see script header)")
        age = cookie_age(cookie)
        if age !== nothing
            println("cookie issued $(age)s ago", age > 3600 ? " ⚠ likely EXPIRED (lifetime ≈ 1 h)" : "")
        end
    end

    if cmd == "probe"
        try
            fetch("GET", ENTRY_URL, cookie)
            println("gate: OK — cookie accepted")
        catch e
            e isa GateBlocked || rethrow()
            showerror(stdout, e); println()
        end
    elseif cmd == "crawl"
        sgnums = eval(Meta.parse(something(getarg(args, "--sgnums"), "1:230")))
        what   = something(getarg(args, "--what"), "all")
        delay  = parse(Float64, something(getarg(args, "--delay"), "2"))
        try
            what ∈ ("vecs", "all")   && crawl_vecs(sgnums, cookie; delay)
            what ∈ ("irreps", "all") && crawl_irreps(sgnums, cookie; delay)
            what ∈ ("pgs", "all")    && crawl_pgs(cookie; delay)
        catch e
            e isa GateBlocked || rethrow()
            # the expected end of a session: the cookie has expired mid-sweep. Everything
            # fetched so far is on disk, so a rerun with a fresh cookie simply continues.
            println(); showerror(stdout, e); println("\n→ partial progress kept; rerun with a fresh cookie.")
            status(sgnums)
            exit(9)
        end
        status(sgnums)
    elseif cmd == "status"
        status()
    else
        error("unknown command $(repr(cmd)); expected probe, crawl, or status")
    end
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
