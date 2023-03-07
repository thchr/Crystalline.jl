struct LabeledKVec{D}
    label :: Symbol
    kv    :: KVec{D}
end
struct Connection{D}
    # assumed sorted such that H ⊂ G for the little group G and H (i.e., kᴳ is a maximal
    # point typically and kᴴ a nonmaximal manifold)
    kᴳ :: LabeledKVec{D}
    kᴴ :: LabeledKVec{D}
end
struct SubductionTable{D}
    num :: Int
    c :: Connection{D}
    irlabsᴳ :: Vector{String}
    irlabsᴴ :: Vector{String}
    table :: Matrix{Int}
    monodromy :: Bool
end
Base.show(io::IO, lk::LabeledKVec) = print(io, lk.label, "=", lk.kv)
function Base.show(io::IO, ::MIME"text/plain", c::Connection)
    print(io, c.kᴳ.label, " ↓ ", c.kᴴ.label)
    print(io, ": ", c.kᴳ.kv, " ↓ ", c.kᴴ.kv)
end
Base.show(io::IO, c::Connection) = print(io, c.kᴳ.label, " ↓ ", c.kᴴ.label)
function Base.show(io::IO,  t::SubductionTable{D}) where D
    summary(io, t)
    println(io, " ⋕", t.num, " (", iuc(t.num, D), ") for ", t.c, " connection:")
    pretty_table(io,
        t.table,
        header = t.irlabsᴴ,
        row_labels = t.irlabsᴳ,
        vlines = [1,],
        hlines = [:begin, 1, :end]
    )
end