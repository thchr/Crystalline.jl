#module ParseMSG

using Crystalline
# --- definitions ---
# to be transferred to Crystalline's definitions
using Crystalline: AbstractGroup
using StaticArrays

struct MSymOperation{D}
    op :: SymOperation{D}
    tr :: Bool
end
Base.show(io::IO, mop::MSymOperation) = (show(io, mop.op); print(io, mop.tr ? "′" : ""))
function Crystalline.compose(
        mop₁::MSymOperation{D},
        mop₂::MSymOperation{D},
        modτ::Bool=true) where D
    return MSymOperation{D}(compose(mop₁.op, mop₂.op, modτ), xor(mop₁.tr, mop₂.tr))
end
Base.:*(mop₁::MSymOperation{D}, mop₂::MSymOperation{D}) where D = compose(mop₁, mop₂)
function Base.isapprox(mop₁::MSymOperation{D}, mop₂::MSymOperation{D}, vs...; kws...) where D
    isapprox(mop₁.op, mop₂.op, vs...; kws...) && mop₁.tr == mop₂.tr
end
struct MSpaceGroup{D} <: AbstractGroup{D}
    # Note: we use BNS numbers, not  OG numbers (see Sec. 4, Acta Cryst. A78, 99 (2022))
    # BNS number: `num[1]`: F space-group associated number
    #             `num[2]`: crystal-system sequential number
    num :: Tuple{Int, Int}
    operations :: Vector{MSymOperation{D}}
end
Crystalline.label(msg::MSpaceGroup) = MSG_BNS_LABELs_D[msg.num]::String
# --- parsing ---
#cd(@__DIR__)
io = open((@__DIR__)*"/../data/operations/msgs/iso-magnetic_table_bns.txt")

# --- read dictionary of operations ---
function read_ops!(io, ops_d, stop_string)
    while (str=readline(io); str ≠ stop_string)
        parts = split(str, isspace; keepempty=false)
        code, xyz = parts[1], parts[2]
        op = SymOperation(xyz)
        ops_d[code] = op
    end
    return ops_d
end

readuntil(io, "Non-hexagonal groups:"); readline(io)
nonhex_ops_d = Dict{String, SymOperation{3}}()
read_ops!(io, nonhex_ops_d, "Hexagonal groups:")
hex_ops_d    = Dict{String, SymOperation{3}}()
read_ops!(io, hex_ops_d, "----------")

# --- read individual groups ---
function parse_opstr(str::AbstractString, ops_d)
    str′ = strip(str,  ['(', ')','\''])
    code, translation_str = split(str′, '|')::Vector{SubString{String}}
    translation_parts = split(translation_str, ',')
    translation_tup = ntuple(Val(3)) do i 
        Crystalline.parsefraction(translation_parts[i])
    end
    translation = SVector{3, Float64}(translation_tup)
    op = SymOperation{3}(ops_d[code].rotation, translation)
    tr = last(str) == '\''
    return MSymOperation(op, tr)
end

function read_group!(io)
    # BNS: number & label
    readuntil(io, "BNS: ")
    num_str = readuntil(io, " ") # numbering in format `sgnum.cnum`
    sgnum, cnum = parse.(Int, split(num_str, "."))::Vector{Int}
    bns_label = replace(readuntil(io, " "), '\''=>'′') # BNS label in IUC-style

    # OG: number and label
    readuntil(io, "OG: ")
    N₁N₂N₃_str = readuntil(io, " ") # numbering in format `N₁.N₂.N₃`
    N₁, N₂, N₃ = parse.(Int, split(N₁N₂N₃_str, "."))::Vector{Int}
    og_label   = replace(readuntil(io, '\n'), '\''=>'′') # OG label in IUC-style

    # operator strings
    is_hex = crystalsystem(sgnum) ∈ ("hexagonal" , "trigonal")
    readuntil(io, "Operators")
    c = read(io, Char)
    if c == ' '
        readuntil(io, "(BNS):")
    elseif c ≠ ':'
        error("unexpected parsing failure")
    end
    ops_str = readuntil(io, "Wyckoff")
    op_strs = split(ops_str, isspace; keepempty=false)

    # parse operator strings to operations
    g = Vector{MSymOperation{3}}(undef, length(op_strs))
    for (i, str) in enumerate(op_strs)
        g[i] = parse_opstr(str, is_hex ? hex_ops_d : nonhex_ops_d)
        isdone = str[end] == '\n'
        isdone && break
    end

    return MSpaceGroup{3}((sgnum, cnum), g), bns_label, og_label, (N₁, N₂, N₃)
end
msgs = MSpaceGroup{3}[]
MSG_BNS_LABELs_D = Dict{Tuple{Int, Int}, String}()
MSG_OG_LABELs_D = Dict{Tuple{Int, Int, Int}, String}()
MSG_BNS2OG_NUMs_D = Dict{Tuple{Int, Int}, Tuple{Int, Int, Int}}()
for i in 1:1651
    msg, bns_label, og_label, (N₁, N₂, N₃) = read_group!(io)
    push!(msgs, msg)
    MSG_BNS_LABELs_D[msg.num] = bns_label
    MSG_OG_LABELs_D[(N₁, N₂, N₃)] = og_label
    MSG_BNS2OG_NUMs_D[msg.num] = (N₁, N₂, N₃)
end

#end # module ParseMSG