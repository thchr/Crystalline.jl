# --- MSymOperation ---
"""
$(TYPEDEF)$(TYPEDFIELDS)
"""
struct MSymOperation{D} <: AbstractOperation{D}
    op :: SymOperation{D}
    tr :: Bool
end
SymOperation{D}(mop::MSymOperation{D}) where D = mop.op
timereversal(mop::MSymOperation) = mop.tr
seitz(mop::MSymOperation) = (s = seitz(mop.op); (mop.tr ? s : s * "′"))
xyzt(mop::MSymOperation) = xyzt(mop.op) # does not feature time-reversal!

function compose(
        mop₁ :: MSymOperation{D},
        mop₂ :: MSymOperation{D},
        modτ :: Bool=true) where D
    return MSymOperation{D}(compose(mop₁.op, mop₂.op, modτ), xor(mop₁.tr, mop₂.tr))
end

Base.:*(mop₁::MSymOperation{D}, mop₂::MSymOperation{D}) where D = compose(mop₁, mop₂)

function Base.isapprox(
        mop₁ :: MSymOperation{D},
        mop₂ :: MSymOperation{D}, 
        vs...; 
        kws...) where D
    return mop₁.tr == mop₂.tr && isapprox(mop₁.op, mop₂.op, vs...; kws...)
end

# --- MSpaceGroup ---
# Notation: there two "main" notations and settings:
#   - BNS (Belov, Neronova, & Smirnova)
#   - OG  (Opechowski & Guccione)
# We use the BNS setting (and notation), not OG (see Sec. 4, Acta Cryst. A78, 99 (2022)). 
# The `num` field of `MSpaceGroup` gives the the two BNS numbers:
#   - `num[1]`: F space-group associated number
#   - `num[2]`: crystal-system sequential number
"""
$(TYPEDEF)$(TYPEDFIELDS)
"""
struct MSpaceGroup{D} <: AbstractGroup{D, MSymOperation{D}}
    num :: Tuple{Int, Int}
    operations :: Vector{MSymOperation{D}}
end
label(msg::MSpaceGroup) = MSG_BNS_LABELs_D[msg.num]::String