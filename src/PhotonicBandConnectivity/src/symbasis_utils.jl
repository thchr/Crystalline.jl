# we bother to optimize this, as it can be a bottleneck; much faster than a naive 
# implementation like `sum(sb[idxs])`
function sum_symbases!(n, sb::SymBasis, idxs)
    Nⁱʳʳ = length(n)
    n .= sb[first(idxs)] # peel off 1st iter & ensure invariance to n's inititialization
    @inbounds for idx in @view idxs[2:end]
        nᴴ = sb[idx]
        for i in 1:Nⁱʳʳ
            n[i] += nᴴ[i]
        end
    end
    return n
end
sum_symbases(sb::SymBasis, idxs) = sum_symbases!(similar(first(sb)), sb, idxs)

function coef2idxs(c::AbstractVector{Int})
    N = sum(c)
    cⁱ = Vector{Int}(undef, N)
    pos₁, pos₂, idx = 0, 0, 0
    while true
        idx  = findnext(≠(0), c, idx+1)
        pos₁ = pos₂+1
        pos₂ = pos₂+c[idx]
        cⁱ[pos₁:pos₂] .= idx
        pos₂ == N && break
    end
    return cⁱ
end

function idxs2coef(cⁱ, Nᴴ) # Nᴴ = length(sb)
    c = zeros(Int, Nᴴ)
    for i in cⁱ
        c[i] += 1
    end
    return c
end
