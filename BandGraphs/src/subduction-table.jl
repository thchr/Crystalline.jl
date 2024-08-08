function SubductionTable(
            c::Connection{D},
            sg::SpaceGroup{D},
            lgirsd::AbstractDict{String, <:AbstractVector{LGIrrep{D}}}
            ) where D

    kᴳ, kᴴ = c.kᴳ, c.kᴴ
    bare_kᴳ_label = rstrip(string(kᴳ.label), '′') # w/o monodromy indicator
    lgirsᴳ = lgirsd[bare_kᴳ_label]

    # determine the free parameters for which kᴳ and kᴴ are compatible/"intersect"
    # NB: `Gp` below is a primitive reciprocal lattice vector s.t. kᴴₚ(αβγ′) + Gp = kᴳₚ
    #     with kᴴₚ and kᴳₚ denoting the primitive settings of kᴴ and kᴳ, respectively.
    cntr = centering(sg)
    bool, αβγ, αβγ′, Gp = Crystalline.can_intersect(primitivize(kᴳ.kv, cntr), 
                                                    primitivize(kᴴ.kv, cntr))
    bool || error("failed to prove compatibility of $(kᴳ) and $(kᴴ)")
    iszero(αβγ) || error("implementation cannot handle nonzero free parameters for G-group")

    # `Gp` (primitive setting here) could in principle contain parts that are actually
    # along `kᴴ`, which would give a different value of αβγ′: we want the reciprocal vector
    # to have no parts along `kᴴ` (i.e., the reciprocal vector should lie in a line/plane
    # orthogonal to the plane/line spanned by `kᴴ`; we fix it up here.
    # [a motivating example is space group 122 (M ↓ Λ)]
    if !iszero(Gp)
        kᴴₚ = primitivize(kᴴ.kv, cntr)
        Q = pseudo_inverse_smith(Matrix(kᴴₚ.free)) * Gp
        if norm(Q) > Crystalline.DEFAULT_ATOL
            ΔGp = kᴴₚ.free * Q
            if all(v -> abs(v - round(v)) < Crystalline.DEFAULT_ATOL, ΔGp)
                # we can & should change `Gp` by a reciprocal lattice vector, if, in fact,
                # such a change would make `Gp` orthogonal to `kᴴₚ`
                Gp_tentative = Gp - ΔGp
                αβγ′_tentative = αβγ′ + Q
                if norm(kᴴₚ(αβγ′_tentative) ⋅ Gp_tentative) < Crystalline.DEFAULT_ATOL
                    Gp = Gp_tentative
                    αβγ′ = αβγ′_tentative
                else
                    # give up ¯\_(ツ)_/¯
                    # can e.g., happen in space group 220 for H ↓ Λ = [1,1,1] ↓ [-α, α, -α]
                    # but apparently does not change the subduction table, so giving up
                    # seems okay, although undeniably iffy - the whole approach here is iffy
                end
            else
                error("there be dragons here!")
                # there be dragons... the `Gp` vector is not actually orthogonal to `kᴴ`
                # ought to be verified that it cannot change the irrep regardless
            end
        end
        # NB: we currently don't use `Gp` for anything - the idea is we shouldn't need to
    end

    lgirsᴴ_tmp = lgirsd[string(kᴴ.label)]
    cosetsᴴ = cosets(sg, group(first(lgirsᴴ_tmp)))
    lgirsᴴ = Crystalline.remap_lgirreps_to_point_in_kstar(lgirsᴴ_tmp, kᴴ.kv, cosetsᴴ)
    
    # compute subduction table of G- and H-group irreps at their αβ′-intersection
    table = [subduction_count(lgirᴳ, lgirᴴ, αβγ′) for lgirᴳ in lgirsᴳ, lgirᴴ in lgirsᴴ]
    
    monodromy = endswith(string(kᴳ.label), '′') ? true : false
    irlabsᴳ, irlabsᴴ = label.(lgirsᴳ), label.(lgirsᴴ)
    if monodromy
        # append '′' to monodromy-affected irrep-labels (subduction tables are already
        # right, since αβγ′ captures the difference from the standard "in-cell" case)
        # our convention is to add `′` to every occurence of the k-label; e.g.,
        # mapping as: "X₁ → X′₁", "X₁X₂ → X′₁X′₂", and "X₁⁺X₂⁻ → X′₁⁺X′₂⁻"
        irlabsᴳ = map(irlabsᴳ) do irlab
            io = IOBuffer()
            start = firstindex(irlab)
            while (r = findnext(bare_kᴳ_label, irlab, start)) !== nothing
                write(io, irlab[start:r[end]], "′")
                start = nextind(irlab, r[end])
            end
            start == firstindex(irlab) && error("missing irrep label-change!")
            write(io, irlab[start:end])
            String(take!(io))
        end
    end

    return SubductionTable{D}(num(sg), c, irlabsᴳ, irlabsᴴ, table, monodromy)
end

function pseudo_inverse_smith(A)
    # compute the pseudo-inverse of a matrix `A` using the Smith normal form
    # (https://en.wikipedia.org/wiki/Pseudoinverse#Using_the_Smith_normal_form)
    F = Crystalline.smith(A)
    A⁺ = F.Tinv * pinv(diagm(F)) * F.Sinv
    return A⁺
end
