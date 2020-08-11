using Crystalline
using LightGraphs, MetaGraphs
using LinearAlgebra: qr

#=
    Many of the methods in this file are generalizations of similar very similar methods in
    Crystalline/src/compatibility.jl.
    The purpose of this file is to allow methods to give an indication of the connectedness
    between high-symmetry points in the Brillouin zone, and how they are connected. 
    
    TODO: These methods should probably eventually be merged with the functionality in 
          Crystalline/src/compatibility.jl. At present, it is just easier to have it 
          separated, because the use-case here is quite specific to analyzing the validity
          of the symmetry vectors extracted by PhotonicBandConnectivity
=#

function is_compatible_kvec(kv::KVec, kv′::KVec, cntr::Char)
    # This method determines whether there is a solution to 
    #   kv + G = kv′(αβγ′)
    # and if so, returns G and αβγ′. kv must be special and kv′ should be nonspecial
    isspecial(kv) || throw(DomainError(kv, "must be special"))
    isspecial(kv′) && return false, nothing, nothing

    k₀, _      = parts(primitivize(kv,  cntr)) 
    k₀′, kabc′ = parts(primitivize(kv′, cntr))

    # least squares solve via QR factorization; equivalent to pinv(kabc)*(k₀-k₀′) but faster
    QR_kabc′ = qr(kabc′, Val(true))
    k₀G = similar(k₀)
    Gspan = (0,-1,1)
    for Gx in Gspan
        for Gy in Gspan
            for Gz in Gspan
                k₀G .= (k₀[1]+Gx, k₀[2]+Gy, k₀[3]+Gz)
                αβγ′ = QR_kabc′\(k₀G-k₀′)
                k′ = k₀′ + kabc′*αβγ′
                # check if least squares solution actually is a solution
                compat_bool = isapprox(k₀G, k′, atol=Crystalline.DEFAULT_ATOL)
                if compat_bool
                    return compat_bool, αβγ′, [Gx,Gy,Gz]
                end
            end
        end
    end

    return false, nothing, nothing
end

function find_compatible_kvec(kv::KVec, kvs′::Vector{KVec}, cntr::Char)
    !isspecial(kv) && throw(DomainError(kv, "input kv must be a special k-point"))

    compat_idxs = Vector{Int64}()
    compat_αβγs = Vector{Vector{Float64}}()
    compat_Gs   = Vector{Vector{Int}}()
    @inbounds for (idx′, kv′) in enumerate(kvs′)
        isspecial(kv′) && continue # must be a line/plane/general point to match a special point kv
        compat_bool, αβγ′,G = is_compatible_kvec(kv, kv′, cntr)
        if compat_bool
            push!(compat_idxs, idx′)
            push!(compat_αβγs, αβγ′)
            push!(compat_Gs, G)
        end
    end

    return compat_idxs, compat_αβγs, compat_Gs
end

# let kvsᴬ refer to a set of special k-points and kvsᴮ to a set of mixed k-points (special
# and non-special)
function connectivity((kvsᴬ, klabsᴬ), (kvsᴮ, klabsᴮ), cntr)

    all(isspecial, kvsᴬ) || throw(DomainError(kvsᴬ, "must only include special points"))

    Nkᴬ = length(kvsᴬ)
    D = dim(first(kvsᴬ))

    # find compatible vectors between A and B
    compat_idxs = Vector{Vector{Int}}(undef, Nkᴬ)
    compat_αβγs = Vector{Vector{Vector{Float64}}}(undef, Nkᴬ)
    compat_Gs = Vector{Vector{Vector{Int}}}(undef, Nkᴬ)
    cgraph = MetaGraph(Nkᴬ) # connectivity graph for special k-vecs
    for (idxᴬ, (kvᴬ, klabᴬ)) in enumerate(zip(kvsᴬ, klabsᴬ))
        compat_idxs[idxᴬ], compat_αβγs[idxᴬ], compat_Gs[idxᴬ] = find_compatible_kvec(kvᴬ, kvsᴮ, cntr)
        set_props!(cgraph, idxᴬ, Dict(:kvec=>kvᴬ, :klab=>klabᴬ))  # TODO: Add idx of kvᴬ in kvsᴮ?
    end

    for (idxᴬ¹, kvᴬ¹) in enumerate(kvsᴬ)
        compat_idxs¹ = compat_idxs[idxᴬ¹]
        for (idxᴬ², kvᴬ²) in enumerate(kvsᴬ)
            idxᴬ² ≤ idxᴬ¹ && continue # avoid double & self-comparisons
            compat_idxs² = compat_idxs[idxᴬ²]

            # find overlap of compat_idxs¹ and compat_idxs² which are lines
            local_idxsᴮ = findall(compat_idxs²) do idxᴮ
                (idxᴮ∈compat_idxs¹) & (Crystalline.nfreeparams(kvsᴮ[idxᴮ]) == 1)
            end
            isempty(local_idxsᴮ) && continue
            idxsᴮ = getindex.(Ref(compat_idxs²), local_idxsᴮ)

            for idxᴮ in idxsᴮ
                println(
                    klabsᴬ[idxᴬ¹], " = ", string(kvsᴬ[idxᴬ¹]), " connected to ",
                    klabsᴬ[idxᴬ²], " = ", string(kvsᴬ[idxᴬ²]), " via ",
                    klabsᴮ[idxᴮ], " = ", string(kvsᴮ[idxᴮ]))
            end
            add_edge!(cgraph, idxᴬ¹, idxᴬ²)
            set_props!(cgraph, Edge(idxᴬ¹, idxᴬ²), 
                Dict(:klabs=>getindex.(Ref(klabsᴮ), idxsᴮ), 
                     :kvecs=>getindex.(Ref(kvsᴮ), idxsᴮ), 
                     :kidxs=>idxsᴮ)
               )
        end
    end
   
    return cgraph
end

## test
sgnum = 230;
sb = compatibility_bases(sgnum)[1];
lgirsvec=realify.(get_lgirreps(sgnum));
kvsᴬ, klabsᴬ = sb.kvs, sb.klabs;
kvsᴮ, klabsᴮ = kvec.(first.(lgirsvec)), klabel.(first.(lgirsvec));
cntr = centering(sgnum);
cg = connectivity((kvsᴬ, klabsᴬ), (kvsᴮ, klabsᴮ), cntr);