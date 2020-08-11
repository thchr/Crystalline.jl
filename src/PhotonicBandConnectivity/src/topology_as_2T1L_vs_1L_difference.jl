"""
    topology_from_2T1L_xor_1L(nᵀ⁺ᴸ::Vector{Int}, nᴸ::Vector{Int}, nontopo_M::Matrix{Int}

Determines a Z₂-type topology a transverse symmetry vector (``T``), defined as the
difference of a transverse + longitudinal (``T+L``) symmetry vector `nᵀ⁺ᴸ` and a 
longitudinal (``L``) symmetry vector `nᴸ`.

Input: 
- `nᵀ⁺ᴸ`: can be computed from an index vector `cⁱ` and a compatibility Hilbert basis 
  `sb::SymBasis` by `sum_symbases(sb, cⁱ)`; 
- `nᴸ = sb[idx¹ᴸ]`: where `idx¹ᴸ` is the 1L-constraints pick in the `sb` basis;
- `nontopo_M`: a matrix representation of the nontopological Hilbert basis `nontopo_sb`,
  obtained from `matrix(nontopo_sb)`.

Output:
- a member of the enum `SymBases.TopologyKind`, either `trivial=0` or `nontrivial=1`.

The idea is quite simple: infer the topology of ``T`` solutions from the "difference" of the
topologies of the ``T+L`` and ``L`` solutions in the xor (`⊻`) sense, i.e. ``T`` is 
nontrivial if the topologies of ``T+L`` and ``L`` differ:

    ┌────────────┬────────────╥────────────╖        ┌───┬───╥───────╖
    │      L     │     T+L    ║     T Z₂   ║        │ x │ y ║ x ⊻ y ║
    ├────────────┼────────────╫────────────╢        ├───┼───╫───────╢
    |  trivial   │  trivial   ║  trivial   ║        | 0 │ 0 ║   0   ║
    |  trivial   │ nontrivial ║ nontrivial ║        | 0 │ 1 ║   1   ║
    | nontrivial │  trivial   ║ nontrivial ║        | 1 │ 0 ║   1   ║
    | nontrivial │ nontrivial ║  trivial   ║        | 1 │ 1 ║   0   ║
    └────────────┴────────────╨────────────╜        └───┴───╨───────╜

This is a meaningful definition, because while neither the ``T+L`` or ``L`` symmetry vectors
in general are unique, their difference is (i.e. ``T``, up to singular ``Γ``-irrep content):
same goes for the associated topology.
"""
function topology_from_2T1L_xor_1L(nᵀ⁺ᴸ::Vector{Int}, nᴸ::Vector{Int}, nontopo_M::Matrix{Int})

    # check topology of L symmetry vector
    is_nontrivialᴸ = nᴸ ∉ eachcol(nontopo_M)

    # check topology of T+L symmetry vector (using `get_solution_topology` from SymBases)
    topologyᵀ⁺ᴸ = get_solution_topology(nᵀ⁺ᴸ, nontopo_M, nothing, nothing) # returns Enum `trivial` or `nontrivial`
    is_nontrivialᵀ⁺ᴸ = topologyᵀ⁺ᴸ == nontrivial

    # infer 2T solution's topology from xor-difference:
    is_nontrivialᵀ = is_nontrivialᴸ ⊻ is_nontrivialᵀ⁺ᴸ

    return is_nontrivialᵀ ? nontrivial : trivial
end


"""
    topology_from_2T1L_xor_1L(nᵀ::Vector{Int}, nᴸ::Vector{Int}, m²ᵀ::Vector{Int}, 
                              Γidxs::AbstractVector{<:Integer}, nontopo_M::Matrix{Int})

Same functionality as the equivalently named 3-variable signature method, but takes the 
transverse (physical; i.e. sans singular Γ-irrep content) symmetry vector `nᵀ`, the
pinned 2T symmetry content `m²ᵀ`, and the Γ indexing vector `Γidx` instead of `nᵀ⁺ᴸ` from 
which it reconstruct the latter. This is often a more convenient access point.
"""
function topology_from_2T1L_xor_1L(nᵀ::Vector{Int}, nᴸ::Vector{Int}, m²ᵀ::Vector{Int}, 
                                   Γidxs::AbstractVector{<:Integer}, nontopo_M::Matrix{Int})

    nᵀ⁺ᴸ = copy(nᵀ)      # reconstruct T+L solution from nᵀ, nᴸ, and m²ᵀ (note the tricky
    nᵀ⁺ᴸ[Γidxs] .+= m²ᵀ  # indexing with Γidxs)
    nᵀ⁺ᴸ .+= nᴸ

    return topology_from_2T1L_xor_1L(nᵀ⁺ᴸ, nᴸ, nontopo_M)
end
