import BandGraphs: subduction_table

function dummy_lgirrep(irlab, irdim, special::Bool=true, D=3)
    klab = filter(isletter, irlab)
    kv = special ? KVec("0,0,0") : KVec("u,0,0")
    lg = LittleGroup{D}(1, kv, klab, [S"x,y,z"])
    return LGIrrep(
        irlab,
        lg,
        [ones(ComplexF64, irdim, irdim)],
        [zeros(Float64, D)],
        REAL, false)
end

function connection(lgirs1, lgirs2)
    k1 = BandGraphs.LabeledKVec(Symbol(klabel(first(lgirs1))), position(lgirs1))
    k2 = BandGraphs.LabeledKVec(Symbol(klabel(first(lgirs2))), position(lgirs2))
    return Connection(k1, k2)
end

function subduction_table(lgirs1, lgirs2, table, monodromy=false)
    return SubductionTable{3}(
            1, connection(lgirs1, lgirs2), label.(lgirs1), label.(lgirs2), table, monodromy)
end

lgirsd = Dict{String, Collection{LGIrrep{3}}}(
    "A" => Collection([dummy_lgirrep("A₁", 1), dummy_lgirrep("A₂", 3)]),
    "B" => Collection([dummy_lgirrep("B₁", 1), dummy_lgirrep("B₂", 2)]),
    "C" => Collection([dummy_lgirrep("C₁", 1), dummy_lgirrep("C₂", 1)]),
    "Δ" => Collection([dummy_lgirrep("Δ₁", 1, false), dummy_lgirrep("Δ₂", 2, false)]),
    "Λ" => Collection([dummy_lgirrep("Λ₁", 1, false), dummy_lgirrep("Λ₂", 1, false)]),
)
lgirsv = [lgirs for lgirs in values(lgirsd) if isspecial(lgirs[1])]
n_str = "[A₂, B₁+B₂, 2C₁+C₂]"
n = parse(SymmetryVector{3}, n_str, lgirsv)

subt_AΔ = subduction_table(lgirsd["A"], lgirsd["Δ"], [1 0; 1 1])
subt_BΔ = subduction_table(lgirsd["B"], lgirsd["Δ"], [1 0; 0 1])
subt_CΛ = subduction_table(lgirsd["C"], lgirsd["Λ"], [1 0; 0 1])
subt_AΛ = subduction_table(lgirsd["A"], lgirsd["Λ"], [1 0; 1 1])

subts = [subt_AΔ, subt_BΔ]

bandg = build_subgraphs(n, subts, lgirsd);

split_bandg = complete_split(bandg, ("A₂", 1); separable_degree=2);
## 
# TODO
faxp, _ = plot_flattened_bandgraph(split_bandg[1]); fax # FAILS due to assumption of at most two-edges-per nonmax node
@test !isnothing(split_bandg)

## Example with a 3-fold degenerate that is `separable_degree=2` separable
criterion = (lgir) -> isspecial(lgir) && irdim(lgir) == 3
separable_degree = 2
sep, _ = findall_separable_vertices(criterion, n, subts, lgirsd; separable_degree = 2)
@test is_separable(sep)