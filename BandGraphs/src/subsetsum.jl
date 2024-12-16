using JuMP, HiGHS

"""
    solve_subset_sum_variant(Sʲs, T, R; 
        verbose::Bool = false, 
        model::Model = DEFAULT_HiGHS_MODEL, 
        preemptive_model_reset::Bool = true)

Solves the following variant of the subset sum problem, using binary integer programming.

## Problem formulation
Given a target list `T` = (T₁, …, T_N), a set of multisets `Sʲs` = {S⁽ʲ⁾}_{j=1}^{J} with 
S⁽ʲ⁾ = {s⁽ʲ⁾₁, …, s⁽ʲ⁾_{M⁽ʲ⁾}}, as well as a multiset `R` = {r₁, …, r_L}, determine if there
exist disjoint subset selections {S⁽ʲ⁾ₙ} and {Rₙ} with S⁽ʲ⁾ₙ ⊂ S⁽ʲ⁾ and Rₙ ⊂ R such that:

    ∑_{s ∈ S⁽ʲ⁾ₙ} s = tₙ + ∑_{r ∈ Rₙ} r

for all j = 1, …, J and n = 1, …, N and with ⋃ₙ S⁽ʲ⁾ₙ = S⁽ʲ⁾} and ⋃ₙ Rₙ = R (and 
S⁽ʲ⁾_n ∩ S⁽ʲ⁾_{n′} = ∅ if n ≠ n′ and similarly for R_n ∩ R_{n′}).

## Output
- First return value (`::Bool`): whether the problem is feasible (i.e., solvable).
- The model (`::Model`), which contains the solution in binary variables `x` (into `Sʲs`)
  and `y` (into `R`).

## Example
```jl
julia> Sʲs = [[1, 3, 2, 1], [3, 2, 2], [1, 1, 2, 2, 1]];
julia> T = [3, 3];
julia> R = [1];
julia> solve_subset_sum_variant(Sʲs, T, R; verbose=true);
Feasible
T   = [3, 3]
R   = [1]
Sʲs = [[1, 3, 2, 1], [3, 2, 2], [1, 1, 2, 2, 1]]

j = 1:
   n = 1: S[1] + S[2] = T[1] + (R[1])
          1 + 3 = 3 + (1)
   n = 2: S[3] + S[4] = T[2]
          2 + 1 = 3
j = 2:
   n = 1: S[2] + S[3] = T[1] + (R[1])
          2 + 2 = 3 + (1)
   n = 2: S[1] = T[2]
          3 = 3
j = 3:
   n = 1: S[1] + S[2] + S[3] = T[1] + (R[1])
          1 + 1 + 2 = 3 + (1)
   n = 2: S[4] + S[5] = T[2]
          2 + 1 = 3
```
"""
function solve_subset_sum_variant(
            Sʲs :: AbstractVector{<:AbstractVector{<:Integer}},
            T   :: AbstractVector{<:Integer},
            R   :: AbstractVector{<:Integer};
            verbose :: Bool = false,
            optimizer = optimizer_with_attributes(HiGHS.Optimizer, "presolve"=>"off")
                        # (HiGHS' presolve is buggy for this problem; disable it)
    )

    J = length(Sʲs)
    N, Mʲs, L = length(T), length.(Sʲs), length(R)

    # --- fast-path checks: necessary conditions for feasibility ---
    # **fast path check #1:** if there is an element `sₘ` of any `Sʲ` in `Sʲs` such that
    # `sₘ` is larger than the largest possible value of `T` + sum(R), then the LHS of
    # ∑_{s ∈ S⁽ʲ⁾ₙ} s = tₙ + ∑_{r ∈ Rₙ} r involving `sₘ` will always be larger than any
    # choice of the RHS (since we assume R and T to be non-negative)
    maxT, sumR = maximum(T), sum(R)
    maxRHS = maxT+sumR
    for Sʲ in Sʲs
        if any(s -> s > maxRHS, Sʲ)
            verbose && printstyled("Infeasible: Sʲ > max(T) + max(R)\n", color=:red)
            return false, nothing
        end
    end

    # **fast path check #2:** the possible right-hand sides ``tₙ + ∑_{r ∈ Rₙ} r`` do not
    # depend on `j` and there are at most |T|^|R| possible variations over `n`. Let's call
    # each set of right-sides values `RHSₙᵖ` with `p = 1, …, |T|^|R|`. Suppose all RHS are
    # odd but only even elements exist in some `Sʲ`, then we can never find a solution.
    # It feels like this could be generalized, but so far, I haven't grokked how.
    if N^(L+1) < 1e4
        #  if total number of rhs-terms is less than 1e4, this is probably faster to check
        if any(Sʲ -> all(iseven, Sʲ), Sʲs) # some sets in Sʲs are entirely even
            rhs = similar(T)
            rhs_only_odd = true
            # loop over permutations of subsets of R (length L) into T-bins (length of N)
            for p in VectorProductIterator(fill(collect(1:N), L))
                rhs .= T # build a new permutaton of the rhs
                for (i, r) in enumerate(p)
                    rhs[r] += R[i]
                end
                all(isodd, rhs) || (rhs_only_odd = false; break)
            end
            if rhs_only_odd
                verbose && printstyled("Infeasible: some Sʲ has only even elements, but there are only odd RHS\n", color=:red)
                return false, nothing
            end
        end
    end

    # --- general & slow (NP) approach: check feasibility by binary integer programming ---
    # set up a model from `optimizer`: do it this way instead of `Model(optimizer)` to save
    # substantial time on instantiation; cf. 
    # https://discourse.julialang.org/t/reusing-a-jump-model-for-several-different-optimization-problems/122031/
    inner = MOI.instantiate(optimizer)
    model = direct_model(inner)
    set_string_names_on_creation(model, false)
    verbose || MOI.set(model, MOI.Silent(), true)

    @variable(model, x[1:sum(Mʲs), 1:N], Bin) # binary "indexing" into `Sʲs`
    if L ≠ 0
        @variable(model, y[1:L, 1:N], Bin) # binary "indexing" into `R`
    end

    # constraints on `x` and `y`
    for n in 1:N
        RHS = L ≠ 0 ? T[n] + sum(y[l,n]*R[l] for l in 1:L) : AffExpr(T[n])
        for j in 1:J
            mʲ_global_idxs = sum(Mʲs[1:j-1])+1:sum(Mʲs[1:j])
            @constraint(model, sum(x[mᵍ,n]*Sʲs[j][m] 
                                   for (m,mᵍ) in enumerate(mʲ_global_idxs)) == RHS)
        end
    end
    for j in 1:J
        mʲ_global_idxs = sum(Mʲs[1:j-1])+1:sum(Mʲs[1:j])
        @constraint(model, sum(x[mᵍ,n] for mᵍ in mʲ_global_idxs for n in 1:N) == Mʲs[j])
        for mᵍ in mʲ_global_idxs
            @constraint(model, sum(x[mᵍ,n] for n in 1:N) == 1)
        end
    end

    # constraints only on `y`
    if L ≠ 0
        @constraint(model, sum(y[l,n] for l in 1:L for n in 1:N) == L)
        for l in 1:L
            @constraint(model, sum(y[l,n] for n in 1:N) == 1)
        end
    end

    optimize!(model)
   
    feasible = is_solved_and_feasible(model)
    if !(feasible || termination_status(model) == INFEASIBLE)
        error(lazy"unexpected termination status: $(termination_status(m))")
    end

    verbose && _print_subset_sum_result(model, Sʲs, T, R)

    return feasible, model
end

function _print_subset_sum_result(io :: IO, model :: Model, Sʲs, T, R)
    J = length(Sʲs)
    N, Mʲs, L = length(T), length.(Sʲs), length(R)

    if is_solved_and_feasible(model)
        printstyled(io, "Feasible\n"; color=:green)
        printstyled(io, "T   = ", T, "\nR   = ", R, "\nSʲs = ", Sʲs, "\n\n"; 
                    color=:light_black)
        x = object_dictionary(model)[:x] # into `Sʲs`
        xv = value.(x)
        yv = if L ≠ 0
            y = object_dictionary(model)[:y] # into `R`
            value.(y)
        else
            Matrix{Float64}(undef, 0, 0)
        end
        for j in 1:J
            mʲ_global_idxs = sum(Mʲs[1:j-1])+1:sum(Mʲs[1:j])
            printstyled(io, "j = ", j, ":\n", color=:light_black)
            for n in 1:N
                printstyled(io, "   n = ", n, ": ", color=:light_black)
                s = ""
                v = ""
                for (m, mᵍ) in enumerate(mʲ_global_idxs)
                    if xv[mᵍ,n] > 0.5
                        s = s * (isempty(s) ? "" : " + ") * "S[$m]"
                        v = v * (isempty(v) ? "" : " + ") * string(Sʲs[j][m])
                    end
                end
                s = *(s, " = T[", string(n), "]")
                v = *(v, " = ", string(T[n]))
                if L ≠ 0
                    for l in 1:L
                        first = true
                        if yv[l,n] > 0.5
                            if first
                                s *= " + ("
                                v *= " + ("
                                first = false
                            else
                                s *= " + "
                                v *= " + "
                            end
                            s *= "R[$l]"
                            v *= string(R[l])
                        end
                        first || (s *= ")"; v *= ")")
                    end
                end
                println(io, s, "\n          ", v)
            end
        end
    elseif termination_status(model) == INFEASIBLE
        printstyled(io, "Infeasible\n"; color=:red)
    end
end
function _print_subset_sum_result(model :: Model, Sʲs, T, R)
    _print_subset_sum_result(stdout, model, Sʲs, T, R)
end

"""
    solve_subset_sum_variant_flexibleT(Sʲs, T, R, Nᵀ;
        verbose::Bool = false,
        model::Model = DEFAULT_HiGHS_MODEL, 
        preemptive_model_reset::Bool = true)

Akin to [`solve_subset_sum_variant`](@ref), but now allowing a subset selection of `T` also,
in `Nᵀ` sets, rather than N = |T| sets.

## Problem formulation
Given a target list `T` = (T₁, …, T_N), a set of multisets `Sʲs` = {S⁽ʲ⁾}_{j=1}^{J} with 
S⁽ʲ⁾ = {s⁽ʲ⁾₁, …, s⁽ʲ⁾_{M⁽ʲ⁾}}, as well as a set `R` = {r₁, …, r_L}, determine if there
exist disjoint subset selections {S⁽ʲ⁾ₙ}, {Tₙ}, and {Rₙ} with S⁽ʲ⁾ₙ ⊂ S⁽ʲ⁾, Tₙ ⊂ T, and
Rₙ ⊂ R such that:

    ∑_{s ∈ S⁽ʲ⁾ₙ} s = ∑_{t ∈ Tₙ} t + ∑_{r ∈ Rₙ} r

for n = 1, …, Nᵀ (with Nᵀ≤N) (and otherwise similar to [`solve_subset_sum_variant`](@ref)).
"""
function solve_subset_sum_variant_flexibleT(
            Sʲs :: AbstractVector{<:AbstractVector{<:Integer}},
            T   :: AbstractVector{<:Integer},
            R   :: AbstractVector{<:Integer},
            Nᵀ  :: Int = length(T);
            verbose :: Bool = false,
            optimizer = optimizer_with_attributes(HiGHS.Optimizer, "presolve"=>"off")
                        # (HiGHS' presolve is buggy for this problem; disable it)
    )

    J = length(Sʲs)
    N, Mʲs, L = length(T), length.(Sʲs), length(R)

    # TODO: fast-path checks as in `solve_subset_sum_variant`

    # set up model (faster than `Model(optimizer)`)
    inner = MOI.instantiate(optimizer)
    model = direct_model(inner)
    set_string_names_on_creation(model, false)
    verbose || MOI.set(model, MOI.Silent(), true)

    @variable(model, x[1:sum(Mʲs), 1:Nᵀ], Bin) # binary "indexing" into `Sʲs`
    @variable(model, z[1:N, 1:Nᵀ], Bin) # binary "indexing" into `T`
    if L ≠ 0
        @variable(model, y[1:L, 1:Nᵀ], Bin) # binary "indexing" into `R`
    end

    # TODO: TEST
    # constraints on `x` and `y` and `z`
    for n in 1:Nᵀ
        RHS_T = sum(z[m,n]*T[m] for m in 1:N)
        RHS_R = L ≠ 0 ? sum(y[l,n]*R[l] for l in 1:L) : AffExpr(0)
        RHS = RHS_T + RHS_R
        for j in 1:J
            mʲ_global_idxs = sum(Mʲs[1:j-1])+1:sum(Mʲs[1:j])
            @constraint(model, sum(x[mᵍ,n]*Sʲs[j][m] 
                                   for (m,mᵍ) in enumerate(mʲ_global_idxs)) == RHS)
        end
    end
    for j in 1:J
        mʲ_global_idxs = sum(Mʲs[1:j-1])+1:sum(Mʲs[1:j])
        @constraint(model, sum(x[mᵍ,n] for mᵍ in mʲ_global_idxs for n in 1:Nᵀ) == Mʲs[j])
        for mᵍ in mʲ_global_idxs
            @constraint(model, sum(x[mᵍ,n] for n in 1:Nᵀ) == 1)
        end
    end

    # constraints only on `z`
    if L ≠ 0
        @constraint(model, sum(z[m,n] for m in 1:N for n in 1:Nᵀ) == N)
        for m in 1:N
            @constraint(model, sum(z[m,n] for n in 1:Nᵀ) == 1)
        end
    end

    # constraints only on `y`
    if L ≠ 0
        @constraint(model, sum(y[l,n] for l in 1:L for n in 1:Nᵀ) == L)
        for l in 1:L
            @constraint(model, sum(y[l,n] for n in 1:Nᵀ) == 1)
        end
    end

    optimize!(model)
   
    feasible = is_solved_and_feasible(model)
    if !(feasible || termination_status(model) == INFEASIBLE)
        error(lazy"unexpected termination status: $(termination_status(m))")
    end

    # TODO: needs to be modified to allow `Nᵀ`
    # verbose && _print_subset_sum_result(model, Sʲs, T, R)

    return feasible, model
end