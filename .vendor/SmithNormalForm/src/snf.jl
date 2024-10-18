# Smith Normal Form

function divisable(y::R, x::R) where {R}
  x == zero(R) && return y == zero(R)
  return div(y,x)*x == y
end

function divide(y::R, x::R) where {R}
    if x != -one(R)
        return R(div(y,x))
    else
        return R(y * x)
    end
end

function rcountnz(X::AbstractMatrix{R}, i) where {R}
    c = 0
    z = zero(R)
    @inbounds for row in eachrow(X)
        if row[i] != z
           c += 1
        end
    end
    return c
end

function ccountnz(X::AbstractMatrix{R}, j) where {R}
    c = 0
    z = zero(R)
    @inbounds for col in eachcol(X)
        if col[j] != z
           c += 1
        end
    end
    return c
end

function rswap!(M::AbstractMatrix, r1::Int, r2::Int)
    r1 == r2 && return M
    @inbounds for col in eachcol(M)
        tmp = col[r1]
        col[r1] = col[r2]
        col[r2] = tmp
    end
    return M
end

function cswap!(M::AbstractMatrix, c1::Int, c2::Int)
    c1 == c2 && return M
    @inbounds for row in eachrow(M)
        tmp = row[c1]
        row[c1] = row[c2]
        row[c2] = tmp
    end
    return M
end

function rowelimination(D::AbstractMatrix{R}, a::R, b::R, c::R, d::R, i::Int, j::Int) where {R}
    @inbounds for col in eachcol(D)
        t = col[i]
        s = col[j]
        col[i] = a*t + b*s
        col[j] = c*t + d*s
    end
    return D
end

function colelimination(D::AbstractMatrix{R}, a::R, b::R, c::R, d::R, i::Int, j::Int) where {R}
    @inbounds for row in eachrow(D)
        t = row[i]
        s = row[j]
        row[i] = a*t + b*s
        row[j] = c*t + d*s
    end
    return D
end

function rowpivot(U::AbstractArray{R,2},
                  Uinv::AbstractArray{R,2},
                  D::AbstractArray{R,2},
                  i, j; inverse=true) where {R}
    for k in reverse(axes(D, 1))
        b = D[k,j]
        (iszero(b) || i == k) && continue
        a = D[i,j]

        g, s, t = bezout(a, b)
        x = divide(a, g)
        y = divide(b, g)

        rowelimination(D, s, t, -y, x, i, k)
        inverse && rowelimination(Uinv, s, t, -y, x, i, k)
        colelimination(U, x, y, -t, s, i, k)
    end
end

function colpivot(V::AbstractArray{R,2},
                  Vinv::AbstractArray{R,2},
                  D::AbstractArray{R,2},
                  i, j; inverse=true) where {R}
    for k in reverse(axes(D, 2))
        b = D[i,k]
        (iszero(b) || j == k) && continue
        a = D[i,j]

        g, s, t = bezout(a, b)
        x = divide(a, g)
        y = divide(b, g)

        colelimination(D, s, t, -y, x, j, k)
        inverse && colelimination(Vinv, s, t, -y, x, j, k)
        rowelimination(V, x, y, -t, s, j, k)
    end
end

function smithpivot(U::AbstractMatrix{R},
                    Uinv::AbstractMatrix{R},
                    V::AbstractMatrix{R},
                    Vinv::AbstractMatrix{R},
                    D::AbstractMatrix{R},
                    i, j; inverse=true) where {R}

    pivot = D[i,j]
    @assert pivot != zero(R) "Pivot cannot be zero"
    while ccountnz(D,i) > 1 || rcountnz(D,j) > 1
        colpivot(V, Vinv, D, i, j, inverse=inverse)
        rowpivot(U, Uinv, D, i, j, inverse=inverse)
    end
end

function init(M::AbstractSparseMatrix{R,Ti}; inverse=true) where {R, Ti}
    D = copy(M)
    rows, cols = size(M)

    U = sparse(R, Ti, I, rows, rows)
    Uinv = inverse ? copy(U) : spzeros(R, Ti, 0, 0)

    V = sparse(R, Ti, I, cols, cols)
    Vinv = inverse ? copy(V) : spzeros(R, Ti, 0, 0)

    return U, V, D, Uinv, Vinv
end

function init(M::AbstractMatrix{R}; inverse=true) where {R}
    D = Matrix{R}(copy(M))
    rows, cols = size(M)

    U = Matrix{R}(I, rows, rows)
    Uinv = inverse ? copy(U) : zeros(R, 0, 0)

    V = Matrix{R}(I, cols, cols)
    Vinv = inverse ? copy(V) : zeros(R, 0, 0)

    return U, V, D, Uinv, Vinv
end

formatmtx(M) =  size(M,1) == 0 ? "[]" : repr(collect(M); context=IOContext(stdout, :compact => true))

function snf(M::AbstractMatrix{R}; inverse=true) where {R}
    U, V, D, Uinv, Vinv = init(M, inverse=inverse)

    T = eltype(eachindex(M))
    t = one(T)
    for j in axes(M, 2) # over column indices
        # @debug "Working on column $j out of $cols" D=formatmtx(D)
        rcountnz(D,j) == 0 && continue

        prow = one(T)
        if D[t,t] != zero(R)
            prow = t
        else
            # the good pivot row for j-th column is the one that has fewest elements
            rsize = typemax(T)
            for i in axes(M, 1) # over row indices
                if D[i,j] != zero(R)
                    c = count(!iszero, view(D, i, :); init=zero(T))
                    if c < rsize
                        rsize = c
                        prow = i
                    end
                end
            end
        end

        # @debug "Pivot Row selected: t = $t, pivot = $prow" D=formatmtx(D)
        rswap!(D, t, prow)
        inverse && rswap!(Uinv, t, prow)
        cswap!(U, t, prow)

        # @debug "Performing the pivot step at (t=$t, j=$j)" D=formatmtx(D)
        smithpivot(U, Uinv, V, Vinv, D, t, j, inverse=inverse)

        cswap!(D, t, j)
        inverse && cswap!(Vinv, t, j)
        rswap!(V, t, j)

        t += one(T)
    end

    # Make sure that d_i is divisible be d_{i+1}.
    r = minimum(size(D))
    one_r = one(typeof(r))
    pass = true
    while pass
        pass = false
        for i in one_r:r-one_r
            ip1 = i + one_r
            divisable(D[ip1,ip1], D[i,i]) && continue
            pass = true
            D[ip1,i] = D[ip1,ip1]

            colelimination(Vinv, one(R), one(R), zero(R), one(R), i, ip1)
            rowelimination(V, one(R), zero(R), -one(R), one(R), i, ip1)

            smithpivot(U, Uinv, V, Vinv, D, i, i, inverse=inverse)
        end
    end

    # To guarantee SNFⱼ = Λⱼ ≥ 0 we absorb the sign of Λ into T and T⁻¹, s.t.
    #    Λ′ = Λ*sign(Λ),   T′ = sign(Λ)*T,    and    T⁻¹′ = T⁻¹*sign(Λ),
    # with the convention that sign(0) = 1. Then we still have that X = SΛT = SΛ′T′
    # and also that Λ = S⁻¹XT⁻¹ ⇒ Λ′ = S⁻¹XT⁻¹′.
    for j in one(T):minimum(size(M))
        Λⱼ = D[j,j]
        if Λⱼ < zero(R)
            @views V[j,:] .*= -one(R)        # T′   = sign(Λ)*T    [rows]
            if inverse
                @views Vinv[:,j] .*= -one(R) # T⁻¹′ = T⁻¹*sign(Λ)  [columns]
            end
            D[j,j] = abs(Λⱼ)                 # Λ′ = Λ*sign(Λ)
        end
    end

    if issparse(D)
        return dropzeros!(U), dropzeros!(V), dropzeros!(D), dropzeros!(Uinv), dropzeros!(Vinv)
    else
        return U, V, D, Uinv, Vinv
    end
end
