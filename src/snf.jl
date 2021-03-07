# Smith Normal Form

function divisable(y::R, x::R ) where {R}
  x == zero(R) && return y == zero(R)
  return div(y,x)*x == y
end

function divide(y::R, x::R) where {R}
    if x != -one(R)
        return div(y,x)
    else
        return y * x
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
    for k in reverse!(findall(!iszero, view(D, :, j)))
        a = D[i,j]
        b = D[k,j]

        i == k && continue

        s, t, g = bezout(a, b)
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
    for k in reverse!(findall(!iszero, view(D, i, :)))
        a = D[i,j]
        b = D[i,k]

        j == k && continue

        s, t, g = bezout(a, b)
        x = divide(a, g)
        y = divide(b, g)

        colelimination(D, s, t, -y, x, j, k)
        inverse && colelimination(Vinv, s, t, -y, x, j, k)
        rowelimination(V, x, y, -t, s, j, k)
    end
end

function smithpivot(U::AbstractArray{R,2},
                    Uinv::AbstractArray{R,2},
                    V::AbstractArray{R,2},
                    Vinv::AbstractArray{R,2},
                    D::AbstractArray{R,2},
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

    U = spzeros(R, rows, rows)
    for i in 1:rows
        U[i,i] = one(R)
    end
    Uinv = inverse ? copy(U) : spzeros(R, 0, 0)

    V = spzeros(R, cols, cols)
    for i in 1:cols
        V[i,i] = one(R)
    end
    Vinv = inverse ? copy(V) : spzeros(R, 0, 0)

    return U, V, D, Uinv, Vinv
end

function init(M::AbstractMatrix{R}; inverse=true) where {R}
    D = copy(M)
    rows, cols = size(M)

    U = zeros(R, rows, rows)
    for i in 1:rows
        U[i,i] = one(R)
    end
    Uinv = inverse ? copy(U) : zeros(R, 0, 0)

    V = zeros(R, cols, cols)
    for i in 1:cols
        V[i,i] = one(R)
    end
    Vinv = inverse ? copy(V) : zeros(R, 0, 0)

    return U, V, D, Uinv, Vinv
end

formatmtx(M) =  size(M,1) == 0 ? "[]" : repr(collect(M); context=IOContext(stdout, :compact => true))

function snf(M::AbstractMatrix{R}; inverse=true) where {R}
    rows, cols = size(M)
    U, V, D, Uinv, Vinv = init(M, inverse=inverse)

    t = 1
    for j in 1:cols
        @debug "Working on column $j out of $cols" D=formatmtx(D)

        rcountnz(D,j) == 0 && continue

        prow = 1
        if D[t,t] != zero(R)
            prow = t
        else
            # Good pivot row for j-th column is the one
            # that have a smallest number of elements
            rsize = typemax(R)
            for i in 1:rows
                if D[i,j] != zero(R)
                    c = count(!iszero, view(D, i, :))
                    if c < rsize
                        rsize = c
                        prow = i
                    end
                end
            end
        end

        @debug "Pivot Row selected: t = $t, pivot = $prow" D=formatmtx(D)
        rswap!(D, t, prow)
        inverse && rswap!(Uinv, t, prow)
        cswap!(U, t, prow)

        @debug "Performing the pivot step at (t=$t, j=$j)" D=formatmtx(D)
        smithpivot(U, Uinv, V, Vinv, D, t, j, inverse=inverse)

        cswap!(D, t, j)
        inverse && cswap!(Vinv, t, j)
        rswap!(V, t, j)

        t += 1

        @logmsg (Base.CoreLogging.Debug-1) "Factorization" D=formatmtx(D) U=formatmtx(U) V=formatmtx(V) U⁻¹=formatmtx(Uinv) V⁻¹=formatmtx(Vinv)
    end

    # Make sure that d_i is divisible be d_{i+1}.
    r = minimum(size(D))
    pass = true
    while pass
        pass = false
        for i in 1:r-1
            divisable(D[i+1,i+1], D[i,i]) && continue
            pass = true
            D[i+1,i] = D[i+1,i+1]

            colelimination(Vinv, one(R), one(R), zero(R), one(R), i, i+1)
            rowelimination(V, one(R), zero(R), -one(R), one(R), i, i+1)

            smithpivot(U, Uinv, V, Vinv, D, i, i, inverse=inverse)
        end
    end

    # To guarantee SNFⱼ = Λⱼ ≥ 0 we absorb the sign of Λ into T and T⁻¹, s.t.
    #    Λ′ = Λ*sign(Λ),   T′ = sign(Λ)*T,    and    T⁻¹′ = T⁻¹*sign(Λ),
    # with the convention that sign(0) = 1. Then we still have that X = SΛT = SΛ′T′
    # and also that Λ = S⁻¹XT⁻¹ ⇒ Λ′ = S⁻¹XT⁻¹′.
    for j in 1:rows
        j > cols && break
        Λⱼ = D[j,j]
        if Λⱼ < 0
            @views V[j,:]    .*= -1 # T′   = sign(Λ)*T    [rows]
            if inverse
                @views Vinv[:,j] .*= -1 # T⁻¹′ = T⁻¹*sign(Λ)  [columns]
            end
            D[j,j] = abs(Λⱼ)        # Λ′ = Λ*sign(Λ)
        end
    end
    @logmsg (Base.CoreLogging.Debug-1) "Factorization" D=formatmtx(D) U=formatmtx(U) V=formatmtx(V) U⁻¹=formatmtx(Uinv) V⁻¹=formatmtx(Vinv)

    if issparse(D)
        return dropzeros!(U), dropzeros!(V), dropzeros!(D), dropzeros!(Uinv), dropzeros!(Vinv)
    else
        return U, V, D, Uinv, Vinv
    end
end
