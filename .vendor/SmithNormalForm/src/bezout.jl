"""
Calculates Bézout coefficients.

This is analogous to `gcdx`, and solves the same problem, but does not necessarily return
identical solutions.
"""
Base.@assume_effects :terminates_locally function bezout(a::R, b::R) where {R}
    rev = a < b
    x, y = rev ? (a,b) : (b,a)

    s0, s1 = oneunit(R), zero(R)
    t0, t1 = zero(R), oneunit(R)

    while !iszero(y)
        q, r = divrem(x, y)
        x, y = y, r
        s0, s1 = s1, s0 - q * s1
        t0, t1 = t1, t0 - q * t1
    end

    s, t = rev ? (s0, t0) : (t0, s0)
    g = x

    if g == a # differs from `gcdx` here (also in `rev` clauses above)
        s = one(R)
        t = zero(R)
    elseif g == -a
        s = -one(R)
        t = zero(R)
    end

    return g, s, t
end
