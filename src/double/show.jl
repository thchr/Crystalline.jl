# --- SU2 ---

function show(io::IO, ::MIME"text/plain", u::SU2)
    rows = _su2_rows(u; brackets=false)
    print(io, "SU(2) element:\n ", rows[1], "\n ", rows[2])
end
show(io::IO, u::SU2) = print(io, "SU2(", u.a, ", ", u.b, ")")

# the two rows of the matrix of `u`, with entries right-aligned by column
function _su2_rows(u::SU2; brackets::Bool=true)
    U = matrix(u)
    strs = [_su2_entry_string(U[i, j]) for i in 1:2, j in 1:2]
    ws = [maximum(textwidth, strs[:, j]) for j in 1:2]
    row(i) = join(lpad.(strs[i, :], ws), "  ")
    brackets || return (row(1), row(2))
    return ("┌ " * row(1) * " ┐", "└ " * row(2) * " ┘")
end

# The real and imaginary parts of the SU(2) element of a crystallographic operation take
# values in {0, ±1/2, ±√2/2, ±√3/2, ±1}: the half-angles are multiples of 30° or 45°, and
# the rotation axes are the directions of `SU2_BINARY_AXES` and the three- and four-fold
# axes. Such values are written exactly; any other value is rounded.
const SU2_SURDS = (1.0 => "1", 0.5 => "1/2", sqrt(2)/2 => "√2/2", sqrt(3)/2 => "√3/2")

# `|x|` as a string: exact if it is among `SU2_SURDS`, rounded otherwise
function _magnitude_string(x::Real)
    for (v, s) in SU2_SURDS
        isapprox(abs(x), v; atol=DEFAULT_ATOL) && return s
    end
    return string(round(abs(x); digits=4))
end

function _su2_entry_string(z::Number)
    r, s = reim(z)
    r = abs(r) < DEFAULT_ATOL ? 0.0 : r
    s = abs(s) < DEFAULT_ATOL ? 0.0 : s
    sr, ss = _magnitude_string(r), _magnitude_string(s)
    iszero(s) && return (signbit(r) ? "-" : "") * (iszero(r) ? "0" : sr)
    # `i` replaces a leading 1 ("i/2"), precedes a root ("i√2/2"), or follows a decimal
    istr = (signbit(s) ? "-" : "") * (startswith(ss, '1') ? "i" * ss[2:end] :
                                      startswith(ss, '√') ? "i" * ss : ss * "i")
    iszero(r) && return istr
    if sr == ss # equal magnitudes: factor out, e.g. "(1-i)√2/2" or "-(1+i)/2"
        sign = signbit(r) ? "-" : ""
        inner = signbit(r) == signbit(s) ? "(1+i)" : "(1-i)"
        return sign * inner * lstrip(sr, '1')
    end
    return (signbit(r) ? "-" : "") * sr * (signbit(s) ? "" : "+") * istr
end

# --- DSymOperation ---
# the spatial operation, followed by its SU(2) element, aligned with the matrix's last rows
function show(io::IO, ::MIME"text/plain", dop::DSymOperation{D}) where D
    _print_operation_header(io, dop)
    get(io, :compact, false) && return nothing
    println(io)
    op_rows = split(sprint(_print_operation_matrix, dop.op; context=io), '\n')
    su2_rows = _su2_rows(dop.su2)
    for (i, op_row) in enumerate(op_rows)
        print(io, op_row)
        j = i - (D - 2)
        j ≥ 1 && printstyled(io, j == 2 ? ", " : "  ", su2_rows[j]; color=:light_black)
        i ≠ D && println(io)
    end
    return nothing
end
