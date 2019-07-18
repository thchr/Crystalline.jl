import Base: show

# Crystalline lattice
struct Crystal{N}
    R::NTuple{N,Vector{Float64}}
end
Crystal(R,type) = Crystal{length(R)}(R,type)
Crystal(R) = Crystal{length(R)}(R,"")
basis(C::Crystal) = C.R
dim(C::Crystal{N}) where N = N
function show(io::IO, ::MIME"text/plain", C::Crystal)
    print(io, "$(dim(C))D Crystal:")
    print(io, " ($(crystalsystem(C)))");
    for (i,R) in enumerate(basis(C))
        print(io, "\n   R$(i): "); print(io, R); 
    end
end
norms(C::Crystal) = norm.(basis(C))
_angle(rA,rB) = acos(dot(rA,rB)/(norm(rA)*norm(rB)))
function angles(C::Crystal{N}) where N
    bvecs = basis(C)
    γ = _angle(bvecs[1], bvecs[2])
    if N == 3
        α = _angle(bvecs[2], bvecs[3])
        β = _angle(bvecs[3], bvecs[1])
        return α,β,γ
    end
    return γ
end

# symmetry operations
struct SymOperation
    shorthand::String
    matrix::Matrix{Float64}
end
SymOperation(s::String) = SymOperation(s,xyzt_op(s))
matrix(op::SymOperation) = op.matrix
shorthand(op::SymOperation) = op.shorthand
dim(op::SymOperation) = size(matrix(op),1)
function show(io::IO, ::MIME"text/plain", op::SymOperation) 
    print(io, "   (", op.shorthand, ")\n")
    Base.print_matrix(IOContext(io, :compact=>true), op.matrix, "   ")
    print(io, '\n')
end
getindex(op::SymOperation, keys...) = matrix(op)[keys...]   # allows direct indexing into an op::SymOperation like op[1,2] to get matrix(op)[1,2]
lastindex(op::SymOperation, d::Int64) = size(matrix(op), d) # allows using `end` in indices
pg(op::SymOperation) = op[:,1:end-1]        # point group part of an operation
translation(op::SymOperation) = op[:,end]   # translation part of an opeation
issymmorph(op::SymOperation) = iszero(translation(op))

# space group
struct SpaceGroup
    num::Integer
    operations::Vector{SymOperation}
    dim::Integer
end
num(sg::SpaceGroup) = sg.num
operations(sg::SpaceGroup) = sg.operations
dim(sg::SpaceGroup) = sg.dim
order(sg::SpaceGroup) = length(operations(sg))
function show(io::IO, ::MIME"text/plain", sg::SpaceGroup)
    Nops = order(sg)
    groupprefix = dim(sg) == 3 ? "Space" : (dim(sg) == 2 ? "Plane" : nothing)
    println(io, groupprefix, " group #", num(sg))
    for (i,op) in enumerate(operations(sg))
        show(io, "text/plain", op)
        if i < Nops; print(io, '\n'); end
    end
end
function show(io::IO, ::MIME"text/plain", sgs::Vector{SpaceGroup})
    Nsgs = order(sgs)
    for (i,sg) in enumerate(sgs); 
        show(io, "text/plain", sg); 
        if i < Nsgs; print(io, '\n'); end
    end
end