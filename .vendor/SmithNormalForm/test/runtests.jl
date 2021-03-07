using LinearAlgebra
using SparseArrays
using Test
using SmithNormalForm

@testset "Bezout" begin
    a = 12
    b = 42

    s,t,g = SmithNormalForm.bezout(a,b)
    @test s*a + t*b == g
    g,s,t = gcdx(a,b)
    @test s*a + t*b == g

    s,t,g = SmithNormalForm.bezout(b,a)
    @test s*b + t*a == g
    g,s,t = gcdx(b,a)
    @test s*b + t*a == g
end

@testset "Smith Normal Form" begin
    M=[ 1  0  0  0  0  0 ;
        0  0  0  0  0  0 ;
       -1  0  1  0  0  1 ;
        0 -1  0 -1  0  0 ;
        0  0  0  1  1 -1 ;
        0  1 -1  0 -1  0 ]

    P, Q, A, Pinv, Qinv = snf(M)
    @test !issparse(A)
    @test P*A*Q == M
    @test inv(P) == Pinv
    @test inv(Q) == Qinv
    @testset "dense diag" for i in 1:4
        @test A[i,i] == 1
    end

    P, Q, A, Pinv, Qinv = snf(dropzeros(sparse(M)))
    @test issparse(A)
    @test P*A*Q == M
    @test inv(collect(P)) == Pinv
    @test inv(collect(Q)) == Qinv
    @testset "sparse diag" for i in 1:4
        @test A[i,i] == 1
    end

    P, Q, A, Pinv, Qinv = snf(sparse(M), inverse=false)
    @test issparse(A)
    @test iszero(Pinv)
    @test iszero(Qinv)
    @test P*A*Q == M
    @testset "sparse diag" for i in 1:4
        @test A[i,i] == 1
    end

    # Factorization
    M =[22  30  64  93  36;
        45  42  22  11  67;
        21   1  35  45  42]
    n, m = size(M)
    F = smith(M)
    @test !issparse(F.SNF)
    @test size(F.S) == (n,n)
    @test size(F.T) == (m,m)
    @test size(F.Sinv) == (n,n)
    @test size(F.Tinv) == (m,m)
    @test F.S*diagm(F)*F.T == M
    @test F.S*F.Sinv == Matrix{eltype(F)}(I, n, n)
    @test F.T*F.Tinv == Matrix{eltype(F)}(I, m, m)

    F = smith(sparse(M), inverse=false)
    @test issparse(F.SNF)
    @test size(F.S) == (n,n)
    @test size(F.T) == (m,m)
    @test F.S*diagm(F)*F.T == M
    @test iszero(F.Sinv)
    @test iszero(F.Tinv)

    # https://en.wikipedia.org/wiki/Smith_normal_form#Example
    M = [2 4 4; -6 6 12; 10 -4 -16]
    F = smith(M, inverse=true)
    @test F.SNF == [2, 6, 12]
    @test F.S*diagm(F)*F.T == M
end
