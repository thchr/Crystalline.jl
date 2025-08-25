using Test
using Crystalline.Jagged: JaggedVector

@testset "JaggedVector" begin
    # tests generated with Claude Sonnet 4
    @testset "Construction" begin
        @testset "Basic constructor with data and offsets" begin
            data = [1, 2, 3, 4, 5]
            offsets = [1, 3, 4, 6]
            jv = JaggedVector(data, offsets)
            
            @test jv isa JaggedVector{Int}
            @test length(jv) == 3
            @test jv[1] == [1, 2]
            @test jv[2] == [3]
            @test jv[3] == [4, 5]
        end
        
        @testset "Type-specified constructor" begin
            data = [1.0, 2.0, 3.0]
            offsets = [1, 2, 4]
            jv = JaggedVector{Float64}(data, offsets)
            
            @test jv isa JaggedVector{Float64}
            @test eltype(jv) === typeof(@view data[1:2])
            @test jv[1] ≈ [1.0]
            @test jv[2] ≈ [2.0, 3.0]
        end
        
        @testset "Empty constructor" begin
            jv_empty = JaggedVector{Int}(Int[], [1])
            @test jv_empty == JaggedVector{Int}() == empty(jv_empty)
            @test JaggedVector() == JaggedVector(Any[], [1])
            @test empty(JaggedVector{Int}([[1,2], [3]]), Float64) == JaggedVector{Float64}()
        end
            
        @testset "Constructor validation" begin
            @test_throws ArgumentError JaggedVector([1, 2, 3], [1, 3, 2, 4]) # unsorted offsets
            @test_throws ArgumentError JaggedVector([1, 2, 3], [2, 3, 4]) # first offset not 1
            @test_throws ArgumentError JaggedVector([1, 2, 3], [1, 2, 3]) # last offset incorrect
            @test_throws ArgumentError JaggedVector([1, 2, 3], [1, 2, 5])
        end
        
        @testset "Constructor from lengths" begin
            lengths = [3, 0, 2, 1]
            jv = JaggedVector{Int}(lengths)
            
            @test length(jv) == 4
            @test length(jv[1]) == 3
            @test length(jv[2]) == 0
            @test length(jv[3]) == 2
            @test length(jv[4]) == 1
            @test length(jv.data) == 6
            
            # test negative length error
            @test_throws ErrorException JaggedVector{Int}([1, -1, 2])
        end
        
        @testset "`zeros` constructor" begin
            lengths = [2, 1, 3]
            jv = zeros(JaggedVector{Float64}, lengths)
            
            @test length(jv) == 3
            @test all(jv[1] .== 0.0)
            @test all(jv[2] .== 0.0)
            @test all(jv[3] .== 0.0)
        end
        
        @testset "`zero` function" begin
            data = [1, 2, 3, 4]
            offsets = [1, 3, 5]
            jv = JaggedVector(data, offsets)
            jv_zero = zero(jv)
            
            @test length(jv_zero) == length(jv)
            @test all(jv_zero.data .== 0)
            @test jv_zero.offsets == jv.offsets
            @test jv_zero.offsets !== jv.offsets # should be a copy
        end
        
        @testset "Constructor from vector of vectors" begin
            vv = [[1, 2, 3], [4], [5, 6]]
            jv = JaggedVector(vv)
            
            @test jv isa JaggedVector{Int}
            @test length(jv) == 3
            @test jv[1] == [1, 2, 3]
            @test jv[2] == [4]
            @test jv[3] == [5, 6]
            
            # test type conversion
            jv_float = JaggedVector{Float64}(vv)
            @test jv_float isa JaggedVector{Float64}
            @test jv_float[1] == [1.0, 2.0, 3.0]
        end
        
        @testset "Constructor from iterable of iterables" begin
            iter = ([1, 2], [3, 4, 5], [6])
            jv = JaggedVector(iter)
            
            @test length(jv) == 3
            @test jv[1] == [1, 2]
            @test jv[2] == [3, 4, 5]
            @test jv[3] == [6]
        end
    end

    @testset "Equivalence with `Vector{Vector{T}}`" begin
        @testset "Basic equivalence" begin
            vs = [[1, 2, 3], [4], [5, 6]]
            jv = JaggedVector(vs)
            
            # test structural equivalence
            @test length(jv) == length(vs)
            @test all(jv[i] == vs[i] for i in eachindex(vs))
            
            # test that elements are views but equal to original vectors
            for i in eachindex(vs)
                @test jv[i] == vs[i]
                @test jv[i] isa SubArray # should be a view, not a copy
            end
            
            # test iteration equivalence
            jv_collected = [copy(v) for v in jv] # copy to get concrete vectors
            @test jv_collected == vs
        end
        
        @testset "Empty and mixed-length vectors" begin
            vs = [Int[], [1], [2, 3, 4, 5], Int[], [6, 7]]
            jv = JaggedVector(vs)
            
            @test length(jv) == length(vs)
            @test all(jv[i] == vs[i] for i in eachindex(vs))
            @test length(jv[1]) == 0
            @test length(jv[2]) == 1
            @test length(jv[3]) == 4
            @test length(jv[4]) == 0
            @test length(jv[5]) == 2
        end
        
        @testset "Type promotion" begin
            vs = [[1, 2], [3, 4]]
            jv_int = JaggedVector{Int}(vs)
            jv_float = JaggedVector{Float64}(vs)
            
            @test eltype(jv_int.data) == Int
            @test eltype(jv_float.data) == Float64
            @test all(jv_int[i] == vs[i] for i in eachindex(vs))
            @test all(jv_float[i] == vs[i] for i in eachindex(vs))
        end
        
        @testset "Round-trip conversion" begin
            vs = [[1, 2, 3], [4], [5, 6]]
            jv = JaggedVector(vs)
            vs_recovered = Vector{Vector{Int}}(jv)
            
            @test vs_recovered == vs
            @test vs_recovered !== vs # should be a new object
            
            # test that modifying recovered doesn't affect original
            vs_recovered[1][1] = 999
            @test vs[1][1] == 1 # `vs` unchanged
            @test jv[1][1] == 1 # `jv` unchanged
        end
    end
    
    @testset "Construction from `AbstractVector{<:AbstractVector}`" begin
        @testset "Vector of views" begin
            data = [1, 2, 3, 4, 5, 6, 7, 8]
            views = [view(data, 1:3), view(data, 4:4), view(data, 5:8)]
            jv = JaggedVector(views)
            
            @test length(jv) == 3
            @test jv[1] == [1, 2, 3]
            @test jv[2] == [4]
            @test jv[3] == [5, 6, 7, 8]
            
            # test type inference
            @test jv isa JaggedVector{Int}
        end
        
        @testset "Mixed `AbstractVector` types" begin
            vec1 = [1, 2, 3]
            vec2 = view([4, 5, 6, 7], 1:2) # [4, 5]
            vec3 = (8:10) # UnitRange
            
            mixed_vecs = AbstractVector{Int}[vec1, vec2, collect(vec3)]
            jv = JaggedVector(mixed_vecs)
            
            @test length(jv) == 3
            @test jv[1] == [1, 2, 3]
            @test jv[2] == [4, 5]
            @test jv[3] == [8, 9, 10]
        end
        
        @testset "Type specification with `AbstractVector`" begin
            views = [view([1.0, 2.0, 3.0], 1:2), view([4.0, 5.0], :)]
            jv_float = JaggedVector{Float64}(views)
            jv_int = JaggedVector{Int}(views)
            
            @test jv_float isa JaggedVector{Float64}
            @test jv_int isa JaggedVector{Int}
            @test jv_int[1] == [1, 2]
            @test jv_int[2] == [4, 5]
        end
        
        @testset "Empty `AbstractVectors`" begin
            empty_view = view([1, 2, 3], 2:1) # empty view
            regular_empty = Int[]
            vecs = AbstractVector{Int}[regular_empty, [1, 2], empty_view, [3]]
            
            jv = JaggedVector(vecs)
            @test length(jv) == 4
            @test isempty(jv[1])
            @test jv[2] == [1, 2]
            @test isempty(jv[3])
            @test jv[4] == [3]
        end
        
        @testset "Conversion method equivalence" begin
            data = [1, 2, 3, 4, 5]
            views = [view(data, 1:2), view(data, 3:3), view(data, 4:5)]
            
            # test both constructor and convert
            jv1 = JaggedVector(views)
            jv2 = convert(JaggedVector{Int}, views)
            
            @test jv1.data == jv2.data
            @test jv1.offsets == jv2.offsets
            @test all(jv1[i] == jv2[i] for i in eachindex(jv1))
        end
    end
       
    @testset "AbstractArray interface" begin
        data = [1, 2, -8, -7, -6, -5]
        offsets = [1, 3, 3, 7] # i.e., defines [[1,2], Int[], [3,4,5,6]]
        jv = JaggedVector(data, offsets)
        
        @testset "Indexing" begin
            @test jv[1] == [1, 2]
            @test jv[2] == []
            @test jv[3] == [-8, -7, -6, -5]
            
            @test_throws BoundsError jv[0]
            @test_throws BoundsError jv[4]
        end
        
        @testset "Size and length" begin
            @test size(jv) == (3,)
            @test length(jv) == 3
        end
        
        @testset "`setindex!`" begin
            jv_copy = copy(jv)
            
            # valid assignment (same length)
            jv_copy[1] = [10, 20]
            @test jv_copy[1] == [10, 20]
            
            # empty to empty
            jv_copy[2] = []
            @test jv_copy[2] == []
            
            # valid assignment for longer vector
            jv_copy[3] = [30, 40, 50, 60]
            @test jv_copy[3] == [30, 40, 50, 60]
            
            # assignment with different lengths than original iterants
            jv_copy = copy(jv) # reset
            jv_copy[1] = [11, 22, 33] # increase length
            @test jv_copy == [[11, 22, 33], Int[], [-8, -7, -6, -5]]
            jv_copy[2] = [-1]         # increase length from 0
            @test jv_copy == [[11, 22, 33], [-1], [-8, -7, -6, -5]]
            jv_copy[1] = [0]          # decrease length of first element
            @test jv_copy == [[0], [-1], [-8, -7, -6, -5]]
            jv_copy[3] = [-10, -20]   # decrease length of last element
            @test jv_copy == [[0], [-1], [-10, -20]]
            jv_copy[3] = Int[]        # decrease length to 0
            @test jv_copy == [[0], [-1], Int[]]
            foreach(1:3) do i; jv_copy[i] = Int[]; end # decrease all lengths to 0
            @test jv_copy == [Int[], Int[], Int[]]

            # out of bounds assignment
            @test_throws BoundsError jv_copy[4] = [1]
        end
        
        @testset "Iteration" begin
            result = Vector{Vector{Int}}()
            for v in jv
                push!(result, copy(v))
            end
            @test result == [[1, 2], Int[], [-8, -7, -6, -5]]
            
            # test iterate function directly
            iter_result = iterate(jv)
            @test iter_result[1] == [1, 2]
            
            iter_result = iterate(jv, iter_result[2])
            @test iter_result[1] == Int[]
            
            iter_result = iterate(jv, iter_result[2])
            @test iter_result[1] == [-8, -7, -6, -5]
            
            @test iterate(jv, iter_result[2]) === nothing

            _vs = [[1,2,3], [7,-1,-20], [0,], Int[], [32,41,-102], Int[]]
            _jv = JaggedVector(_vs)
            @test mapreduce(sum, +, _vs) == mapreduce(sum, +, _jv)
            @test sum(Iterators.flatten(_jv)) == mapreduce(sum, +, _vs)
        end
        
        @testset "`similar`" begin
            jv_sim = similar(jv)
            @test jv_sim isa JaggedVector{Int}
            @test size(jv_sim) == size(jv)
            @test jv_sim.offsets == jv.offsets
            @test jv_sim.offsets !== jv.offsets # should be a copy
            
            jv_sim_float = similar(jv, Float64)
            @test jv_sim_float isa JaggedVector{Float64}
            @test eltype(jv_sim_float.data) == Float64
        end
        
        @testset "checkbounds" begin
            @test checkbounds(Bool, jv, 1) == true
            @test checkbounds(Bool, jv, 3) == true
            @test checkbounds(Bool, jv, 0) == false
            @test checkbounds(Bool, jv, 4) == false
        end
    end
    
    @testset "Mutation operations" begin
        @testset "`push!`" begin
            jv = JaggedVector([[1, 2], [3]])
            push!(jv, [4, 5, 6])
            
            @test length(jv) == 3
            @test jv[3] == [4, 5, 6]
            
            push!(jv, [7], [8, 9]) # multiple push
            @test length(jv) == 5
            @test jv[4] == [7]
            @test jv[5] == [8, 9]
        end
        
        @testset "`append!`" begin
            jv = JaggedVector([[1, 2]])
            append!(jv, [[3, 4], [5]])
            
            @test length(jv) == 3
            @test jv[2] == [3, 4]
            @test jv[3] == [5]
            
            append!(jv, [[6]], [[7, 8], [9]]) # multiple append!
            @test length(jv) == 6
            @test jv[6] == [9]
            
            len_before = length(jv) # empty append!
            append!(jv, [])
            @test length(jv) == len_before
        end
        
        @testset "`pop!`" begin
            jv = JaggedVector([[1, 2], [3, 4, 5]])
            popped = pop!(jv)
            
            @test popped == [3, 4, 5]
            @test length(jv) == 1
            @test jv[1] == [1, 2]
        end
        
        @testset "Unsupported operations" begin
            jv = JaggedVector([[1, 2], [3]])
            
            @test_throws ErrorException pushfirst!(jv, [0])
            @test_throws ErrorException popfirst!(jv)
            @test_throws ErrorException popat!(jv, 1)
            @test_throws ErrorException deleteat!(jv, 1)
            @test_throws ErrorException keepat!(jv, [1])
        end
    end
    
    @testset "Arithmetic operations" begin
        jv1 = JaggedVector([[1, 2], [3], [4, 5, 6]])
        jv2 = JaggedVector([[2, 3], [1], [1, 1, 1]])
        
        @testset "Addition" begin
            jv_sum = jv1 + jv2
            @test jv_sum[1] == [3, 5]
            @test jv_sum[2] == [4]
            @test jv_sum[3] == [5, 6, 7]
            
            # mismatched offsets
            jv3 = JaggedVector([[1, 2, 3], [4]])
            @test_throws ErrorException jv1 + jv3
        end
        
        @testset "Subtraction" begin
            jv_diff = jv1 - jv2
            @test jv_diff[1] == [-1, -1]
            @test jv_diff[2] == [2]
            @test jv_diff[3] == [3, 4, 5]
            
            # unary minus
            jv_neg = -jv1
            @test jv_neg[1] == [-1, -2]
            @test jv_neg[2] == [-3]
            @test jv_neg[3] == [-4, -5, -6]
        end
        
        @testset "Scalar multiplication and division" begin
            jv_scaled = 2 * jv1
            @test jv_scaled[1] == [2, 4]
            @test jv_scaled[2] == [6]
            @test jv_scaled[3] == [8, 10, 12]
            
            jv_scaled2 = jv1 * 3
            @test jv_scaled2[1] == [3, 6]
            @test jv_scaled2[2] == [9]
            @test jv_scaled2[3] == [12, 15, 18]
            
            jv_div = jv_scaled / 2
            @test jv_div[1] == [1.0, 2.0]
            @test jv_div[2] == [3.0]
            @test jv_div[3] == [4.0, 5.0, 6.0]
        end
    end
    
    @testset "Utility functions" begin
        jv = JaggedVector([[1, 2], [3], [4, 5, 6]])
        
        @testset "`parent`" begin
            @test parent(jv) === jv.data
            @test parent(jv) == [1, 2, 3, 4, 5, 6]
            @test collect(Iterators.flatten(jv)) == parent(jv)
        end
        
        @testset "`reduce(vcat, ...)`" begin
            flattened = reduce(vcat, jv)
            @test flattened == [1, 2, 3, 4, 5, 6]
            @test flattened !== jv.data # should be a copy
            
            # test with empty JaggedVector
            jv_empty = JaggedVector{Int}(Int[], [1])
            @test reduce(vcat, jv_empty) == Int[]
        end
        
        @testset "convert to `Vector{Vector{T}}`" begin
            vv = convert(Vector{Vector{Int}}, jv)
            @test vv == [[1, 2], [3], [4, 5, 6]]
            
            # modify the converted result shouldn't affect original
            vv[1][1] = 999
            @test jv[1][1] == 1 # should be unchanged
            
            # test type conversion
            jv_float = JaggedVector([[1.0, 2.0], [3.0]])
            vv_float = convert(Vector{Vector{Float64}}, jv_float)
            @test vv_float isa Vector{Vector{Float64}}
            @test vv_float == [[1.0, 2.0], [3.0]]
        end
        
        @testset "`copy`" begin
            jv_copy = copy(jv)
            @test jv_copy.data == jv.data
            @test jv_copy.offsets == jv.offsets
            @test jv_copy.data !== jv.data # should be different objects
            @test jv_copy.offsets !== jv.offsets
            
            # modifying copy shouldn't affect original
            jv_copy.data[1] = 999
            @test jv.data[1] == 1
        end
        
        @testset "`IndexStyle`" begin
            @test Base.IndexStyle(JaggedVector{Int}) == IndexLinear()
        end

        @testset "`empty!`" begin
            jv_copy = copy(jv)
            empty!(jv)
            @test jv == JaggedVector{Int}()
        end
    end
    
    @testset "Edge cases" begin
        @testset "Empty JaggedVector" begin
            jv_empty = JaggedVector{Int}(Int[], [1])
            @test length(jv_empty) == 0
            @test isempty(jv_empty)
        end
        
        @testset "JaggedVector with empty sub-vectors" begin
            jv = JaggedVector([Int[1, 2], Int[], Int[3]])
            @test length(jv) == 3
            @test length(jv[1]) == 2
            @test length(jv[2]) == 0
            @test length(jv[3]) == 1
            @test jv[2] == Int[]
        end
        
        @testset "Single element JaggedVector" begin
            jv = JaggedVector([Int[1, 2, 3]])
            @test length(jv) == 1
            @test jv[1] == [1, 2, 3]
        end
        
        @testset "All empty sub-vectors" begin
            lengths = [0, 0, 0]
            jv = JaggedVector{Int}(lengths)
            @test length(jv) == 3
            @test all(isempty(jv[i]) for i in 1:3)
            @test length(jv.data) == 0
        end
    end
    
    @testset "Advanced features & implementation details" begin
        @testset "`@propagate_inbounds` behavior" begin
            # test that bounds checking works correctly
            jv = JaggedVector([[1, 2], [3, 4, 5]])
            
            # these should work
            @test jv[1] == [1, 2]
            @test jv[2] == [3, 4, 5]
            
            # this should throw bounds error
            @test_throws BoundsError jv[3]
        end
        
        @testset "Memory layout and efficiency" begin
            vs = [[1, 2, 3], [4, 5], [6]]
            jv = JaggedVector(vs)
            
            # test that data is stored contiguously
            @test jv.data == [1, 2, 3, 4, 5, 6]
            @test jv.offsets == [1, 4, 6, 7]
            
            # test that views share memory with data
            view1 = jv[1]
            jv.data[1] = 999
            @test view1[1] == 999 # view should reflect the change
        end
        
        @testset "Constructor from general iterables" begin
            # test with generators and other iterables
            gen = ([i, i+1] for i in 1:3)
            jv = JaggedVector{Int}(gen)
            
            @test length(jv) == 3
            @test jv[1] == [1, 2]
            @test jv[2] == [2, 3]
            @test jv[3] == [3, 4]
            
            # test automatic type inference
            gen2 = ([i, i+1] for i in 1:2)
            jv2 = JaggedVector(gen2)
            @test jv2 isa JaggedVector{Int}
        end
        
        @testset "Performance optimizations in convert" begin
            # test the optimized Vector{Vector{T}} conversion
            jv = JaggedVector([[1, 2], [3, 4, 5], [6]])
            vv = convert(Vector{Vector{Int}}, jv)
            
            # ensure it's a proper copy, not sharing memory
            vv[1][1] = 999
            @test jv[1][1] == 1
            
            # test the structure
            @test length(vv) == 3
            @test vv[1] == [999, 2] # modified
            @test vv[2] == [3, 4, 5]
            @test vv[3] == [6]
        end
    end

    @testset "`Float` to `Int` conversion" begin
        vs_float = [[1.0, 2.0], [3.0, 4.0]]
        jv_int = JaggedVector{Int}(vs_float)
        @test jv_int isa JaggedVector{Int}
        @test jv_int[1] == [1, 2]
        @test jv_int[2] == [3, 4]
    end
end