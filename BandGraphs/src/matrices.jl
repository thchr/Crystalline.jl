# ---------------------------------------------------------------------------------------- #

# build adjacency matrix as a block matrix from subgraphs
function assemble_adjacency(subgraphs, partitions)
    Nirs_total  = last(partitions[end].iridxs)
    Nirs_each    = [length(p.lgirs) for p in partitions]

    A = BlockArray{Int}(zeros(Int, Nirs_total, Nirs_total), Nirs_each, Nirs_each)
    for subgraph in subgraphs
        i, j = subgraph.p_max.kidx, subgraph.p_nonmax.kidx
        # fill in block (i,j) of the adjacency matrix
        A[Block(i), Block(j)] .= subgraph.A
        A[Block(j), Block(i)] .= subgraph.A' # symmetric matrix (undirected graph)
    end

    return A
end

# ---------------------------------------------------------------------------------------- #

function assemble_degree(subgraphs, partitions)
    # Diagonal matrix D, whose entries Dᵢᵢ are equal to the degree of the ith node (i.e.,
    # the number of edges incident upon it). We define the degree of node nᵢ as on p. 9 of
    # Vergniory's 2017 PRB, i.e. as `Pₖ * irdim(lgirₖᵢ)` with `Pₖ` denoting the number of
    # distinct k-manifolds connected to the k-point of the partition that nᵢ belongs to and
    # `lgirₖᵢ` denoting the irrep associated with nᵢ.
    # In practice, this is equal to the number of edges incident upon each node, if each
    # edge is counted with the dimensionality of the connecting non-maximal irrep.
    # These choices are e.g., motivated by the desire to have the rows of the Laplacian
    # matrix sum to zero.

    v = mapreduce(irdims, vcat, partitions)    # irrep dimension for each node
    manifold_multiplicity = mapreduce(vcat, partitions) do p
        c = count(subgraphs) do subgraph
            if p.maximal  # maximal partition
                subgraph.p_max == p
            else          # nonmaximal partition
                subgraph.p_nonmax == p
            end
        end
        fill(c, length(p))
    end
    return Diagonal(v .* manifold_multiplicity)
end

# ---------------------------------------------------------------------------------------- #

function assemble_laplacian(subgraphs, partitions)
    A = assemble_adjacency(subgraphs, partitions)
    D = assemble_degree(subgraphs, partitions)
    return D - A
end

# ---------------------------------------------------------------------------------------- #