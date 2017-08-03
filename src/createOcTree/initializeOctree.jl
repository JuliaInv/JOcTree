export initializeOctree

"""
    S = initializeOctree(n, bsz)

    Initialize octree mesh to the coarsest cells

    Input:

        n::Vector{Int64} - number of underlying cells
        bsz::Int64=0     - cell size (0 for largest)

    Output:

        S::SparseArray3D

"""
function initializeOctree( n::Vector{Int64}, bsz::Int64=0)

    if !(ispow2(n[1]) && ispow2(n[2]) && ispow2(n[3]))
        error("n must be power of 2.")
    end

    if bsz == 0
        nb = minimum(n)  # largest cell size
    else
        if !ispow2(bsz) || bsz < 0
            error("bsz must be power of 2.")
       elseif any( n .< bsz )
           error("bsz is too large.")
       end
       nb = bsz
    end

    i,j,k = ndgrid(1:nb:n[1], 1:nb:n[2], 1:nb:n[3])
    b = fill(nb, length(i))

    S = sparse3(vec(i),vec(j),vec(k), b, n)
    return S
end
