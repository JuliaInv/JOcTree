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
initializeOctree(n::Tuple{Integer,Integer,Integer}, bsz::Integer=0) =
                  initializeOctree(Int,Int,n,bsz)
function initializeOctree{Tn<:Integer,Tn2<:Integer}(::Type{Tn}, ::Type{Tn2},
                              n::Tuple{Integer,Integer,Integer}, bsz::Integer=0)

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
    b = fill(Tn(nb), length(i))
    i = convert(Vector{Tn2},vec(i))
    j = convert(Vector{Tn2},vec(j))
    k = convert(Vector{Tn2},vec(k))
    S = sparse3(i,j,k, b, (Tn2(n[1]),Tn2(n[2]),Tn2(n[3])))
    return S
end
