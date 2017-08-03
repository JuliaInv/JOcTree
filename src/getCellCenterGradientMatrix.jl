
export getCellCenterGradientMatrix

function getCellCenterGradientMatrix(Mesh::OcTreeMeshFV)

    CN = getCellNumbering(Mesh)
    FXN, FYN, FZN = getFaceNumbering(Mesh)
    i, j, k, bsz = find3(Mesh.S)

    upper, lower, left, right, front, back = getSizeOfNeighbors(Mesh.S)

    ######################################################################
    #### GRAD ON X-FACES
    ######################################################################
    ii = Int[]
    jj = Int[]
    vv = Float64[]

    # cells with the same size
    I = (upper .== 1)
    if any(I)
        append!(ii, sub2ind(size(FXN), i[I], j[I], k[I]))
        append!(jj, sub2ind(size(CN), i[I], j[I], k[I]))
        append!(vv, 1./bsz[I])
    end

    I = (lower .== 1)
    if any(I)
        append!(ii, sub2ind(size(FXN), i[I]+bsz[I], j[I], k[I]))
        append!(jj, sub2ind(size(CN), i[I], j[I], k[I]))
        append!(vv, -1./bsz[I])
    end

    # smaller cells
    I = (upper .== 0.5)
    if any(I)
        b2 = div.(bsz[I], 2)
        b3 = (1/6)*2./bsz[I]
        for i1 in [0, 1], i2 in [0, 1]
            append!(ii, sub2ind(size(FXN), i[I], j[I]+(b2*i1), k[I]+(b2*i2)))
            append!(jj, sub2ind(size(CN), i[I], j[I], k[I]))
            append!(vv, 4*b3)

            append!(ii, sub2ind(size(FXN), i[I], j[I]+(b2*i1), k[I]+(b2*i2)))
            append!(jj, sub2ind(size(CN), i[I]-b2, j[I], k[I]))
            append!(vv, -b3)

            append!(ii, sub2ind(size(FXN),i[I], j[I]+(b2*i1), k[I]+(b2*i2)))
            append!(jj, sub2ind(size(CN), i[I]-b2, j[I]+b2, k[I] ))
            append!(vv, -b3)

            append!(ii, sub2ind(size(FXN),i[I], j[I]+(b2*i1), k[I]+ (b2*i2)))
            append!(jj, sub2ind(size(CN), i[I]-b2, j[I], k[I]+b2 ))
            append!(vv, -b3)

            append!(ii, sub2ind(size(FXN),i[I], j[I]+(b2*i1), k[I]+(b2*i2)))
            append!(jj, sub2ind(size(CN), i[I]-b2, j[I]+b2, k[I]+b2))
            append!(vv, -b3)
        end
    end

    I = (lower .== 0.5)
    if any(I)
        b2 = div.(bsz[I],2)
        b3 = (1/6)./bsz[I]*2
        for i1 in [0, 1], i2 = [0, 1]
            append!(ii, sub2ind(size(FXN), i[I]+bsz[I], j[I]+(b2*i1), k[I]+(b2*i2)))
            append!(jj, sub2ind(size(CN), i[I], j[I], k[I]))
            append!(vv, -4*b3)

            append!(ii, sub2ind(size(FXN), i[I]+bsz[I], j[I]+(b2*i1), k[I]+(b2*i2)))
            append!(jj, sub2ind(size(CN), i[I]+bsz[I], j[I], k[I]))
            append!(vv, b3)

            append!(ii, sub2ind(size(FXN), i[I]+bsz[I], j[I]+(b2*i1), k[I]+(b2*i2)))
            append!(jj, sub2ind(size(CN), i[I]+bsz[I], j[I]+b2, k[I]))
            append!(vv, b3)

            append!(ii, sub2ind(size(FXN), i[I]+bsz[I], j[I]+(b2*i1), k[I]+ (b2*i2)))
            append!(jj, sub2ind(size(CN), i[I]+bsz[I], j[I], k[I]+b2))
            append!(vv, b3)

            append!(ii, sub2ind(size(FXN), i[I]+bsz[I], j[I]+(b2*i1), k[I]+(b2*i2)))
            append!(jj, sub2ind(size(CN), i[I]+bsz[I], j[I]+b2, k[I]+b2))
            append!(vv, b3)
        end
    end

    GX = sparse(FXN.SV[ii], CN.SV[jj], vv, nnz(FXN), nnz(CN))

    ######################################################################
    #### GRAD ON Y-FACES
    ######################################################################
    ii = Int[]
    jj = Int[]
    vv = Float64[]

    # cells with the same size
    I = (left .== 1)
    if any(I)
        append!(ii, sub2ind(size(FYN), i[I], j[I], k[I]))
        append!(jj, sub2ind(size(CN), i[I], j[I], k[I]))
        append!(vv, 1./bsz[I])
    end

    I = (right .== 1)
    if any(I)
        append!(ii, sub2ind(size(FYN), i[I], j[I]+bsz[I], k[I]))
        append!(jj, sub2ind(size(CN), i[I], j[I], k[I]))
        append!(vv, -1./bsz[I])
    end

    # smaller cells
    I = (left .== 0.5)
    if any(I)
        b2 = div.(bsz[I], 2)
        b3 = (1/6)*2./bsz[I]
        for i1 in [0, 1], i2 in [0, 1]
            append!(ii, sub2ind(size(FYN), i[I]+(b2*i1), j[I], k[I]+(b2*i2)))
            append!(jj, sub2ind(size(CN), i[I], j[I], k[I]))
            append!(vv, 4*b3)

            append!(ii, sub2ind(size(FYN), i[I]+(b2*i1) , j[I], k[I]+(b2*i2)))
            append!(jj, sub2ind(size(CN), i[I], j[I]-b2, k[I]))
            append!(vv, -b3)

            append!(ii, sub2ind(size(FYN), i[I]+(b2*i1), j[I], k[I]+(b2*i2)))
            append!(jj, sub2ind(size(CN), i[I]+b2, j[I]-b2, k[I] ))
            append!(vv, -b3)

            append!(ii, sub2ind(size(FYN), i[I]+(b2*i1), j[I], k[I]+ (b2*i2)))
            append!(jj, sub2ind(size(CN), i[I], j[I]-b2, k[I]+b2 ))
            append!(vv, -b3)

            append!(ii, sub2ind(size(FYN), i[I]+(b2*i1), j[I], k[I]+(b2*i2)))
            append!(jj, sub2ind(size(CN), i[I]+b2, j[I]-b2, k[I]+b2))
            append!(vv, -b3)
        end
    end

    I = (right .== 0.5)
    if any(I)
        b2 = div.(bsz[I], 2)
        b3 = (1/6)*2./bsz[I]
        for i1 in [0, 1], i2 in [0, 1]
            append!(ii, sub2ind(size(FYN), i[I]+(b2*i1), j[I]+bsz[I], k[I]+(b2*i2)))
            append!(jj, sub2ind(size(CN), i[I], j[I], k[I]))
            append!(vv, -4*b3)

            append!(ii, sub2ind(size(FYN), i[I] +(b2*i1), j[I]+bsz[I], k[I]+(b2*i2)))
            append!(jj, sub2ind(size(CN), i[I], j[I]+bsz[I], k[I]))
            append!(vv, b3)

            append!(ii, sub2ind(size(FYN), i[I]+(b2*i1), j[I]+bsz[I], k[I]+(b2*i2)))
            append!(jj, sub2ind(size(CN), i[I]+b2, j[I]+bsz[I], k[I] ))
            append!(vv, b3)

            append!(ii, sub2ind(size(FYN), i[I]+(b2*i1), j[I]+bsz[I], k[I]+ (b2*i2)))
            append!(jj, sub2ind(size(CN), i[I], j[I]+bsz[I], k[I]+b2 ))
            append!(vv, b3)

            append!(ii, sub2ind(size(FYN), i[I]+(b2*i1), j[I]+bsz[I], k[I]+(b2*i2)))
            append!(jj, sub2ind(size(CN), i[I]+b2, j[I]+bsz[I], k[I]+b2))
            append!(vv, b3)
        end
    end

    GY = sparse(FYN.SV[ii], CN.SV[jj], vv, nnz(FYN), nnz(CN))

    ######################################################################
    #### GRAD ON Z-FACES
    ######################################################################
    ii = Int[]
    jj = Int[]
    vv = Float64[]

    # cells with the same size
    I = (front .== 1)
    if any(I)
        append!(ii, sub2ind(size(FZN), i[I], j[I], k[I]))
        append!(jj, sub2ind(size(CN), i[I], j[I], k[I]))
        append!(vv, 1./bsz[I])
    end

    I = (back .== 1)
    if any(I)
        append!(ii, sub2ind(size(FZN), i[I], j[I], k[I]+bsz[I]))
        append!(jj, sub2ind(size(CN), i[I], j[I], k[I]))
        append!(vv, -1./bsz[I])
    end

    # smaller cells
    I = (front .== 0.5)
    if any(I)
        b2 = div.(bsz[I],2)
        b3 = (1/6)*2./bsz[I]
        for i1 in [0,1], i2 in [0,1]
            append!(ii, sub2ind(size(FZN), i[I]+b2*i1, j[I]+(b2*i2), k[I]))
            append!(jj, sub2ind(size(CN), i[I], j[I], k[I]))
            append!(vv, 4*b3)

            append!(ii, sub2ind(size(FZN), i[I]+b2*i1 , j[I]+(b2*i2), k[I]))
            append!(jj, sub2ind(size(CN), i[I], j[I], k[I]-b2))
            append!(vv, -b3)

            append!(ii, sub2ind(size(FZN), i[I]+b2*i1, j[I]+(b2*i2), k[I]))
            append!(jj, sub2ind(size(CN), i[I]+b2, j[I], k[I]-b2))
            append!(vv, -b3)

            append!(ii, sub2ind(size(FZN), i[I]+b2*i1, j[I]+ (b2*i2), k[I]))
            append!(jj, sub2ind(size(CN), i[I], j[I]+b2, k[I]-b2))
            append!(vv, -b3)

            append!(ii, sub2ind(size(FZN), i[I]+b2*i1, j[I]+(b2*i2), k[I]))
            append!(jj, sub2ind(size(CN), i[I]+b2, j[I]+b2, k[I]-b2))
            append!(vv, -b3)
        end
    end

    I = (back .== 0.5)
    if any(I)
        b2 = div.(bsz[I],2)
        b3 = (1/6)*2./bsz[I]
        for i1 in [0,1], i2 in [0,1]
            append!(ii, sub2ind(size(FZN), i[I]+b2*i1, j[I]+ (b2*i2), k[I]+bsz[I]))
            append!(jj, sub2ind(size(CN), i[I], j[I], k[I]))
            append!(vv, -4*b3)

            append!(ii, sub2ind(size(FZN), i[I] +b2*i1, j[I]+(b2*i2), k[I]+bsz[I]))
            append!(jj, sub2ind(size(CN), i[I], j[I], k[I]+bsz[I]))
            append!(vv, b3)

            append!(ii, sub2ind(size(FZN), i[I]+b2*i1, j[I]+(b2*i2), k[I]+bsz[I]))
            append!(jj, sub2ind(size(CN), i[I]+b2, j[I], k[I]+bsz[I]))
            append!(vv, b3)

            append!(ii, sub2ind(size(FZN), i[I]+b2*i1, j[I]+ (b2*i2), k[I]+bsz[I]))
            append!(jj, sub2ind(size(CN), i[I], j[I]+b2, k[I]+bsz[I]))
            append!(vv, b3)

            append!(ii, sub2ind(size(FZN), i[I]+b2*i1, j[I]+(b2*i2), k[I]+bsz[I]))
            append!(jj, sub2ind(size(CN), i[I]+b2, j[I]+b2, k[I]+bsz[I]))
            append!(vv, b3)
        end
    end

    GZ = sparse(FZN.SV[ii], CN.SV[jj], vv, nnz(FZN), nnz(CN))

    ##### put it all together
    GX = (1./Mesh.h[1])*GX
    GY = (1./Mesh.h[2])*GY
    GZ = (1./Mesh.h[3])*GZ

    nx = nnz(FXN)
    ny = nnz(FYN)
    nz = nnz(FZN)

    GRAD = spzeros( nx+ny+nz, nnz(CN))
    GRAD[1:nx, :] = GX
    GRAD[nx+1:nx+ny, :] = GY
    GRAD[nx+ny+1:nx+ny+nz, :] = GZ

    return GRAD
end
