export getCellCenteredGrid, getNodalGrid, getEdgeGrids, getFaceGrids

function getCellCenteredGrid(M::OcTreeMesh)
    # Xc = getCellCenteredGrid(M::OcTreeMesh)

    i,j,k,bsz = find3(M.S)
    X = i+bsz/2
    Y = j+bsz/2
    Z = k+bsz/2

    X = (X .- 1) * M.h[1] .+ M.x0[1]
    Y = (Y .- 1) * M.h[2] .+ M.x0[2]
    Z = (Z .- 1) * M.h[3] .+ M.x0[3]

    return [X Y Z]
end

function getNodalGrid(M::OcTreeMesh)
    # X = getNodalGrid(M::OcTreeMesh)

    N      = getNodalNumbering(M.S)
    X,Y,Z  = find3(N)

    X = (X .- 1) * M.h[1] .+ M.x0[1]
    Y = (Y .- 1) * M.h[2] .+ M.x0[2]
    Z = (Z .- 1) * M.h[3] .+ M.x0[3]

    return [X Y Z]
end

function getEdgeGrids(M::OcTreeMesh)
    #  EX,EY,EZ = getEdgeGrids(M::OcTreeMesh)

    EX,EY,EZ = getEdgeSize(M)
    i,j,k,esz = find3(EX)
    EXX = (i+esz/2 .-1) * M.h[1] .+ M.x0[1]
    EXY = (j       .-1) * M.h[2] .+ M.x0[2]
    EXZ = (k       .-1) * M.h[3] .+ M.x0[3]


    i,j,k,esz = find3(EY)
    EYX = (i       .-1) * M.h[1] .+ M.x0[1]
    EYY = (j+esz/2 .-1) * M.h[2] .+ M.x0[2]
    EYZ = (k       .-1) * M.h[3] .+ M.x0[3]

    i,j,k,esz = find3(EZ)
    EZX = (i       .-1) * M.h[1] .+ M.x0[1]
    EZY = (j       .-1) * M.h[2] .+ M.x0[2]
    EZZ = (k+esz/2 .-1) * M.h[3] .+ M.x0[3]

    return [EXX EXY EXZ],[EYX EYY EYZ],[EZX EZY EZZ]
end

function getFaceGrids(M::OcTreeMesh)
    #  FX,FY,FZ = getFaceGrids(M::OcTreeMesh)

    FX,FY,FZ = getFaceSize(M)
    i,j,k,fsz = find3(FX)
    FXX = (i       .-1) * M.h[1] .+ M.x0[1]
    FXY = (j+fsz/2 .-1) * M.h[2] .+ M.x0[2]
    FXZ = (k+fsz/2 .-1) * M.h[3] .+ M.x0[3]


    i,j,k,fsz = find3(FY)
    FYX = (i+fsz/2 .-1) * M.h[1] .+ M.x0[1]
    FYY = (j       .-1) * M.h[2] .+ M.x0[2]
    FYZ = (k+fsz/2 .-1) * M.h[3] .+ M.x0[3]

    i,j,k,fsz = find3(FZ)
    FZX = (i+fsz/2 .-1) * M.h[1] .+ M.x0[1]
    FZY = (j+fsz/2 .-1) * M.h[2] .+ M.x0[2]
    FZZ = (k       .-1) * M.h[3] .+ M.x0[3]

    return [FXX FXY FXZ],[FYX FYY FYZ],[FZX FZY FZZ]
end
