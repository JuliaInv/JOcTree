export getFaceGrids

function getFaceGrids(M::OcTreeMesh)
#  FX,FY,FZ = getFaceGrids(M::OcTreeMesh)

FX,FY,FZ = getFaceSize(M.S)
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
