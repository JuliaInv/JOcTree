export getEdgeGrids

function getEdgeGrids(M::OcTreeMesh)
#  EX,EY,EZ = getEdgeGrids(M::OcTreeMesh)

EX,EY,EZ = getEdgeSize(M.S)
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
