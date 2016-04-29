export getCellCenteredGrid


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
