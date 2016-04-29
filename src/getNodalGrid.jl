export getNodalGrid

function getNodalGrid(M::OcTreeMesh)
# X = getNodalGrid(M::OcTreeMesh)

N      = getNodalNumbering(M.S)
X,Y,Z  = find3(N)

X = (X .- 1) * M.h[1] .+ M.x0[1]
Y = (Y .- 1) * M.h[2] .+ M.x0[2]
Z = (Z .- 1) * M.h[3] .+ M.x0[3]

return [X Y Z]
end
