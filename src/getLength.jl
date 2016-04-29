export getLength

function getLength(Mesh::OcTreeMesh)
# Mesh.L = getLength(Mesh::OcTreeMesh) computes edge lengths l, returns sdiag(l)
    if isempty(Mesh.L)
        l1, l2, l3 = getEdgeSize(Mesh)
        l1 = nonzeros(l1)*Mesh.h[1]
        l2 = nonzeros(l2)*Mesh.h[2]
        l3 = nonzeros(l3)*Mesh.h[3]
        Mesh.L   = blkdiag(blkdiag(spdiagm(l1,0),spdiagm(l2,0)),spdiagm(l3,0))
    end
return Mesh.L
end