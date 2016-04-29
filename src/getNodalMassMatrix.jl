export getNodalMassMatrix

function getNodalMassMatrix(M::OcTreeMesh,sigma::Vector)
    An = getNodeToCellCenteredMatrix(M)
    V = getVolume(M)
    Mn = spdiagm(An'*(V*sigma))
  return Mn
end