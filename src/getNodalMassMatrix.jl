export getNodalMassMatrix

function getNodalMassMatrix(M::OcTreeMesh,sigma::Vector)
    An = getNodalAverageMatrix(M)
    V = getVolume(M)
    Mn = spdiagm(An'*(V*sigma))
  return Mn
end