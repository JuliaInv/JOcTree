export getNodalMassMatrix

function getNodalMassMatrix(M::OcTreeMesh,sigma::Vector)
    An = getNodalAverageMatrix(M)
    v = getVolume(M)
    Mn = spdiagm(An'*(v.*sigma))
  return Mn
end
