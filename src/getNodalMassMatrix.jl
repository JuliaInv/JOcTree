export getNodalMassMatrix

function getNodalMassMatrix(M::OcTreeMesh,sigma::Vector)
    An = getNodalAverageMatrix(M)
    v = getVolumeVector(M)
    Mn = spdiagm(An'*(v.*sigma))
  return Mn
end
