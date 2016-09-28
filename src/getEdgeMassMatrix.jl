export getEdgeMassMatrix, getdEdgeMassMatrix, getEdgeMassMatrixNoWeight


function getEdgeMassMatrix(M::OcTreeMeshFV,sigma::Vector)
    # For octree meshes.
	 P = getEdgeAverageMatrix(M)
	 M = P'* kron(speye(24),sdiag(sigma)) *P
	 return M
end

function getEdgeMassMatrixNoWeight(M::OcTreeMeshFV,sigma::Vector)
    # For octree meshes.
	 P = getEdgeMassMatrixAnisotropic(M.S,M.h)
	 M = P'* kron(speye(24),sdiag(sigma)) *P
	 return M
end



function getEdgeMassMatrix(M::OcTreeMeshFV,sigma::Array{Float64,2})
    return getEdgeMassMatrix(M.S,M.h,sigma)
end


function getEdgeMassMatrix(S::SparseArray3D,h,sigma::Array{Float64,2})
	P = getEdgeMassMatrixAnisotropic(S,h)
	S11 = sdiag(sigma[:,1])
	S22 = sdiag(sigma[:,2])
	S33 = sdiag(sigma[:,3])
	S12 = sdiag(sigma[:,4])
	S13 = sdiag(sigma[:,5])
	S23 = sdiag(sigma[:,6])

	Sig = [
	 			S11 S12 S13;
				S12 S22 S23;
				S13 S23 S33;
			]
	M = P'* kron(speye(8),Sig) * P
	return M
end

#-----------------------------------------------------

function getdEdgeMassMatrix(M::OcTreeMeshFV,v::Vector)
   # Derivative
	 P = getEdgeAverageMatrix(M)
	 dM = P'* sdiag(P*v) *  kron(ones(24),speye(nnz(M.S)))
	return dM
end


# function getdEdgeMassMatrix(S::SparseArray3D,h,v::Vector)
# #  M = getdEdgeMassMatrix(S,h,sigma)
# 	P  = getEdgeMassMatrixAnisotropic(S,h)
# 	dM = P'* sdiag(P*v) *  kron(ones(24),speye(nnz(S)))
# 	return dM
# end


