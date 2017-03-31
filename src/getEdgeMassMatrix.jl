export getEdgeMassMatrix, getdEdgeMassMatrix, getEdgeMassMatrixNoWeight,
       DerivativeTimesVector, DerivativeTrTimesVector


function getEdgeMassMatrix(M::OcTreeMeshFV,sigma::Vector)
    if isempty(M.Pe)
        M.Pe = getEdgeMassMatrixAnisotropic(M.S,M.h)
    end
	M = M.Pe' * kron(speye(24),sdiag(sigma)) * M.Pe
	return M
end

function getEdgeMassMatrixNoWeight(M::OcTreeMeshFV,sigma::Vector)
    P = getEdgeMassMatrixAnisotropicNoWeight(M.S,M.h)    
	M = P'* kron(speye(24),sdiag(sigma)) *P
	return M
end



function getEdgeMassMatrix(M::OcTreeMeshFV,sigma::Array{Float64,2})
    return getEdgeMassMatrix(M.S,M.h,sigma)
end


function getEdgeMassMatrix(S::SparseArray3D,h,sigma::Array{Float64,2})
    if isempty(M.Pe)
        M.Pe = getEdgeMassMatrixAnisotropic(M.S,M.h)
    end
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
	M = M.Pe'* kron(speye(8),Sig) * M.Pe
	return M
end

#-----------------------------------------------------

function getdEdgeMassMatrix(M::OcTreeMeshFV,v::Vector)
    # Derivative
    if isempty(M.Pe)
        M.Pe = getEdgeMassMatrixAnisotropic(M.S,M.h)
    end
    dM = M.Pe'* sdiag(M.Pe*v) *  kron(ones(24),speye(nnz(M.S)))
    return dM
end


function DerivativeTimesVector(M::OcTreeMeshFV, v::Vector,
                               x::Vector)
   # Derivative (getdEdgeMassMatrix) times a vector(x)
   if isempty(M.Pe)
      M.Pe = getEdgeMassMatrixAnisotropic(M.S,M.h)
   end

  # dM = M.Pe'* sdiag(M.Pe*v) * kron(ones(24,1),speye(nnz(M.S)))

   pv = M.Pe * v

   dMw = M.Pe' * (pv .* repmat(x,24))
   
   return dMw
end  # function DerivativeTimesVector

function DerivativeTrTimesVector(M::OcTreeMeshFV, v::Vector,
                                 x::Vector)
   # Derivative (getdEdgeMassMatrix) transpose times a vector(x)
   if isempty(M.Pe)
      M.Pe = getEdgeMassMatrixAnisotropic(M.S,M.h)
   end

  # dM' = kron(ones(24,1),speye(nnz(M.S)))' * sdiag(M.Pe*v) * M.Pe

   pv = M.Pe * v
   dd = pv .* (M.Pe * conj(x)) 
   pv = [] 

   ns = nnz(M.S)
   dMTw = dd[1:ns]

   j = ns
   for i = 2:24
      @inbounds dMTw += dd[j+1 : j+ns]
      j += ns
   end  # i
   dMTw = conj(dMTw)

   return dMTw
end  # function DerivativeTrTimesVector

