export getVolume

function getVolume(M::OcTreeMesh)
	if isempty(M.V)
		i,j,k,bsz = find3(M.S)
		h         = M.h
		M.V 		 = sdiag(bsz.^3*prod(h))
	end
	return M.V
end
