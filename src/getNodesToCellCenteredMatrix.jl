export getNodeToCellCenteredMatrix, getNodalAverageMatrix

function getNodalAverageMatrix(M::OcTreeMesh)
	if isempty(M.An)
		M.An, = getNodeToCellCenteredMatrix(M)
	end
	return M.An
end

function getNodeToCellCenteredMatrix(M::OcTreeMesh)

G = copy(getNodalGradientMatrix(M))
G.nzval = ones(nnz(G))
A = getEdgeToCellCenteredMatrix(M.S)[1]
A = A*G
D = spdiags(1./sum(A,2), [0], size(A,1.), size(A,1.))
A = D*A
return A

end
