export getNodeToCellCenteredMatrix, getNodalAverageMatrix

function getNodalAverageMatrix(M::OcTreeMesh)
	if isempty(M.An)
		M.An, = getNodeToCellCenteredMatrix(M)
	end
	return M.An
end

function getNodeToCellCenteredMatrix(M::OcTreeMesh)

Gtmp = getNodalGradientMatrix(M)
G    = SparseMatrixCSC(copy(Gtmp.m),copy(Gtmp.n),copy(Gtmp.colptr),
                       copy(Gtmp.rowval),ones(nnz(Gtmp)))
A = getEdgeToCellCenteredMatrix(M.S)[1]
A = A*G
D = spdiagm(vec(1./sum(A,2)), 0, size(A,1.), size(A,1.))
A = D*A
return A

end
