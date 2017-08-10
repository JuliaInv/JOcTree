export getEdgeToCellCenteredMatrix, getEdgeAverageMatrix

function getEdgeAverageMatrix(M::OcTreeMesh)
	if isempty(M.Ae)
		Ae,   = getEdgeToCellCenteredMatrix(M.S)
		M.Ae  = Ae
		M.Aet = Ae'
	end
	return M.Ae,M.Aet
end

function getEdgeToCellCenteredMatrix(S::SparseArray3D)
#  getEdgeToCellCenteredMatrix(S::SparseArray3D)
#

n   = S.sz;
nex = n + [0, 1, 1]
ney = n + [1, 0, 1]
nez = n + [1, 1, 0]

i,j,k,bsz = find3(S)
ex,ey,ez = getEdgeNumbering(S)

nx = nnz(ex)
ny = nnz(ey)
nz = nnz(ez)
nc = nnz(S)

Px1 = ex.SV[sub2ind(nex,i,j,k),1]
Px2 = ex.SV[sub2ind(nex,i,j+bsz,k),1]
Px3 = ex.SV[sub2ind(nex,i,j,k+bsz),1]
Px4 = ex.SV[sub2ind(nex,i,j+bsz,k+bsz),1]

Py1 = ey.SV[sub2ind(ney,i,j,k),1]
Py2 = ey.SV[sub2ind(ney,i+bsz,j,k),1]
Py3 = ey.SV[sub2ind(ney,i,j,k+bsz),1]
Py4 = ey.SV[sub2ind(ney,i+bsz,j,k+bsz),1]

Pz1 = ez.SV[sub2ind(nez,i,j,k),1]
Pz2 = ez.SV[sub2ind(nez,i+bsz,j,k),1]
Pz3 = ez.SV[sub2ind(nez,i,j+bsz,k),1]
Pz4 = ez.SV[sub2ind(nez,i+bsz,j+bsz,k),1]


sp1(Q) = sparse(1:Base.nnz(Q),Base.nonzeros(Q),ones(Base.nnz(Q)),nc,nx)
sp2(Q) = sparse(1:Base.nnz(Q),Base.nonzeros(Q),ones(Base.nnz(Q)),nc,ny)
sp3(Q) = sparse(1:Base.nnz(Q),Base.nonzeros(Q),ones(Base.nnz(Q)),nc,nz)

Ax = 1/4 * (sp1(Px1) + sp1(Px2) + sp1(Px3) + sp1(Px4))
Ay = 1/4 * (sp2(Py1) + sp2(Py2) + sp2(Py3) + sp2(Py4))
Az = 1/4 * (sp3(Pz1) + sp3(Pz2) + sp3(Pz3) + sp3(Pz4))

Av = (1/3)*[Ax Ay Az]

return Av, Ax, Ay, Az

end
