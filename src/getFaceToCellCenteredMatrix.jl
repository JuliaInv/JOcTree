export getFaceToCellCenteredMatrix, getFaceAverageMatrix

function getFaceAverageMatrix(M::OcTreeMesh)
	if isempty(M.Af)
		M.Af, = getFaceToCellCenteredMatrix(M.S)
	end
	return M.Af
end

function getFaceToCellCenteredMatrixAnisotropic(S::SparseArray3D)
#  M = getFaceMassMatrix(S)
#

n = S.sz;
nex = n + [1, 0, 0]
ney = n + [0, 1, 0]
nez = n + [0, 0, 1]

i,j,k,bsz = find3(S)
ex,ey,ez = getFaceNumbering(S)

nx = nnz(ex)
ny = nnz(ey)
nz = nnz(ez)
nc = nnz(S)

Px1 = ex.SV[sub2ind(nex,i,j,k),1]
Px2 = ex.SV[sub2ind(nex,i+bsz,j,k),1]

Py1 = ey.SV[sub2ind(ney,i,j,k),1]
Py2 = ey.SV[sub2ind(ney,i,j+bsz,k),1]

Pz1 = ez.SV[sub2ind(nez,i,j,k),1]
Pz2 = ez.SV[sub2ind(nez,i,j,k+bsz),1]


sp1(Q) = sparse(1:Base.nnz(Q),Base.nonzeros(Q),ones(Base.nnz(Q)),nc,nx)
sp2(Q) = sparse(1:Base.nnz(Q),Base.nonzeros(Q),ones(Base.nnz(Q)),nc,ny)
sp3(Q) = sparse(1:Base.nnz(Q),Base.nonzeros(Q),ones(Base.nnz(Q)),nc,nz)

Ax = 0.5*(sp1(Px1) + sp1(Px2));
Ay = 0.5*(sp2(Py1) + sp2(Py2));
Az = 0.5*(sp3(Pz1) + sp3(Pz2));

Af = [Ax Ay Az]

return Af, Ax, Ay, Az

end

function getFaceToCellCenteredMatrix(S)
# [A, A1, A2, A3] = getFaceToCellCenteredMatrix(S)

A      = getDivergenceMatrixRec(S,[1,1,1])
FX,FY,FZ = getFaceSize(S)
HF       = sdiag([nonzeros(FX);nonzeros(FY);nonzeros(FZ)])

for ii = 1:nnz(A)
   A.nzval[ii] = 1.0
end

A  = A*(HF.^2);

NF1 = nnz(FX); NF2 = nnz(FY); NF3 = nnz(FZ);
A1 = A[:,1:NF1]; A2 = A[:,(1+NF1):(NF1+NF2)];  A3 = A[:,(1+NF1+NF2):end]; 

# make sure it sums to 1.
W1 = sdiag(1./vec(sum(A1,2)));
W2 = sdiag(1./vec(sum(A2,2)));
W3 = sdiag(1./vec(sum(A3,2)));

A1 = W1*A1; A2 = W2*A2; A3 = W3*A3;

A = [A1  A2   A3]

return A, A1, A2, A3

end
