export getFaceToCellCenteredMatrix, getFaceAverageMatrix

function getFaceAverageMatrix(M::OcTreeMesh)
	if isempty(M.Af)
		Tf = eltype(M.h)
		M.Af, = getFaceToCellCenteredMatrix(Tf,M.S)
	end
	return M.Af
end

function getFaceToCellCenteredMatrixAnisotropic(S::SparseArray3D)
#  M = getFaceMassMatrix(S)
#

n = S.sz;
nex = (n[1]+1,n[2],n[3])
ney = (n[1],n[2]+1,n[3])
nez = (n[1],n[2],n[3]+1)

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

Tn = eltype(S.SV.nzval)
Tn1 = one(Tn)
sp1(Q) = sparse(Tn1:Tn(Base.nnz(Q)),Tn(Base.nonzeros(Q)),ones(Base.nnz(Q)),nc,nx)
sp2(Q) = sparse(Tn1:Tn(Base.nnz(Q)),Tn(Base.nonzeros(Q)),ones(Base.nnz(Q)),nc,ny)
sp3(Q) = sparse(Tn1:Tn(Base.nnz(Q)),Tn(Base.nonzeros(Q)),ones(Base.nnz(Q)),nc,nz)

Ax = 0.5*(sp1(Px1) + sp1(Px2))
Ay = 0.5*(sp2(Py1) + sp2(Py2))
Az = 0.5*(sp3(Pz1) + sp3(Pz2))

Af = [Ax Ay Az]

return Af, Ax, Ay, Az

end

function getFaceToCellCenteredMatrix(::Type{Tf},S) where Tf <: Number
# [A, A1, A2, A3] = getFaceToCellCenteredMatrix(S)

A      = getDivergenceMatrixRec(S,ones(Tf,3))
FX,FY,FZ = getFaceSize(S)

HF = vcat(nonzeros(FX), nonzeros(FY), nonzeros(FZ))

fill!(A.nzval, one(Tf))

HF = HF.^2
A = MTimesDiag(A, HF)

NF1 = nnz(FX); NF2 = nnz(FY); NF3 = nnz(FZ);
A1 = A[:,1:NF1]; A2 = A[:,(1+NF1):(NF1+NF2)];  A3 = A[:,(1+NF1+NF2):end];

# make sure it sums to 1.
W1 = one(Tf) ./ vec(sum(A1,2))
W2 = one(Tf) ./ vec(sum(A2,2))
W3 = one(Tf) ./ vec(sum(A3,2))

#A1 = W1*A1; A2 = W2*A2; A3 = W3*A3;
A1 = DiagTimesM(W1, A1)
A2 = DiagTimesM(W2, A2)
A3 = DiagTimesM(W3, A3)

A = [A1  A2   A3]

return A, A1, A2, A3

end
