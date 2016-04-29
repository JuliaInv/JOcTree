export getMatricesFEM, getLocalElementMatrices, getMassMatrixFEM, getDiffMassMatrixFEM

function getMatricesFEM(M::OcTreeMeshFEM,sigma)
# K, M, Msig = getMatricesFEM(S,h,sigma)
# 
 
S = M.S
h = M.h

Nx,Ny,Nz  = getEdgeNumbering(M)
        
ii,jj,kk,bsz = find3(S)

nex = nnz(Nx)
ney = nnz(Ny)
nez = nnz(Nz)
nc  = nnz(S)
ne  = nex + ney + nez

Ke,Me,Ge = getLocalElementMatrices(h)

IE = zeros(nc,12);
# identify edges
IE[:,1]   = Nx.SV[sub2ind(size(Nx),ii,jj,kk)]
IE[:,2]   = Nx.SV[sub2ind(size(Nx),ii,jj+bsz,kk)]
IE[:,3]   = Nx.SV[sub2ind(size(Nx),ii,jj,kk+bsz)]
IE[:,4]   = Nx.SV[sub2ind(size(Nx),ii,jj+bsz,kk+bsz)]
IE[:,5]   = Ny.SV[sub2ind(size(Ny),ii,jj,kk)] + nex
IE[:,6]   = Ny.SV[sub2ind(size(Ny),ii+bsz,jj,kk)] + nex
IE[:,7]   = Ny.SV[sub2ind(size(Ny),ii,jj,kk+bsz)] + nex
IE[:,8]   = Ny.SV[sub2ind(size(Ny),ii+bsz,jj,kk+bsz)] + nex
IE[:,9]   = Nz.SV[sub2ind(size(Nz),ii,jj,kk)] + nex + ney
IE[:,10]  = Nz.SV[sub2ind(size(Nz),ii+bsz,jj,kk)] + nex + ney
IE[:,11]  = Nz.SV[sub2ind(size(Nz),ii,jj+bsz,kk)] + nex + ney
IE[:,12]  = Nz.SV[sub2ind(size(Nz),ii+bsz,jj+bsz,kk)] + nex + ney

IE = round(Int64,IE)
ii = vec((kron(ones(Int64,1,12),IE))')
jj = vec((kron(IE,ones(Int64,1,12)))')
kc = vec(kron(bsz,Ke[:]))
km = vec(kron(bsz.^3,Me[:]))
ks = vec(kron(sigma.*(bsz.^3),Me[:]))

K  = sparse(ii,jj,vec(full(kc)),ne,ne)
M  = sparse(ii,jj,vec(full(km)),ne,ne)

Msig = sparse(ii,jj,vec(full(ks)),ne,ne)


return K, M, Msig

end

function getMassMatrixFEM(M::OcTreeMeshFEM,sigma)
# K, M, Msig = getMatricesFEM(S,h,sigma)
# 
 
S = M.S
h = M.h

Nx,Ny,Nz  = getEdgeNumbering(M)
        
ii,jj,kk,bsz = find3(S)

nex = nnz(Nx)
ney = nnz(Ny)
nez = nnz(Nz)
nc  = nnz(S)
ne  = nex + ney + nez

Ke,Me,Ge = getLocalElementMatrices(h)

IE = zeros(nc,12);
# identify edges
IE[:,1]  = Nx.SV[sub2ind(size(Nx),ii,jj,kk)]
IE[:,2]  = Nx.SV[sub2ind(size(Nx),ii,jj+bsz,kk)]
IE[:,3]  = Nx.SV[sub2ind(size(Nx),ii,jj,kk+bsz)]
IE[:,4]  = Nx.SV[sub2ind(size(Nx),ii,jj+bsz,kk+bsz)]
IE[:,5]  = Ny.SV[sub2ind(size(Ny),ii,jj,kk)] + nex
IE[:,6]  = Ny.SV[sub2ind(size(Ny),ii+bsz,jj,kk)] + nex
IE[:,7]  = Ny.SV[sub2ind(size(Ny),ii,jj,kk+bsz)] + nex
IE[:,8]  = Ny.SV[sub2ind(size(Ny),ii+bsz,jj,kk+bsz)] + nex
IE[:,9]  = Nz.SV[sub2ind(size(Nz),ii,jj,kk)] + nex + ney
IE[:,10]  = Nz.SV[sub2ind(size(Nz),ii+bsz,jj,kk)] + nex + ney
IE[:,11]  = Nz.SV[sub2ind(size(Nz),ii,jj+bsz,kk)] + nex + ney
IE[:,12]  = Nz.SV[sub2ind(size(Nz),ii+bsz,jj+bsz,kk)] + nex + ney

IE = round(Int64,IE)
ii = vec((kron(ones(Int64,1,12),IE))')
jj = vec((kron(IE,ones(Int64,1,12)))')
ks = vec(kron(sigma.*(bsz.^3),Me[:]))

Msig = sparse(ii,jj,vec(full(ks)),ne,ne)

return Msig

end


function getDiffMassMatrixFEM(M::OcTreeMeshFEM,u::Vector)


	S = M.S; h = M.h
	
	Nx,Ny,Nz  = getEdgeNumbering(M)
        
	ii,jj,kk,bsz = find3(S)

	nex = nnz(Nx)
	ney = nnz(Ny)
	nez = nnz(Nz)
	nc  = nnz(S)
	ne  = nex + ney + nez

	Ke,Me,Ge = getLocalElementMatrices(h)

	IE = zeros(nc,12)

	# identify edges
	IE[:,1]  = Nx.SV[sub2ind(size(Nx),ii,jj,kk)]
	IE[:,2]  = Nx.SV[sub2ind(size(Nx),ii,jj+bsz,kk)]
	IE[:,3]  = Nx.SV[sub2ind(size(Nx),ii,jj,kk+bsz)]
	IE[:,4]  = Nx.SV[sub2ind(size(Nx),ii,jj+bsz,kk+bsz)]
	IE[:,5]  = Ny.SV[sub2ind(size(Ny),ii,jj,kk)] + nex
	IE[:,6]  = Ny.SV[sub2ind(size(Ny),ii+bsz,jj,kk)] + nex
	IE[:,7]  = Ny.SV[sub2ind(size(Ny),ii,jj,kk+bsz)] + nex
	IE[:,8]  = Ny.SV[sub2ind(size(Ny),ii+bsz,jj,kk+bsz)] + nex
	IE[:,9]  = Nz.SV[sub2ind(size(Nz),ii,jj,kk)] + nex + ney
	IE[:,10]  = Nz.SV[sub2ind(size(Nz),ii+bsz,jj,kk)] + nex + ney
	IE[:,11]  = Nz.SV[sub2ind(size(Nz),ii,jj+bsz,kk)] + nex + ney
	IE[:,12]  = Nz.SV[sub2ind(size(Nz),ii+bsz,jj+bsz,kk)] + nex + ney

   IE = round(Int64,IE)
	U  = diagm(sparse(bsz.^3))*u[IE]

	ii = vec(IE')
	jj = kron([1:nnz(S);],ones(Int64,12))
	vv = vec(Me*U')

	dMsig = sparse(ii,jj,vv,ne,nc)

	return dMsig
	
end

function getLocalElementMatrices(h)
# [K,M] = getElementMatrices(a,b,c)
# Compute stiffness and mass matrix for brick element of size a x b x c.
#
# The edge degrees of freedom are enumerated lexicographically as follows:
#
# x-edges:
#    +---4---+
#   /|      /|
#  +---3---+ |
#  | |     | |   z
#  | +--2--|-+   |  y
#  |/      |/    |/
#  +---1---+     +--- x
#
# y-edges:
#    +-------+
#   7|      8|
#  +-------+ |
#  | |     | |   z
#  | +-----|-+   |  y
#  |5      |6    |/
#  +-------+     +--- x
#
# z-edges:
#    +-------+
#   /|      /|
#  +11-----+12
#  | |     | |   z
#  9 +----10-+   |  y
#  |/      |/    |/
#  +-------+     +--- x
#


a = h[1]; b = h[2]; c = h[3]

# unit cube edge mass matrix times 36
Mei = [
   4  2  2  1
   2  4  1  2
   2  1  4  2
   1  2  2  4]

Me = kron(speye(3),Mei)

# unit cube face mass matrix times 6
Mfi = [
  2  1
  1  2]

Mf = kron(speye(3),Mfi)

# topological curl operator (edge integral to face integral degrees of
# freedom)
CURL = [
   0  0  0  0  1  0 -1  0 -1  0  1  0
   0  0  0  0  0  1  0 -1  0 -1  0  1
  -1  0  1  0  0  0  0  0  1 -1  0  0
   0 -1  0  1  0  0  0  0  0  0  1 -1
   1 -1  0  0 -1  1  0  0  0  0  0  0
   0  0  1 -1  0  0 -1  1  0  0  0  0]

CURL = sparse(CURL)

E = diagm(sparse([a a a a b b b b c c c c]))    # edge length
F = diagm(sparse(1./[b*c b*c a*c a*c a*b a*b])) # face area
V = a*b*c                           # cell volume

# scale by metric
CURL = F * (CURL * E)
Me   = Me * (V / 36.0)
Mf   = Mf * (V / 6.0)

# assemble
K = CURL' * Mf * CURL
M = Me

# Topological Gradient
GRAD = spzeros(12,8)
GRAD[1,[1;2;]]  = [-1 1]
GRAD[2,[3;4;]]  = [-1 1]
GRAD[3,[5;6;]]  = [-1 1]
GRAD[4,[7;8;]]  = [-1 1]
GRAD[5,[1;3;]]  = [-1 1]
GRAD[6,[2;4;]]  = [-1 1]
GRAD[7,[5;7;]]  = [-1 1]
GRAD[8,[6;8;]]  = [-1 1]
GRAD[9,[1;5;]]  = [-1 1]
GRAD[10,[2;6;]] = [-1 1]
GRAD[11,[3;7;]] = [-1 1]
GRAD[12,[4;8;]] = [-1 1]

GRAD = E\GRAD

return K,M,GRAD

end