export getEdgeMassMatrixAnisotropic

function getEdgeMassMatrixAnisotropic(S, h)
# P = getEdgeMassMatrixAnisotropic2(S, h)
# Compute a matrix P such that the edge function mass matrix
# for a general symmetric anisotropic coefficient can be computed by
#   M = P' * C * P
# Here,
#       | diag(cxx) diag(cxy) diag(cxz) |
#   C = | diag(cxy) diag(cyy) diag(cyz) |
#       | diag(cxz) diag(cyz) diag(czz) |
# expands the six components cxx, cyy, czz, cxy, cxz, cyz of the symmetric
# coefficient tensor for all numCell cells to a [3*numCells]^2 matrix.
#
# To illustrate the integration, consider the 2D case (QuadTree) and, in
# particular, the cell C1 of size 2 * h which has two right neighbors of
# of size h:
#
#            ex2
#      +-------------+------+------+
#      |q3         q4|      |      |
#      |             | ey3  |      |
#      |             |      |      |
#  ey1 |      C1     +------+------+
#      |             |      |      |
#      |             | ey2  |      |
#      |q1         q2|      |      |
#      +-------------+------+------+
#            ex1
#
# The integration of the mass matrix for the cell is performed in three
# steps:
# 1. Interpolate the edge function, which contains the x- and y-field
#    components ex1, ex2, ey1, ey2, ey3 staggered on the x- and y- edges,
#    to the four quadrature points qi (i = 1...4)
#      Ex(qi) = Px(ex1,ex2)
#      Ey(qi) = Py(ey1,ey2,ey3)
# 2. Compute the inner product of the collocated field, weighted by the
#    tensor C1
#      I1 = (Ex(qi), Ey(qi)) * C1 * (Ex(qi), Ey(qi))^T   (i = 1...4)
# 3. Integrate using numerical quadrature
#      M1 = (I1 + I2 + I3 + I4) / 4 * (2 * h)^2
#
# We use the following nearest neighbour interpolation
#   Ex(q1) = ex1, Ey(q1) = ey1
#   Ex(q2) = ex1, Ey(q2) = ey2
#   Ex(q3) = ex2, Ey(q3) = ey1
#   Ex(q4) = ex2, Ey(q4) = ey3
#
# The mass matrix has up to 9 non-zero entries per row/column.
#
# Note that this integration method neglects hanging edges. The resulting
# mass matrix M for an isotropic medium differs from the mass matrix Me
# which is obtained by the standard integration method by
#   Ae = getEdgeToCellCenterMatrix(S); # edge to cell center average
#   v  = prod(h) * nonzeros(S).^3;     # cell volume
#   Me = sdiag( (v .* c)' * Ae );
# where c is the isotropic coefficient for every cell. Let
#   C  = blkdiag(sdiag(c), sdiag(c), sdiag(c));
#   M  = P' * C * P;
# Then,
#   norm(M - Me) > 0
# unless the OcTree is degenerate (a regular mesh).

# C. Schwarzbach, April 2014

# edge numbering (implies location of edge)

EX,EY,EZ  = getEdgeNumbering(S)

i,j,k,bsz = find3(S)
mx,my,mz  = size(S)
nex = nnz(EX)
ney = nnz(EY)
nez = nnz(EZ)

# locate edge numbers for quadrature points

# x-edges
i0 = i;
i1 = floor(Integer,i + bsz / 2);
j0 = j;
j1 = j + bsz;
k0 = k;
k1 = k + bsz;

ex000 = getNodesFromIndices(EX.SV,[mx;my+1;mz+1;],i0,j0,k0) 
ex100 = getNodesFromIndices(EX.SV,[mx;my+1;mz+1;],i1,j0,k0) 
ex010 = getNodesFromIndices(EX.SV,[mx;my+1;mz+1;],i0,j1,k0) 
ex110 = getNodesFromIndices(EX.SV,[mx;my+1;mz+1;],i1,j1,k0) 
ex001 = getNodesFromIndices(EX.SV,[mx;my+1;mz+1;],i0,j0,k1) 
ex101 = getNodesFromIndices(EX.SV,[mx;my+1;mz+1;],i1,j0,k1) 
ex011 = getNodesFromIndices(EX.SV,[mx;my+1;mz+1;],i0,j1,k1) 
ex111 = getNodesFromIndices(EX.SV,[mx;my+1;mz+1;],i1,j1,k1) 

ex100 = merge(ex100, ex000);
ex110 = merge(ex110, ex010);
ex101 = merge(ex101, ex001);
ex111 = merge(ex111, ex011);

# y-edges
i0 = i;
i1 = i + bsz;
j0 = j;
j1 = floor(Integer,j + bsz / 2);
k0 = k;
k1 = k + bsz;

ey000 = getNodesFromIndices(EY.SV,[mx+1;my;mz+1;],i0,j0,k0)
ey100 = getNodesFromIndices(EY.SV,[mx+1;my;mz+1;],i1,j0,k0)
ey010 = getNodesFromIndices(EY.SV,[mx+1;my;mz+1;],i0,j1,k0)
ey110 = getNodesFromIndices(EY.SV,[mx+1;my;mz+1;],i1,j1,k0)
ey001 = getNodesFromIndices(EY.SV,[mx+1;my;mz+1;],i0,j0,k1)
ey101 = getNodesFromIndices(EY.SV,[mx+1;my;mz+1;],i1,j0,k1)
ey011 = getNodesFromIndices(EY.SV,[mx+1;my;mz+1;],i0,j1,k1)
ey111 = getNodesFromIndices(EY.SV,[mx+1;my;mz+1;],i1,j1,k1)

ey010 = merge(ey010, ey000);
ey110 = merge(ey110, ey100);
ey011 = merge(ey011, ey001);
ey111 = merge(ey111, ey101);

# z-edges
i0 = i;
i1 = i + bsz;
j0 = j;
j1 = j + bsz;
k0 = k;
k1 = floor(Integer,k + bsz / 2);

ez000 = getNodesFromIndices(EZ.SV,[mx+1;my+1;mz;],i0,j0,k0)
ez100 = getNodesFromIndices(EZ.SV,[mx+1;my+1;mz;],i1,j0,k0)
ez010 = getNodesFromIndices(EZ.SV,[mx+1;my+1;mz;],i0,j1,k0)
ez110 = getNodesFromIndices(EZ.SV,[mx+1;my+1;mz;],i1,j1,k0)
ez001 = getNodesFromIndices(EZ.SV,[mx+1;my+1;mz;],i0,j0,k1)
ez101 = getNodesFromIndices(EZ.SV,[mx+1;my+1;mz;],i1,j0,k1)
ez011 = getNodesFromIndices(EZ.SV,[mx+1;my+1;mz;],i0,j1,k1)
ez111 = getNodesFromIndices(EZ.SV,[mx+1;my+1;mz;],i1,j1,k1)
        
ez001 = merge(ez001, ez000);
ez101 = merge(ez101, ez100);
ez011 = merge(ez011, ez010);
ez111 = merge(ez111, ez110);

# cell indices

nc = length(bsz);
c  = [1:nc;]

# set values for each cell to square root of cell volume and quadrature
# weight; product of P' * C * P gives the correct scaling by cell volume
# and quadrature weight
uc = sqrt(nonzeros(S).^3 * prod(h) / 8);

# integrate using edge to cell corner interpolation
Px = sparse(c, ex000, uc, nc, nex); # Px1
Py = sparse(c, ey000, uc, nc, ney); # Py1
Pz = sparse(c, ez000, uc, nc, nez); # Pz1
P1 = blkdiag(Px, Py, Pz);

Px = sparse(c, ex100, uc, nc, nex); # Px2
Py = sparse(c, ey100, uc, nc, ney); # Py2
Pz = sparse(c, ez100, uc, nc, nez); # Pz2
P2 = blkdiag(Px, Py, Pz);

Px = sparse(c, ex010, uc, nc, nex); # Px3
Py = sparse(c, ey010, uc, nc, ney); # Py3
Pz = sparse(c, ez010, uc, nc, nez); # Pz3
P3 = blkdiag(Px, Py, Pz);

Px = sparse(c, ex110, uc, nc, nex); # Px4
Py = sparse(c, ey110, uc, nc, ney); # Py4
Pz = sparse(c, ez110, uc, nc, nez); # Pz4
P4 = blkdiag(Px, Py, Pz);

Px = sparse(c, ex001, uc, nc, nex); # Px5
Py = sparse(c, ey001, uc, nc, ney); # Py5
Pz = sparse(c, ez001, uc, nc, nez); # Pz5
P5 = blkdiag(Px, Py, Pz);

Px = sparse(c, ex101, uc, nc, nex); # Px6
Py = sparse(c, ey101, uc, nc, ney); # Py6
Pz = sparse(c, ez101, uc, nc, nez); # Pz6
P6 = blkdiag(Px, Py, Pz);

Px = sparse(c, ex011, uc, nc, nex); # Px7
Py = sparse(c, ey011, uc, nc, ney); # Py7
Pz = sparse(c, ez011, uc, nc, nez); # Pz7
P7 = blkdiag(Px, Py, Pz);

Px = sparse(c, ex111, uc, nc, nex); # Px8
Py = sparse(c, ey111, uc, nc, ney); # Py8
Pz = sparse(c, ez111, uc, nc, nez); # Pz8
P8 = blkdiag(Px, Py, Pz);

#P = [P1; P2; P3; P4; P5; P6; P7; P8;]
# Following line faster than line above because hcat + transpose is
# faster than vcat.
P = [P1' P2' P3' P4' P5' P6' P7' P8']' 
return P

end

function merge(a, b)
# copy entries of b into a for zero entries of a
c = copy(a)
idz    = c .== 0;
c[idz] = b[idz]

return c

end

function getNodesFromIndices(sv,mm,i0::Vector{Int},j0::Vector{Int},k0::Vector{Int})
	
	jj = sub2ind([mm[1],mm[2],mm[3]],i0,j0,k0);
	v = zeros(Int64,length(jj));
	for i=1:length(v); v[i] = sv[jj[i]]; end;
	
	return v
	
end