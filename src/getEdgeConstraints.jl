export getEdgeConstraints

function getEdgeConstraints(M::OcTreeMesh)
	if isempty(M.Ne)
		T = eltype(M.h)
		if all(M.S.SV.nzval.==M.S.SV.nzval[1]) # uniform mesh
			N = typeof(M.nc)
			M.Ne = speye(T,N,sum(M.ne))
			M.Qe = speye(T,N,sum(M.ne))
            M.activeEdges = collect(N,1:sum(M.ne))
		else
			M.Ne,M.Qe, Ce, M.activeEdges = getEdgeConstraints(T,M.S)
		end
	end
	return M.Ne,M.Qe, M.activeEdges
end

function getEdgeConstraints(::Type{Tf},S::SparseArray3D) where Tf
# N,Q,C,p = getEdgeConstraints
# Eliminates the hanging edges
#   C ... constraint matrix
#   N ... null space matrix (interpolation: coarse to fine)
#   Q ... projection matrix (restriction: fine to coarse)
#                     which satisfy Q * N = I
#                     (note that Q ~= N')
#   p ............... lookup table for new edge enumeration:
#                     new face number = p(old face number)
#                     (returns zero for deleted faces)
# requires reg2 OcTree
#

i0,j0,k0,bsz = find3(S)
Tn2 = eltype(i0)
bsz = convert(Vector{Tn2},bsz)

i1 = i0 + div.(bsz, Tn2(2))
i2 = i0 + bsz
j1 = j0 + div.(bsz, Tn2(2))
j2 = j0 + bsz
k1 = k0 + div.(bsz, Tn2(2))
k2 = k0 + bsz

upper,lower,left,right,front,back = getNumberOfNeighbors(S)
nn = [upper; lower; left; right; front; back]
if ~all( (nn .== 0) .| (nn .== 1) .| (nn .== 4) )
  error("Implemented only for regularized OcTree meshes")
end
ENX,ENY,ENZ = getEdgeNumbering(S)

Tn  = eltype(ENX.SV.nzval)
ex1 = Tn[]; ex2 = Tn[]; ex3 = Tn[]; ex4 = Tn[]; ex5 = Tn[]; ex6 = Tn[];
ey1 = Tn[]; ey2 = Tn[]; ey3 = Tn[]; ey4 = Tn[]; ey5 = Tn[]; ey6 = Tn[];
ez1 = Tn[]; ez2 = Tn[]; ez3 = Tn[]; ez4 = Tn[]; ez5 = Tn[]; ez6 = Tn[];


#    ^ k z
#   /
#  /
# +-------> j y
# |
# |
# v x


# y-edges
#        +--ey5-----------+---ey6--------+
#       /                /              /
#    ez2               ez4            ez6
#     /                /              /
#    +--ey3-----------+------ey4-----+
#   /                /              /
#  ez1             ez3            ez5
# /                /              /
# +--ey1----------+---ey2--------+

I = (upper.==4)  # find  "bigger" cells

if any(I)
    append!(ey1, ENY.SV[sub2ind(size(ENY), i0[I], j0[I], k0[I])])
    append!(ey2, ENY.SV[sub2ind(size(ENY), i0[I], j1[I], k0[I])])
    append!(ey3, ENY.SV[sub2ind(size(ENY), i0[I], j0[I], k1[I])])
    append!(ey4, ENY.SV[sub2ind(size(ENY), i0[I], j1[I], k1[I])])
    append!(ey5, ENY.SV[sub2ind(size(ENY), i0[I], j0[I], k2[I])])
    append!(ey6, ENY.SV[sub2ind(size(ENY), i0[I], j1[I], k2[I])])

    append!(ez1, ENZ.SV[sub2ind(size(ENZ), i0[I], j0[I], k0[I])])
    append!(ez2, ENZ.SV[sub2ind(size(ENZ), i0[I], j0[I], k1[I])])
    append!(ez3, ENZ.SV[sub2ind(size(ENZ), i0[I], j1[I], k0[I])])
    append!(ez4, ENZ.SV[sub2ind(size(ENZ), i0[I], j1[I], k1[I])])
    append!(ez5, ENZ.SV[sub2ind(size(ENZ), i0[I], j2[I], k0[I])])
    append!(ez6, ENZ.SV[sub2ind(size(ENZ), i0[I], j2[I], k1[I])])
end

I = (lower.==4) # find  "bigger" cells

if any(I)
    append!( ey1 , ENY.SV[sub2ind(size(ENY), i2[I], j0[I], k0[I])])
    append!( ey2 , ENY.SV[sub2ind(size(ENY), i2[I], j1[I], k0[I])])
    append!( ey3 , ENY.SV[sub2ind(size(ENY), i2[I], j0[I], k1[I])])
    append!( ey4 , ENY.SV[sub2ind(size(ENY), i2[I], j1[I], k1[I])])
    append!( ey5 , ENY.SV[sub2ind(size(ENY), i2[I], j0[I], k2[I])])
    append!( ey6 , ENY.SV[sub2ind(size(ENY), i2[I], j1[I], k2[I])])

    append!( ez1 , ENZ.SV[sub2ind(size(ENZ), i2[I], j0[I], k0[I])])
    append!( ez2 , ENZ.SV[sub2ind(size(ENZ), i2[I], j0[I], k1[I])])
    append!( ez3 , ENZ.SV[sub2ind(size(ENZ), i2[I], j1[I], k0[I])])
    append!( ez4 , ENZ.SV[sub2ind(size(ENZ), i2[I], j1[I], k1[I])])
    append!( ez5 , ENZ.SV[sub2ind(size(ENZ), i2[I], j2[I], k0[I])])
    append!( ez6 , ENZ.SV[sub2ind(size(ENZ), i2[I], j2[I], k1[I])])
end

################################

I = (left.==4) # find  "bigger" cells

if any(I)
    append!(ex1, ENX.SV[sub2ind(size(ENX), i0[I], j0[I], k0[I])])
    append!(ex2, ENX.SV[sub2ind(size(ENX), i1[I], j0[I], k0[I])])
    append!(ex3, ENX.SV[sub2ind(size(ENX), i0[I], j0[I], k1[I])])
    append!(ex4, ENX.SV[sub2ind(size(ENX), i1[I], j0[I], k1[I])])
    append!(ex5, ENX.SV[sub2ind(size(ENX), i0[I], j0[I], k2[I])])
    append!(ex6, ENX.SV[sub2ind(size(ENX), i1[I], j0[I], k2[I])])

    append!(ez1, ENZ.SV[sub2ind(size(ENZ), i0[I], j0[I], k0[I])])
    append!(ez2, ENZ.SV[sub2ind(size(ENZ), i0[I], j0[I], k1[I])])
    append!(ez3, ENZ.SV[sub2ind(size(ENZ), i1[I], j0[I], k0[I])])
    append!(ez4, ENZ.SV[sub2ind(size(ENZ), i1[I], j0[I], k1[I])])
    append!(ez5, ENZ.SV[sub2ind(size(ENZ), i2[I], j0[I], k0[I])])
    append!(ez6, ENZ.SV[sub2ind(size(ENZ), i2[I], j0[I], k1[I])])
end

I = (right.==4) # find  "bigger" cells

if any(I)
    append!(ex1, ENX.SV[sub2ind(size(ENX), i0[I], j2[I], k0[I])])
    append!(ex2, ENX.SV[sub2ind(size(ENX), i1[I], j2[I], k0[I])])
    append!(ex3, ENX.SV[sub2ind(size(ENX), i0[I], j2[I], k1[I])])
    append!(ex4, ENX.SV[sub2ind(size(ENX), i1[I], j2[I], k1[I])])
    append!(ex5, ENX.SV[sub2ind(size(ENX), i0[I], j2[I], k2[I])])
    append!(ex6, ENX.SV[sub2ind(size(ENX), i1[I], j2[I], k2[I])])

    append!(ez1, ENZ.SV[sub2ind(size(ENZ), i0[I], j2[I], k0[I])])
    append!(ez2, ENZ.SV[sub2ind(size(ENZ), i0[I], j2[I], k1[I])])
    append!(ez3, ENZ.SV[sub2ind(size(ENZ), i1[I], j2[I], k0[I])])
    append!(ez4, ENZ.SV[sub2ind(size(ENZ), i1[I], j2[I], k1[I])])
    append!(ez5, ENZ.SV[sub2ind(size(ENZ), i2[I], j2[I], k0[I])])
    append!(ez6, ENZ.SV[sub2ind(size(ENZ), i2[I], j2[I], k1[I])])
end

################################

I = (front.==4) # find  "bigger" cells

if any(I)
    append!(ex1, ENX.SV[sub2ind(size(ENX), i0[I], j0[I], k0[I])])
    append!(ex2, ENX.SV[sub2ind(size(ENX), i1[I], j0[I], k0[I])])
    append!(ex3, ENX.SV[sub2ind(size(ENX), i0[I], j1[I], k0[I])])
    append!(ex4, ENX.SV[sub2ind(size(ENX), i1[I], j1[I], k0[I])])
    append!(ex5, ENX.SV[sub2ind(size(ENX), i0[I], j2[I], k0[I])])
    append!(ex6, ENX.SV[sub2ind(size(ENX), i1[I], j2[I], k0[I])])

    append!(ey1, ENY.SV[sub2ind(size(ENY), i0[I], j0[I], k0[I])])
    append!(ey2, ENY.SV[sub2ind(size(ENY), i0[I], j1[I], k0[I])])
    append!(ey3, ENY.SV[sub2ind(size(ENY), i1[I], j0[I] ,k0[I])])
    append!(ey4, ENY.SV[sub2ind(size(ENY), i1[I], j1[I], k0[I])])
    append!(ey5, ENY.SV[sub2ind(size(ENY), i2[I], j0[I], k0[I])])
    append!(ey6, ENY.SV[sub2ind(size(ENY), i2[I], j1[I], k0[I])])
end

I = (back.==4) # find  "bigger" cells

if any(I)
    append!(ex1, ENX.SV[sub2ind(size(ENX), i0[I], j0[I], k2[I])])
    append!(ex2, ENX.SV[sub2ind(size(ENX), i1[I], j0[I], k2[I])])
    append!(ex3, ENX.SV[sub2ind(size(ENX), i0[I], j1[I], k2[I])])
    append!(ex4, ENX.SV[sub2ind(size(ENX), i1[I], j1[I], k2[I])])
    append!(ex5, ENX.SV[sub2ind(size(ENX), i0[I], j2[I], k2[I])])
    append!(ex6, ENX.SV[sub2ind(size(ENX), i1[I], j2[I], k2[I])])

    append!(ey1, ENY.SV[sub2ind(size(ENY), i0[I], j0[I], k2[I])])
    append!(ey2, ENY.SV[sub2ind(size(ENY), i0[I], j1[I], k2[I])])
    append!(ey3, ENY.SV[sub2ind(size(ENY), i1[I], j0[I], k2[I])])
    append!(ey4, ENY.SV[sub2ind(size(ENY), i1[I], j1[I], k2[I])])
    append!(ey5, ENY.SV[sub2ind(size(ENY), i2[I], j0[I], k2[I])])
    append!(ey6, ENY.SV[sub2ind(size(ENY), i2[I], j1[I], k2[I])])
end

#Convert from n x 1 arrays to vectors
ex1 = vec(ex1); ex2 = vec(ex2); ex3 = vec(ex3)
ex4 = vec(ex4); ex5 = vec(ex5); ex6 = vec(ex6)
ey1 = vec(ey1); ey2 = vec(ey2); ey3 = vec(ey3)
ey4 = vec(ey4); ey5 = vec(ey5); ey6 = vec(ey6)
ez1 = vec(ez1); ez2 = vec(ez2); ez3 = vec(ez3)
ez4 = vec(ez4); ez5 = vec(ez5); ez6 = vec(ez6)

#e7 are the edges neither in e1, e2, e3, e4, e5 or e6
nx  = Tn(nnz(ENX)); ny = Tn(nnz(ENY));  nz = Tn(nnz(ENZ))
Tn1 = one(Tn)
ex7 = setdiff([Tn1:nx;], [ex1; ex2; ex3; ex4; ex5; ex6])
ey7 = setdiff([Tn1:ny;], [ey1; ey2; ey3; ey4; ey5; ey6])
ez7 = setdiff([Tn1:nz;], [ez1; ez2; ez3; ez4; ez5; ez6])

# look up table for new edge numbering
bx = falses(nx)
by = falses(ny)
bz = falses(nz)

bx[[ex1; ex5; ex7]] = true;
by[[ey1; ey5; ey7]] = true;
bz[[ez1; ez5; ez7]] = true;
b = [bx; by; bz];

mx = sum(bx);
my = sum(by);
mz = sum(bz);

px = zeros(Tn,nx)
py = zeros(Tn,ny)
pz = zeros(Tn,nz)
p  = zeros(Tn,nx+ny+nz)

px[bx] = 1:mx;
py[by] = 1:my;
pz[bz] = 1:mz;
p[b]   = 1:mx+my+mz;

###################################################################
##
## Build interpolation matrix
##
###################################################################
#        e1    e5    e2    e6    e3   e4
#        I     O     -I     O    O    O
#        O     I      O    -I    O    O
#       I/2    I/2    O     O   -I    O
#       I/2    I/2    O     O    O   -I
#
# Null space matrix
#
#    e1 I1
#    e5      I1
#    e2 I
#    e6      I
#    e3 I/2  I/2
#    e4 I/2  I/2


# X edges
# remove duplicates in [ex1; ex5]
ex15 = [ex1; ex5]
ex26 = [ex2; ex6]
ex15 = unique(ex15)
i15  = indexin(ex15,[ex1;ex5])
ex26 = ex26[i15]

k  = Tn(length(ex15))
m  = Tn(length(ex1))
n  = Tn(length(ex7))

# constraint matrix
i  = [Tn1:k; k+(Tn1:Tn(2*m)); k+(Tn1:Tn(2*m)); Tn1:k; k+(Tn1:Tn(2*m))]
j  = [ex15; ex1; ex1; ex5; ex5; ex26; ex3; ex4]
v  = [ones(Tf,k); 0.5*ones(Tf,4*m); -ones(Tf,k+2*m)]
Cx = sparse(i,j,v,k+2*m,nx)

# null space matrix
i  = [ex15; ex26; ex3; ex4; ex3; ex4; ex7]
j  = px[[ex15; ex15; ex1; ex1; ex5; ex5; ex7]]
v  = [ones(Tf,2*k); 0.5*ones(Tf,4*m); ones(Tf,n)]
Nx = sparse(i,j,v,nx,mx)

# projection matrix
i  = px[[ex15; ex7; ex15]]
j  = [ex15; ex7; ex26];
v  = [0.5*ones(Tf,k); ones(Tf,n); 0.5*ones(Tf,k)]
Qx = sparse(i,j,v,mx,nx)

# Y edges
# remove duplicates in [ey1; ey5]
ey15 = [ey1; ey5];
ey26 = [ey2; ey6];
ey15 = unique(ey15);
i15  = indexin(ey15,[ey1;ey5])
ey26 = ey26[i15]

k  = Tn(length(ey15))
m  = Tn(length(ey1))
n  = Tn(length(ey7))

# constraint matrix
i  = [Tn1:k; k+(Tn1:Tn(2*m)); k+(Tn1:Tn(2*m)); Tn1:k; k+(Tn1:Tn(2*m))]
j  = [ey15; ey1; ey1; ey5; ey5; ey26; ey3; ey4]
v  = [ones(Tf,k); 0.5*ones(Tf,4*m); -ones(Tf,k+2*m)]
Cy = sparse(i,j,v,k+2*m,ny)

# null space matrix
i  = [ey15; ey26; ey3; ey4; ey3; ey4; ey7]
j  = py[[ey15; ey15; ey1; ey1; ey5; ey5; ey7]]
v  = [ones(Tf,2*k); 0.5*ones(Tf,4*m); ones(Tf,n)]
Ny = sparse(i,j,v,ny,my)

# projection matrix
i  = py[[ey15; ey7; ey15]]
j  = [ey15; ey7; ey26]
v  = [0.5*ones(Tf,k); ones(Tf,n); 0.5*ones(Tf,k)]
Qy = sparse(i,j,v,my,ny)

## Z edges
# remove duplicates in [ez1; ez5]
ez15 = [ez1; ez5];
ez26 = [ez2; ez6];
ez15 = unique(ez15);
i15  = indexin(ez15,[ez1;ez5])
ez26 = ez26[i15]

k  = Tn(length(ez15))
m  = Tn(length(ez1))
n  = Tn(length(ez7))

# constraint matrix
i  = [Tn1:k; k+(Tn1:Tn(2*m)); k+(Tn1:Tn(2*m)); Tn1:k; k+(Tn1:Tn(2*m))]
j  = [ez15; ez1; ez1; ez5; ez5; ez26; ez3; ez4]
v  = [ones(Tf,k); 0.5*ones(Tf,4*m); -ones(Tf,k+2*m)]
Cz = sparse(i,j,v,k+2*m,nz)

# null space matrix
i  = [ez15; ez26; ez3; ez4; ez3; ez4; ez7]
j  = pz[[ez15; ez15; ez1; ez1; ez5; ez5; ez7]]
v  = [ones(Tf,2*k); 0.5*ones(Tf,4*m); ones(Tf,n)]
Nz = sparse(i,j,v,nz,mz)

# projection matrix
i  = pz[[ez15; ez7; ez15]]
j  = [ez15; ez7; ez26];
v  = [0.5*ones(Tf,k); ones(Tf,n); 0.5*ones(Tf,k)]
Qz = sparse(i,j,v,mz,nz)

## Put it all together

C = blkdiag(Cx,Cy,Cz)
N = blkdiag(Nx,Ny,Nz)
Q = blkdiag(Qx,Qy,Qz)

return N,Q,C,p

end
