#This file contains functions for getting edge lengths, face sizes and cell volumes. 
# 
# Functions:
# 
# getEdgeSize --> returns 3 sparseArray3Ds with x,y,z edge sizes as multiples of base mesh size.
# getLength --> returns diagonal matrix of edge lengths
# getFaceSize --> returns 3 sparseArray3Ds with x,y,z face sizes as multiples of base mesh size.
#                 It doesn't have a corresponding function that retursn face area, as far as I can tell
#                 
# getVolume --> Returns diagonal matrix of cell volumes

export getVolume,getFaceSize,getEdgeSize,getLength

"""
Mesh.L = getLength(Mesh::OcTreeMesh) computes edge lengths l, returns spdiagm(l)
"""
function getLength(Mesh::OcTreeMesh)
    if isempty(Mesh.L)
        l1, l2, l3 = getEdgeSize(Mesh)
        l1 = nonzeros(l1)*Mesh.h[1]
        l2 = nonzeros(l2)*Mesh.h[2]
        l3 = nonzeros(l3)*Mesh.h[3]
        Mesh.L = spdiagm([l1;l2;l3])
    end
    return Mesh.L
end

"""
Mesh.V = getVolume(M::OcTreeMesh) returns diagonal matrix of cell volumes
"""
function getVolume(M::OcTreeMesh)
    if isempty(M.V)
    	i,j,k,bsz = find3(M.S)
    	h         = M.h
    	M.V 		 = spdiagm(bsz.^3*prod(h))
    end
    return M.V
end


function getEdgeSize(M::OcTreeMesh)
   if nnz(M.EX) != M.ne[1] || nnz(M.EY) != M.ne[2] || nnz(M.EZ) != M.ne[3]
		M.EX,M.EY,M.EZ = getEdgeSize(M.S)
	end
   return M.EX,M.EY,M.EZ
end


function getEdgeSize(S::SparseArray3D)

m1,m2,m3 = S.sz
I,J,K,BSZ = find3(S)

# allocate
EXi = Int64[]; EXj = Int64[]; EXk = Int64[]; 
EYi = Int64[]; EYj = Int64[]; EYk = Int64[]; 
EZi = Int64[]; EZj = Int64[]; EZk = Int64[]; 
EV  = Int64[];

for bsz = 2.^round(Int64,[0:log2(maximum(nonzeros(S)));])
	ind = find(BSZ.==bsz)
	i = I[ind]; j = J[ind]; k = K[ind]
	
	if ~isempty(i)
		EV =  [EV;fill(bsz,4*length(i));]
		
		# get X edges
		ii  = [  i      ;   i       ;   i       ;    i        ;]
		jj  = [  j      ;   j.+bsz  ;   j       ;   j.+bsz    ;]
		kk  = [  k      ;   k       ;   k.+bsz  ;   k.+bsz    ;]
		EXi = [EXi; ii;]
		EXj = [EXj; jj;]
		EXk = [EXk; kk;]
		
		# get Y edges
		ii    = [  i      ;   i.+bsz   ;   i       ;   i.+bsz   ;]
		jj    = [  j      ;   j        ;   j       ;   j        ;]
		kk    = kk
		EYi = [EYi; ii;]
		EYj = [EYj; jj;]
		EYk = [EYk; kk;]
		
		# get Z edges
		ii    = ii
		jj    = [  j       ;    j       ;   j.+bsz   ;   j.+bsz  ; ]
		kk    = [  k       ;    k       ;   k        ;   k       ; ]
		EZi = [EZi; ii;]
		EZj = [EZj; jj;]
		EZk = [EZk; kk;]
	end
end

EX  = sparse3(EXi,EXj,EXk,EV,[m1;m2+1;m3+1;],(x,y)->x)
EY  = sparse3(EYi,EYj,EYk,EV,[m1+1;m2;m3+1;],(x,y)->x)
EZ  = sparse3(EZi,EZj,EZk,EV,[m1+1;m2+1;m3;],(x,y)->x)


return EX, EY, EZ
end

function getFaceSize(M::OcTreeMesh)
   if nnz(M.FX) != M.nf[1] || nnz(M.FY) != M.nf[2] || nnz(M.FZ) != M.nf[3]
		M.FX,M.FY,M.FZ = getFaceSize(M.S)
	end
   return M.FX,M.FY,M.FZ
end

function getFaceSize(S)

m1,m2,m3 = S.sz

I,J,K,BSZ = find3(S)

FXi = Int64[]; FXj = Int64[]; FXk = Int64[]; FXv = Int64[];
FYi = Int64[]; FYj = Int64[]; FYk = Int64[]; FYv = Int64[];
FZi = Int64[]; FZj = Int64[]; FZk = Int64[]; FZv = Int64[];

for bsz = 2.^round(Int64,[0:log2(maximum(nonzeros(S)));])
    ind = find(BSZ.==bsz)
    i = I[ind]; j = J[ind]; k = K[ind]

    if ~isempty(i)
			v = fill(bsz,2*length(i)) 
			#indices of upper (i,j,k) and lower (i+bsz,j,k) faces
			ii = [ i ; i.+bsz ;]; jj = [j;j;]; kk = [k;k;];
			FXi = [FXi; ii;]
			FXj = [FXj; jj;]
			FXk = [FXk; kk;]
			FXv = [FXv; v;]
			
			# indices of left (i,j,k) and right (i,j+bsz,k) faces
			ii  = [ i ; i ;]; jj = [j ; j.+bsz; ]; kk = kk;
			FYi = [FYi; ii;]
			FYj = [FYj; jj;]
			FYk = [FYk; kk;]
			FYv = [FYv; v;];
			
			# indices of front (i,j,k) and back (i,j,k+bsz) faces
			ii = ii; jj = [j;j;]; kk =  [k ; k.+bsz;];
			FZi = [FZi; ii;]
			FZj = [FZj; jj;]
			FZk = [FZk; kk;]
			FZv = [FZv; v;];
    end
end
FX  = sparse3(FXi,FXj,FXk,FXv,[m1+1;m2;m3;],(x,y)->x)
FY  = sparse3(FYi,FYj,FYk,FYv,[m1;m2+1;m3;],(x,y)->x)
FZ  = sparse3(FZi,FZj,FZk,FZv,[m1;m2;m3+1;],(x,y)->x)

return FX, FY, FZ

end