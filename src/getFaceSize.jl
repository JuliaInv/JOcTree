export getFaceSize

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
