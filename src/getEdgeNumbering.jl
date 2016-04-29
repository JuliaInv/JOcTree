export getEdgeNumbering

function getEdgeNumbering(M::OcTreeMesh)
   if isempty(M.NEX) || isempty(M.NEY) || isempty(M.NEZ)
		M.NEX,M.NEY,M.NEZ = getEdgeNumbering(M.S)
	end
   return M.NEX,M.NEY,M.NEZ
end

function getEdgeNumbering(S)
# [EX,EY,EZ] = getEdgeNumbering(S)
# Numbering of the edges of the structure

m1,m2,m3 = S.sz
i,j,k,bsz = find3(S)


sizeEX      = (m1,m2+1,m3+1)
ii          = [  i      ;   i       ;   i       ;   i       ;]
jj          = [  j      ;   j+bsz   ;   j       ;   j+bsz   ;]
kk          = [  k      ;   k       ;   k+bsz   ;   k+bsz   ;]
ii,jj,kk  = ind2sub(sizeEX,sort(unique(sub2ind(sizeEX,ii,jj,kk))))  # make em unique
EX          = sparse3(ii,jj,kk,1:length(ii),[m1;m2+1;m3+1;])

sizeEY      = (m1+1,m2,m3+1)
ii          = [  i      ;   i+bsz   ;   i       ;   i+bsz   ;]
jj          = [  j      ;   j       ;   j       ;   j       ;]
kk          = [  k      ;   k       ;   k+bsz   ;   k+bsz   ;]
ii,jj,kk  = ind2sub(sizeEY,sort(unique(sub2ind(sizeEY,ii,jj,kk))))  # make em unique
EY          = sparse3(ii,jj,kk,1:length(ii),[m1+1;m2;m3+1;]);

sizeEZ      = (m1+1,m2+1,m3)
ii          = [  i       ;    i+bsz   ;   i       ;   i+bsz   ;]
jj          = [  j       ;    j       ;   j+bsz   ;   j+bsz   ;]
kk          = [  k       ;    k       ;   k       ;   k       ;]
ii,jj,kk  = ind2sub(sizeEZ,sort(unique(sub2ind(sizeEZ,ii,jj,kk))))  # make em unique
EZ          = sparse3(ii,jj,kk,1:length(ii),[m1+1;m2+1;m3;])

return EX, EY, EZ

end
