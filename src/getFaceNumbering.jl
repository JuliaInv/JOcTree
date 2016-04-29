export getFaceNumbering

function getFaceNumbering(M::OcTreeMesh)
   if isempty(M.NFX) || isempty(M.NFY) || isempty(M.NFZ)
		M.NFX,M.NFY,M.NFZ = getFaceNumbering(M.S)
	end
   return M.NFX,M.NFY,M.NFZ
end



function getFaceNumbering(S)

m1,m2,m3 = S.sz
i,j,k,bsz = find3(S)

##     mark upper    mark lower
##      |                |
##      v                v
ii = [  i          ; i+bsz ;]
jj = [  j          ; j     ;]
kk = [  k          ; k     ;]

sizeFX = [m1+1 m2 m3]
I = sort(unique(sub2ind((sizeFX[1],sizeFX[2],sizeFX[3]),ii,jj,kk)))   ## create unique sorted linear indices
ii,jj,kk = ind2sub((sizeFX[1],sizeFX[2],sizeFX[3]),I)         ## linear indices to nd indices
FX = sparse3(ii,jj,kk,1:length(ii),[m1+1;m2;m3;])


##     mark left     mark right
##      |                |
##      v                v
ii = [  i          ; i    ;]
jj = [  j          ; j+bsz;]
kk = [  k          ; k    ;]
sizeFY = [m1 m2+1 m3]
I = sort(unique(sub2ind((sizeFY[1],sizeFY[2],sizeFY[3]),ii,jj,kk)))   ## create unique sorted linear indices
ii,jj,kk = ind2sub((sizeFY[1],sizeFY[2],sizeFY[3]),I)         ## linear indices to nd indices
FY = sparse3(ii,jj,kk,1:length(ii),[m1;m2+1;m3;])


##     mark front    mark back
##      |                |
##      v                v
ii = [  i          ; i     ;]
jj = [  j          ; j     ;]
kk = [  k          ; k+bsz ;]

sizeFZ = [m1 m2 m3+1]
I = sort(unique(sub2ind((sizeFZ[1],sizeFZ[2],sizeFZ[3]),ii,jj,kk)))   ## create unique sorted linear indices
ii,jj,kk = ind2sub((sizeFZ[1],sizeFZ[2],sizeFZ[3]),I)         ## linear indices to nd indices
FZ = sparse3(ii,jj,kk,1:length(ii),[m1;m2;m3+1;])

return FX, FY, FZ

end
