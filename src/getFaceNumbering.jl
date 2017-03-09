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
ii = vcat(  i          , i+bsz )
jj = vcat(  j          , j     )
kk = vcat(  k          , k     )

sizeFX = (m1+1, m2, m3)
I = sort(unique(sub2ind(sizeFX, ii,jj,kk)))   ## create unique sorted linear indices
ii,jj,kk = ind2sub(sizeFX, I)         ## linear indices to nd indices
FX = sparse3(ii,jj,kk,1:length(ii), collect(sizeFX))


##     mark left     mark right
##      |                |
##      v                v
ii = vcat(  i          , i    )
jj = vcat(  j          , j+bsz)
kk = vcat(  k          , k    )
sizeFY = (m1, m2+1, m3)
I = sort(unique(sub2ind(sizeFY, ii,jj,kk)))   ## create unique sorted linear indices
ii,jj,kk = ind2sub(sizeFY, I)         ## linear indices to nd indices
FY = sparse3(ii,jj,kk,1:length(ii), collect(sizeFY))


##     mark front    mark back
##      |                |
##      v                v
ii = vcat(  i          , i     )
jj = vcat(  j          , j     )
kk = vcat(  k          , k+bsz )

sizeFZ = (m1, m2, m3+1)
I = sort(unique(sub2ind(sizeFZ, ii,jj,kk)))   ## create unique sorted linear indices
ii,jj,kk = ind2sub(sizeFZ, I)         ## linear indices to nd indices
FZ = sparse3(ii,jj,kk,1:length(ii), collect(sizeFZ))

return FX, FY, FZ

end  # function getFaceNumbering
