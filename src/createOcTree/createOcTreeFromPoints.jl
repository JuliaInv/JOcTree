export createOcTreeFromPoints

function createOcTreeFromPoints(n,nb,Xi)
#S = createOcTreeFromPoints(n,X)


# initialize the coarsest
i,j,k = ndgrid(1:nb:n[1],1:nb:n[2],1:nb:n[3])
S = sparse3(vec(i),vec(j),vec(k),ones(Int,prod(size(i)))*nb,n)

while minimum(nonzeros(S)) > 1
    # find the closest cell
    X,Y,Z = getCellCenteredGrid(S)
    tau = zeros(nnz(S))
    for jj = 1:size(Xi,1)
       d = (X-Xi[jj,1]).^2 + (Y-Xi[jj,2]).^2 + (Z-Xi[jj,3]).^2
       i = indmin(d)
       tau[i] = 1
    end
    S = refineOcTree(S,tau,0.9)
    S = regularizeOcTree(S)
end

return S

end
