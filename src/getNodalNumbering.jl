export getNodalNumbering

function getNodalNumbering(M::OcTreeMesh)
	return getNodalNumbering(M.S)
end


function getNodalNumbering(S)
# N = getEdgeNumbering(S)
# Numbering of the nodes of an OcTree structure

m1,m2,m3 = S.sz
i,j,k,bsz = find3(S); bsz = round(Int64,bsz)

nind = [i        j       k;
        i        j       k+bsz;
        i        j+bsz   k;
        i        j+bsz   k+bsz;
        i+bsz    j       k;
        i+bsz    j       k+bsz;
        i+bsz    j+bsz   k;
        i+bsz    j+bsz   k+bsz ]

Ni = sparse3(nind[:,1],nind[:,2],nind[:,3],nind[:,3],[m1+1,m2+1,m3+1])


i,j,k = find3(Ni)
N = sparse3(i,j,k,1:length(i), [m1+1,m2+1,m3+1]);

return N

end
