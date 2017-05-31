export getNodalGradientRec, getNodalGradientMatrix

function getNodalGradientMatrix(M::OcTreeMeshFV)
# M.Grad = getNodalGradientMatrix(M::OcTreeMesh)
        if isempty(M.Grad)
			M.Grad = getNodalGradientRec(M)
		end
        return M.Grad
end


function getNodalGradientRec(M)
#function getNodalGradientRec(S,h)
# GRAD = getGradientMatrixRec(S,h)
# G : nodal -> edges


S = M.S; h = M.h
N           = getNodalNumbering(M)
EX,EY,EZ, ENX,ENY,ENZ = getEdgeSizeNumbering(M)

i,j,k,esz = find3(EX)  #; esz = round(Int64,esz)

ii = [ nonzeros(ENX)              ; nonzeros(ENX)                     ]
jj = [ N.SV[sub2ind(N.sz,i,j,k),1]  ; N.SV[sub2ind(N.sz ,i+esz,j,k),1]]
vv = [ -ones(length(i))             ; ones(length(i))                 ]

G1 = sparse(ii, jj, vv, nnz(EX), nnz(N))

i,j,k,esz = find3(EY) #;  esz = round(Int64,esz)
ii = [ nonzeros(ENY)              ; nonzeros(ENY)                      ]
jj = [ N.SV[sub2ind(N.sz,i,j,k),1]  ; N.SV[sub2ind(N.sz ,i,j+esz,k),1] ]
vv = [ -ones(length(i))             ; ones(length(i))                  ]

G2 = sparse(ii, jj, vv, nnz(EY), nnz(N))


i,j,k,esz = find3(EZ) #;  esz = round(Int64,esz)
ii = [ nonzeros(ENZ)              ; nonzeros(ENZ)                     ]
jj = [ N.SV[sub2ind(N.sz,i,j,k),1]; N.SV[sub2ind(N.sz  ,i,j,k+esz),1] ]
vv = [ -ones(size(i))             ; ones(size(i))                     ]

G3 = sparse(ii, jj, vv, nnz(EZ), nnz(N));

G = [G1 ; G2 ; G3];

ESZi = 1. ./ vcat(nonzeros(EX)*h[1],
                  nonzeros(EY)*h[2],
                  nonzeros(EZ)*h[3] )

GRAD = DiagTimesM(ESZi, G)

return GRAD

end