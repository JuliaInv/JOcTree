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
T   = eltype(h)
Tn2 = eltype(S.SV.nzind)
N           = getNodalNumbering(M)
EX,EY,EZ, ENX,ENY,ENZ = getEdgeSizeNumbering(M)

i,j,k,esz = find3(EX)  #; esz = round(Int64,esz)

ii = [ nonzeros(ENX)              ; nonzeros(ENX)                     ]
jj = [ N.SV[sub2ind(N.sz,i,j,k)]    ; N.SV[sub2ind(N.sz ,convert(Vector{Tn2},i+esz),j,k)] ]
vv = [ -ones(T,length(i))             ; ones(T,length(i))                 ]

G1 = sparse(ii, jj, vv, nnz(EX), nnz(N))

i,j,k,esz = find3(EY) #;  esz = round(Int64,esz)
ii = [ nonzeros(ENY)              ; nonzeros(ENY)                      ]
jj = [ N.SV[sub2ind(N.sz,i,j,k)]    ; N.SV[sub2ind(N.sz ,i,convert(Vector{Tn2},j+esz),k)] ]
vv = [ -ones(T,length(i))             ; ones(T,length(i))                  ]

G2 = sparse(ii, jj, vv, nnz(EY), nnz(N))


i,j,k,esz = find3(EZ) #;  esz = round(Int64,esz)
ii = [ nonzeros(ENZ)              ; nonzeros(ENZ)                     ]
jj = [ N.SV[sub2ind(N.sz,i,j,k)]  ; N.SV[sub2ind(N.sz  ,i,j,convert(Vector{Tn2},k+esz))] ]
vv = [ -ones(T,size(i))             ; ones(T,size(i))                     ]

G3 = sparse(ii, jj, vv, nnz(EZ), nnz(N));

G = [G1 ; G2 ; G3];

ESZi = one(T) ./ vcat(nonzeros(EX)*h[1],
                  nonzeros(EY)*h[2],
                  nonzeros(EZ)*h[3] )

GRAD = DiagTimesM(ESZi, G)

return GRAD

end
