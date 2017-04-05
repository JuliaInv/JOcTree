export getNodalGradientRec, getNodalGradientMatrix

function getNodalGradientMatrix(M::OcTreeMeshFV)
# M.Grad = getNodalGradientMatrix(M::OcTreeMesh)
   if isempty(M.Grad)
      M.Grad = getNodalGradientRec(M)
   end
   return M.Grad
end

using MaxwellUtils.DiagTimesM

function getNodalGradientRec(M)
#function getNodalGradientRec(S,h)
# GRAD = getGradientMatrixRec(S,h)
# G : nodal -> edges


S = M.S; h = M.h
N           = getNodalNumbering(M)
#ENX,ENY,ENZ = getEdgeNumbering(M)
#EX,EY,EZ    = getEdgeSize(M)
EX,EY,EZ, ENX,ENY,ENZ = getEdgeSizeNumbering(M)

i,j,k,esz = find3(EX)  #; esz = round(Int64,esz)

ii = [ nonzeros(ENX)                ; nonzeros(ENX)                   ]
jj = [ N.SV[sub2ind(N.sz,i,j,k),1]  ; N.SV[sub2ind(N.sz ,i+esz,j,k),1]]
vv = [ -ones(length(i))             ; ones(length(i))                 ]

#ii = round(Int64,vec(ii)); jj = round(Int64,vec(full(jj))); vv = vec(full(vv))

G1 = sparse(ii, jj, vv, nnz(EX), nnz(N))

i,j,k,esz = find3(EY) #;  esz = round(Int64,esz)
ii = [ nonzeros(ENY)                ; nonzeros(ENY)                    ]
jj = [ N.SV[sub2ind(N.sz,i,j,k),1]  ; N.SV[sub2ind(N.sz ,i,j+esz,k),1] ]
vv = [ -ones(length(i))             ; ones(length(i))                  ]

#ii = round(Int64,vec(ii)); jj = round(Int64,vec(full(jj))); vv = vec(full(vv))

G2 = sparse(ii, jj, vv, nnz(EY), nnz(N))


i,j,k,esz = find3(EZ) #;  esz = round(Int64,esz)
ii = [ nonzeros(ENZ)              ; nonzeros(ENZ)                     ]
jj = [ N.SV[sub2ind(N.sz,i,j,k),1]; N.SV[sub2ind(N.sz  ,i,j,k+esz),1] ]
vv = [ -ones(size(i))             ; ones(size(i))                     ]

#ii = round(Int64,vec(ii)); jj = round(Int64,vec(full(jj))); vv = vec(full(vv))

G3 = sparse(ii, jj, vv, nnz(EZ), nnz(N));

G = [G1 ; G2 ; G3];


<<<<<<< HEAD
ESZi = 1. ./ vcat(nonzeros(EX)*h[1],
                  nonzeros(EY)*h[2],
                  nonzeros(EZ)*h[3] )

#GRAD = ESZi * G
GRAD = DiagTimesM(ESZi, G)
=======
ESZi = 1./[ nonzeros(EX)*h[1];
            nonzeros(EY)*h[2];
            nonzeros(EZ)*h[3] ]

G = unsafe_scale!(ESZi,G)
>>>>>>> 623eeedc401ce98db0adfc8ff8bec159758cdff3

return G

end

import Base.scale!

"""
unsafe_scale! performs A = spdiagm(b)*A
modifying A in place and avoiding overhead
of converting b into a sparse matrix.
It's unsafe because it doesn't check that
b and A are of compatible sizes.
"""
function unsafe_scale!{T,N}(b::Vector{T}, A::SparseMatrixCSC{T,N})
    m, n = size(A)
    Anzval = A.nzval
    Arowval = A.rowval
    for col = 1:n, p = A.colptr[col]:(A.colptr[col+1]-1)
        @inbounds Anzval[p] = Anzval[p] * b[Arowval[p]]
    end
    return A
end