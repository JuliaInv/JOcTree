
export splitCells

function splitCells(S::SparseArray3D, idx::Vector)
# Split cells idx into 2*2*2=8 smaller cells.

i,j,k,bsz = find3(S)
n = S.sz

S = splitCells( i,j,k,bsz, n, idx)

return S
end  # function splitCells

#-----------------------------

function splitCells( i,j,k,bsz, n::Tuple, idx::Vector)
# Split cells idx into 2*2*2=8 smaller cells.

old_mesh_size = length(i)
nidx = length(idx)  # # of cells to split

mesh_size = old_mesh_size + 7*nidx   # new mesh size
Tn2 = eltype(i)
Tn  = eltype(bsz)
ii  = Array{Tn2}(mesh_size)
jj  = Array{Tn2}(mesh_size)
kk  = Array{Tn2}(mesh_size)
vv  = Array{Tn}(mesh_size)

ii[1:old_mesh_size] = i
jj[1:old_mesh_size] = j
kk[1:old_mesh_size] = k
vv[1:old_mesh_size] = bsz

for icel = 1:nidx

   id = idx[icel]

   i1   = i[id]
   j1   = j[id]
   k1   = k[id]
   bsz12 = div(bsz[id], 2)

   vv[id] = bsz12  # original cell becomes smaller

   # Add 7 small cells
   id1 = old_mesh_size + (icel-1)*7 + 1
   id2 = id1 + 7-1

   i1b = i1 + bsz12
   j1b = j1 + bsz12
   k1b = k1 + bsz12

   ii[id1:id2] = vcat( i1b, i1,  i1b,  i1,  i1b, i1,  i1b )
   jj[id1:id2] = vcat( j1,  j1b, j1b,  j1,  j1,  j1b, j1b )
   kk[id1:id1+2] = k1
   kk[id1+3:id2] = k1b

   vv[id1:id2] = bsz12

end # icel

S = sparse3(ii,jj,kk,vv,n)

return S
end  # function splitCells
