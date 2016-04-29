
export getInterfaceWeights

function getInterfaceWeights(
                     M::OcTreeMeshFV,
                     itopo::Array{Int64,2}, # # of SURFACE cells
                     surfweight::Vector{Float64},   # values to assign to the surface cells
                     C=[] )  # alpha values

if any( surfweight .<= 0 )
   error("surfweight .<= 0")
end

FX,FY,FZ = getFaceSize(M)

nXfaces = nnz(FX)
nYfaces = nnz(FY)
nZfaces = nnz(FZ)

if length(C) == 0
   alpha1 = 1.; alpha2 = 1.; alpha3 = 1.;
elseif length(C) == 3
   alpha1 = C[1]; alpha2 = C[2]; alpha3 = C[3]
else
   error("bad regparams in getInterfaceWeights")
end

wx = fill( float(alpha1), nXfaces)
wy = fill( float(alpha2), nYfaces)
wz = fill( float(alpha3), nZfaces)  # no weights on the Z faces

nlayers = length(surfweight)


# X faces
ii,jj,kk,fsz = find3(FX)

for i = 1:nXfaces

   if ii[i]==1 || ii[i]==M.n[1]+1
      continue  # mesh edge
   end

   itp = min( itopo[ii[i]-1,jj[i]],
              itopo[ii[i]  ,jj[i]] )

   heightbot = kk[i]
   height = heightbot + fsz[i] - 1
   layers_down = itp - height + 1

   if layers_down >= 1 && layers_down <= nlayers
      wx[i] *= surfweight[layers_down]
   end
end  # i


# Y faces
ii,jj,kk,fsz = find3(FY)

for i = 1:nYfaces

   if jj[i]==1 || jj[i]==M.n[2]+1
      continue  # mesh edge
   end

   itp = min( itopo[ii[i], jj[i]-1],
              itopo[ii[i], jj[i]  ] )

   heightbot = kk[i]
   height = heightbot + fsz[i] - 1
   layers_down = itp - height + 1

   if layers_down >= 1 && layers_down <= nlayers
      wy[i] *= surfweight[layers_down]
   end
end  # i


W = vcat( wx, wy, wz )
return W   
end # function getInterfaceWeights
