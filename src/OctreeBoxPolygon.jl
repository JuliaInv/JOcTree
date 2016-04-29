
export OctreeBoxPolygon, OctreeBoxPolygonTopo, OctreePolygonBelowSurf


function OctreeBoxPolygon( S::SparseArray3D,
                           h::Vector{Float64},          # (3) underlying cell size
                           x0::Vector{Float64},         # corner coordinates
                           x::Vector{Float64}, y::Vector{Float64},  # polygon points
                           k1,k2,
                           cellsize)
# Sf( (x,y), k1:k2 ) = cellsize
# Set cells inside the polygon (x,y) to cellsize.

if  k1>k2 || k1<1 || k2>S.sz[3]
   error("k1>k2 ...")
end

const dx = h[1]
const dy = h[2]

min_x = minimum(x)
max_x = maximum(x)
min_y = minimum(y)
max_y = maximum(y)


while true
   i,j,k,bsz = find3(S)
   nns = nnz(S)
   tau = zeros(nns)

   for ic = 1:nns
      ii = i[ic]
      jj = j[ic]
      kk = k[ic]
      bb = bsz[ic]

      if bb <= cellsize
         continue  # cell is already small
      end


      if kk+bb-1 >= k1 && kk <= k2
         # X,Y points making up the cell
         xx0 = x0[1] + (ii-1)*dx
         yy0 = x0[2] + (jj-1)*dy
         xx1 = xx0 + bb*dx
         yy1 = yy0 + bb*dy

         poly_in_cell = ( min_x > xx0 && max_x < xx1 &&
                          max_y > yy0 && min_y < yy1 )  ||
                        ( min_y > yy0 && max_y < yy1 &&
                          max_x > xx0 && min_x < xx1 )

         if poly_in_cell
            # Cell gets split if the polygon is inside the cell.
   	      tau[ic] = 1  # large cell to split

         elseif insidePolygon(x,y, xx0,yy0) || insidePolygon(x,y, xx0,yy1) ||
                insidePolygon(x,y, xx1,yy0) || insidePolygon(x,y, xx1,yy1)
            # Cell gets split if any of the corners are inside the polygon.
   	      tau[ic] = 1  # large cell to split
         end

      end  # kk+bb-1 >= k1 && kk <= k2

   end  # ic


   if all(tau .== 0)
      break  # we are done
   end

   S = refineOcTree(S,tau,0.9)

end  # while

return S
end # function OctreeBoxPolygon

#-----------------------------------------------------------------

function OctreeBoxPolygonTopo( S::SparseArray3D,
                           h::Vector{Float64},          # (3) underlying cell size
                           x0::Vector{Float64},         # corner coordinates
                           x::Vector{Float64}, y::Vector{Float64},  # polygon points
                           itopo::Array{Int64,2},       # # of SURFACE cells
                           cellsize)
# Sf( (x,y), itopo ) = cellsize
# Set cells inside the polygon (x,y) to cellsize.

const dx = h[1]
const dy = h[2]

min_x = minimum(x)
max_x = maximum(x)
min_y = minimum(y)
max_y = maximum(y)


while true
   i,j,k,bsz = find3(S)
   nns = nnz(S)
   tau = zeros(nns)

   for ic = 1:nns
      ii = i[ic]
      jj = j[ic]
      kk = k[ic]
      bb = bsz[ic]

      if bb <= cellsize
      	continue  # cell is already small
      end


      k1 = minimum(itopo[ii:ii+bb-1, jj:jj+bb-1])
      k2 = maximum(itopo[ii:ii+bb-1, jj:jj+bb-1])

      if kk+bb-1 >= k1 && kk <= k2
         # X,Y points making up the cell
         xx0 = x0[1] + (ii-1)*dx
         yy0 = x0[2] + (jj-1)*dy
         xx1 = xx0 + bb*dx
         yy1 = yy0 + bb*dy

         poly_in_cell = ( min_x > xx0 && max_x < xx1 &&
                          max_y > yy0 && min_y < yy1 )  ||
                        ( min_y > yy0 && max_y < yy1 &&
                          max_x > xx0 && min_x < xx1 )

         if poly_in_cell
            # Cell gets split if the polygon is inside the cell.
            tau[ic] = 1  # large cell to split

         elseif insidePolygon(x,y, xx0,yy0) || insidePolygon(x,y, xx0,yy1) ||
                insidePolygon(x,y, xx1,yy0) || insidePolygon(x,y, xx1,yy1)
            # Cell gets split if any of the corners are inside the polygon.
            tau[ic] = 1  # large cell to split
         end

      end  # kk+bb-1 >= k1 && kk <= k2

   end  # ic


   if all(tau .== 0)
      break  # we are done
   end

   S = refineOcTree(S,tau,0.9)

end  # while

return S
end # function OctreeBoxPolygonTopo

#-----------------------------------------------------------------------------------

function OctreePolygonBelowSurf( S::SparseArray3D,
                           h::Vector{Float64},          # (3) underlying cell size
                           x0::Vector{Float64},         # corner coordinates
                           x::Vector{Float64}, y::Vector{Float64},  # polygon points
                           itopo::Array{Int64,2},       # # of SURFACE cells
                           depth::Float64,              # how far down to go
                           cellsize)
# Sf( (x,y), ndepth:itopo ) = cellsize
# Set cells inside the polygon (x,y) to cellsize.

dx = h[1]
dy = h[2]
dz = h[3]

ndepth = cld(depth, dz)  # how many fine cells down to go

min_x = minimum(x)
max_x = maximum(x)
min_y = minimum(y)
max_y = maximum(y)


while true
   i,j,k,bsz = find3(S)
   nns = nnz(S)
   tau = zeros(nns)

   for ic = 1:nns
      ii = i[ic]
      jj = j[ic]
      kk = k[ic]
      bb = bsz[ic]

      if bb <= cellsize
      	continue  # cell is already small
      end


      k1 = minimum(itopo[ii:ii+bb-1, jj:jj+bb-1]) - ndepth
      k2 = maximum(itopo[ii:ii+bb-1, jj:jj+bb-1])

      if kk+bb-1 >= k1 && kk <= k2
         # X,Y points making up the cell
         xx0 = x0[1] + (ii-1)*dx
         yy0 = x0[2] + (jj-1)*dy
         xx1 = xx0 + bb*dx
         yy1 = yy0 + bb*dy

         poly_in_cell = ( min_x > xx0 && max_x < xx1 &&
                          max_y > yy0 && min_y < yy1 )  ||
                        ( min_y > yy0 && max_y < yy1 &&
                          max_x > xx0 && min_x < xx1 )

         if poly_in_cell
            # Cell gets split if the polygon is inside the cell.
            tau[ic] = 1  # large cell to split

         elseif insidePolygon(x,y, xx0,yy0) || insidePolygon(x,y, xx0,yy1) ||
                insidePolygon(x,y, xx1,yy0) || insidePolygon(x,y, xx1,yy1)
            # Cell gets split if any of the corners are inside the polygon.
            tau[ic] = 1  # large cell to split
         end

      end  # kk+bb-1 >= k1 && kk <= k2

   end  # ic


   if all(tau .== 0)
      break  # we are done
   end

   S = refineOcTree(S,tau,0.9)

end  # while

return S
end # function OctreePolygonBelowSurf
