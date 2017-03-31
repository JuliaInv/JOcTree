export getEdgeIntegralOfPolygonalChain

using Base.BLAS

function getEdgeIntegralOfPolygonalChain(mesh::OcTreeMesh, polygon::Array{Float64,2};
                                         normalize=false)
# s = getEdgeIntegralPolygonalChain(mesh,polygon)
# s = getEdgeIntegralPolygonalChain(mesh,polygon,normalize)
#
# Compute the integral of a piecewise linear edge grid basis projected onto
# the edges of a polygonal chain. This function can be used to evaluate the
# source term of a line current carried by the polygonal chain.
#
# The piecewise linear edge grid basis consists of the functions
#   phix_i (i = 1 ... numXEdges)
#   phiy_i (i = 1 ... numYEdges)
#   phiz_i (i = 1 ... numZEdges)
# where phix_i, phiy_i, phiz_i = 1 at the i-th x/y/z-directed edge
#   and phix_i, phiy_i, phiz_i = 0 at all other edges.
#
# INPUT
# mesh ...... OcTree mesh
# polygon ... vertices of polygonal chain as numVertices x 3 array
# normalize . divide line integral by length of integration path (boolean, optional)
# 
# OUTPUT
# s ......... source vector
# 
# For a closed current loop, specify polygon such that
#   polygon[1,:] == polygon[end,:]
#
# NOTE: The edge integration ignores hanging edges. If the polygon traverses cells
#       with hanging edges, it is imperative that hanging edges are eliminated
#       (getEdgeConstraints).

# Christoph Schwarzbach, January 2016

# OcTree mesh
S        = mesh.S
hx       = mesh.h[1]
hy       = mesh.h[2]
hz       = mesh.h[3]
x0       = mesh.x0[1]
y0       = mesh.x0[2]
z0       = mesh.x0[3]
nx,ny,nz = mesh.n

EX, EY, EZ = getEdgeNumbering(mesh)

# number of line segments
np = size(polygon, 1) - 1

# convert line segments to OcTree index coordinates
px = (polygon[:,1] .- x0) ./ hx .+ 1.0
py = (polygon[:,2] .- y0) ./ hy .+ 1.0
pz = (polygon[:,3] .- z0) ./ hz .+ 1.0

# check that all polygon vertices are inside the mesh
j = (px .< 1.0) | (px .> nx + 1.0) |
    (py .< 1.0) | (py .> ny + 1.0) |
    (pz .< 1.0) | (pz .> nz + 1.0)
n = sum(j)
if n > 0
   if n == 1
      msg = @sprintf("%d polygon vertex is outside the mesh\n", n)
   else
      msg = @sprintf("%d polygon vertices are outside the mesh\n", n)
   end
   for i = 1:np+1
      if j[i]
         msg = @sprintf("%s       vertex #%d: (%.3f,%.3f,%.3f)\n",
            msg, i, polygon[i,1], polygon[i,2], polygon[i,3])
       end
   end
   error(msg)
end

# check that all polygon vertices are contained in interior cells of mesh
j = (px .< 2.0) | (px .> nx) |
    (py .< 2.0) | (py .> ny) |
    (pz .< 2.0) | (pz .> nz)
n = sum(j)
if n > 0
   if n == 1
      msg = @sprintf("%d polygon vertex is too close to the boundary\n", n)
   else
      msg = @sprintf("%d polygon vertices are too close to the boundary\n", n)
   end
   for i = 1:np+1
      if j[i]
         msg = @sprintf("%s       vertex #%d: (%.3f,%.3f,%.3f)\n",
            msg, i, polygon[i,1], polygon[i,2], polygon[i,3])
       end
   end
   error(msg)
end


# unit cube to actual cell size
scaleX = hx
scaleY = hy
scaleZ = hz

if normalize
   scale = getNormalizeVal(polygon)
   scaleX /= scale
   scaleY /= scale
   scaleZ /= scale
end  # normalize


# fine mesh nodal grid
x = float([1:nx+1;])
y = float([1:ny+1;])
z = float([1:nz+1;])

# allocate space for source vector
nnX = nnz(EX)
nnY = nnz(EY)
nnZ = nnz(EZ)
ss = spzeros(nnX+nnY+nnZ,1)

# allocate local variables
sxloc = zeros(4)
syloc = zeros(4)
szloc = zeros(4)
kx    = zeros(Int64,4)
ky    = zeros(Int64,4)
kz    = zeros(Int64,4)
tol   = 0.0

# integrate each line segment
for ip = 1:np
   
   # start and end vertices
   ax = px[ip]
   ay = py[ip]
   az = pz[ip]
   bx = px[ip+1]
   by = py[ip+1]
   bz = pz[ip+1]
   
   # find intersection with mesh planes
   dx  = bx - ax
   dy  = by - ay
   dz  = bz - az
   d   = sqrt(dx*dx + dy*dy + dz*dz)
   tol = d * eps()
   # x-planes
   if bx > ax + tol
      tx = ([ceil(ax):floor(bx);] .- ax) ./ dx
   elseif ax > bx + tol
      tx = ([floor(ax):-1.0:ceil(bx);] .- ax) ./ dx
   else
      tx = zeros(0)
   end
   # y-planes
   if by > ay + tol
      ty = ([ceil(ay):floor(by);] .- ay) ./ dy
   elseif ay > by + tol
      ty = ([floor(ay):-1.0:ceil(by);] .- ay) ./ dy
   else
      ty = zeros(0)
   end
   # z-planes
   if bz > az + tol
      tz = ([ceil(az):floor(bz);] .- az) ./ dz
   elseif az > bz + tol
      tz = ([floor(az):-1.0:ceil(bz);] .- az) ./ dz
   else
      tz = zeros(0)
   end
   t = unique(sort([0.0; tx; ty; tz; 1.0;]))
   
   # number of cells that the current line segment crosses
   nq = length(t) - 1
   
   # integrate by cells
   for iq = 1:nq
      
      # center of line segment
      tc = 0.5 * (t[iq] + t[iq+1])
      
      # locate cell id
      ix = floor(Integer,ax + tc * dx)
      iy = floor(Integer,ay + tc * dy)
      iz = floor(Integer,az + tc * dz)
      ix,iy,iz,bsz = findBlocks(S,ix,iy,iz)
      
      # integration limits in local coordinates
      axloc = ax + t[iq]   * dx - x[ix]
      ayloc = ay + t[iq]   * dy - y[iy]
      azloc = az + t[iq]   * dz - z[iz]
      bxloc = ax + t[iq+1] * dx - x[ix]
      byloc = ay + t[iq+1] * dy - y[iy]
      bzloc = az + t[iq+1] * dz - z[iz]
      
      # basis functions are defined on cube of size bsz^3
      b = bsz * 1.0
      
      # integrate
      getStraightLineCurrentIntegral!(
         b, b, b, axloc, ayloc, azloc, bxloc, byloc, bzloc,
         sxloc, syloc, szloc)
      
      # find edge numbers
      jx = ix + bsz
      jy = iy + bsz
      jz = iz + bsz
      kx[1] = EX[ix,iy,iz]
      kx[2] = EX[ix,jy,iz]
      kx[3] = EX[ix,iy,jz]
      kx[4] = EX[ix,jy,jz]
      ky[1] = EY[ix,iy,iz]
      ky[2] = ky[1] + 1  # = EY[jx,iy,iz]
      ky[3] = EY[ix,iy,jz]
      ky[4] = ky[3] + 1  # = EY[jx,iy,jz]
      kz[1] = EZ[ix,iy,iz]
      kz[2] = kz[1] + 1  # = EZ[jx,iy,iz]
      kz[3] = EZ[ix,jy,iz]
      kz[4] = kz[3] + 1  # = EZ[jx,jy,iz]
      
      # add to source vector
      for ii = 1:4
         ss[kx[ii]]         += sxloc[ii] * scaleX
         ss[ky[ii]+nnX]     += syloc[ii] * scaleY
         ss[kz[ii]+nnX+nnY] += szloc[ii] * scaleZ
      end  # ii
      
   end  # iq
end  # ip


return ss
end # function getEdgeIntegralOfPolygonalChain


function getStraightLineCurrentIntegral!(hx,hy,hz,ax,ay,az,bx,by,bz,sx,sy,sz)
# s = getStraightLineCurrentIntegral(hx,hy,hz,ax,ay,az,bx,by,bz)
# Compute integral round(Int64,W . J dx^3) in brick of size hx x hy x hz
# where W denotes the 12 local bilinear edge basis functions
# and where J prescribes a unit line current
# between points (ax,ay,az) and (bx,by,bz).

# length of line segment
lx = bx - ax
ly = by - ay
lz = bz - az
l  = sqrt(lx*lx + ly*ly + lz*lz);

if l == 0.0
   fill!(sx, 0.0)
   fill!(sy, 0.0)
   fill!(sz, 0.0)
   return sx, sy, sz
end

# integration using Simpson's rule: need end points and mid point of interval
x1 = ax / hx
x2 = 0.5 * (ax + bx) / hx
x3 = bx / hx

y1 = ay / hy
y2 = 0.5 * (ay + by) / hy
y3 = by / hy

z1 = az / hz
z2 = 0.5 * (az + bz) / hz
z3 = bz / hz

u1 = 1.0 - x1
u2 = 1.0 - x2
u3 = 1.0 - x3

v1 = 1.0 - y1
v2 = 1.0 - y2
v3 = 1.0 - y3

w1 = 1.0 - z1
w2 = 1.0 - z2
w3 = 1.0 - z3

# integrate edge basis functions
sx[1] = v1 * w1 + 4.0 * v2 * w2 + v3 * w3
sx[2] = y1 * w1 + 4.0 * y2 * w2 + y3 * w3
sx[3] = v1 * z1 + 4.0 * v2 * z2 + v3 * z3
sx[4] = y1 * z1 + 4.0 * y2 * z2 + y3 * z3

sy[1] = u1 * w1 + 4.0 * u2 * w2 + u3 * w3
sy[2] = x1 * w1 + 4.0 * x2 * w2 + x3 * w3
sy[3] = u1 * z1 + 4.0 * u2 * z2 + u3 * z3
sy[4] = x1 * z1 + 4.0 * x2 * z2 + x3 * z3

sz[1] = u1 * v1 + 4.0 * u2 * v2 + u3 * v3
sz[2] = x1 * v1 + 4.0 * x2 * v2 + x3 * v3
sz[3] = u1 * y1 + 4.0 * u2 * y2 + u3 * y3
sz[4] = x1 * y1 + 4.0 * x2 * y2 + x3 * y3

# sx *= (lx / 6.0)
# sy *= (ly / 6.0)
# sz *= (lz / 6.0)

BLAS.scal!(4, lx / 6.0, sx, 1)
BLAS.scal!(4, ly / 6.0, sy, 1)
BLAS.scal!(4, lz / 6.0, sz, 1)

return  #sx, sy, sz

end  # function getStraightLineCurrentIntegral!

#-----------------------------------------------------

function getNormalizeVal( polygon::Array{Float64,2} )
# number of line segments
np = size(polygon, 1) - 1
   
   if all(polygon[1,:] .== polygon[np+1,:])
      
      # closed polygon: divide by enclosed area
      a  = 0.0
      px = polygon[2:np,1] .- polygon[1,1]
      py = polygon[2:np,2] .- polygon[1,2]
      pz = polygon[2:np,3] .- polygon[1,3]
      for ip = 1:np-2
         cx = py[ip] * pz[ip+1] - pz[ip] * py[ip+1]
         cy = pz[ip] * px[ip+1] - px[ip] * pz[ip+1]
         cz = px[ip] * py[ip+1] - py[ip] * px[ip+1]
         a += sqrt(cx * cx + cy * cy + cz * cz)
      end
      a *= 0.5
      
   else
      
      # open polygon: divide by length
      a = 0.0
      for ip = 1:np
         dx = polygon[ip+1,1] - polygon[ip,1]
         dy = polygon[ip+1,2] - polygon[ip,2]
         dz = polygon[ip+1,3] - polygon[ip,3]
         a += sqrt(dx * dx + dy * dy + dz * dz)
      end
      
   end
   
   return a
end  # getNormalizeVal
