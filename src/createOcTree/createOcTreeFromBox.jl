export createOcTreeFromBox

function createOcTreeFromBox(
  x0::AbstractFloat, y0::AbstractFloat, z0::AbstractFloat, 
  nx::Int, ny::Int, nz::Int,
  hx::AbstractFloat, hy::AbstractFloat, hz::AbstractFloat,
  xa::AbstractFloat, xb::AbstractFloat,
  ya::AbstractFloat, yb::AbstractFloat,
  za::AbstractFloat, zb::AbstractFloat,
  nf::Array{Int,1}, nc::Array{Int,1}, 
  bsz::Int = 1)
# S = createOcTreeFromBox( ...
#   x0, y0, z0, nx, ny, nz, hx, hy, hz, xa, xb, ya, yb, za, zb, nf, nc)
#
# Create OcTree with fine cells in box [xa,xb] x [ya,yb] x [za,zb]
# and pad by increasingly larger cells to fill the domain.
#
# Input:
# x0, y0, z0   coordinate origin of OcTree in meters
# nx, ny, nz   number of base mesh cells in x/y/z direction (powers of 2)
# hx, hy, hz   extension of a base mesh cell in meters
# xa, xb       x-limits of box that will be filled with fine cells
# ya, yb       y-limits of box that will be filled with fine cells
# za, zb       z-limits of box that will be filled with fine cells
# nf           minimum thickness of layer of fine cells around box
# nc           minimum thickness of padding layers of coarser cells
#
# Output:
# S            OcTree mesh (sparse3)
#
# Note: The padding by nf and nc cells may require more cells than
#       specified to fill up coarser cells.

# Christoph Schwarzbach, January 2014

# check input
if nextpow2(nx) != nx || nextpow2(ny) != ny || nextpow2(nz) != nz
  error("nx, ny, nz must be powers of 2.")
end

if !((length(nf) == 6) && all(nf .>= 0))
  error("nf must be an all non-negative vector of length 6")
end
if !((length(nc) == 6) && all(nc .> 0))
  error("nc must be an all positive vector of length 6")
end

if bsz < 1 || nextpow2(bsz) != bsz
	error("bsz must be greater or equal 1 and a power of 2")
end

# maximum number of cell sizes
nbsz = minimum(round(Int64,log2([nx, ny, nz]))) - round(Int64,log2(bsz))
if nbsz < 1
  error("bsz is too large")
end

# initialize (k,j,i,bsz) array
S = Array(Int, (0,4))

# convert box to integer grid
xa = (xa - x0) / hx + 1.0
xb = (xb - x0) / hx + 1.0
ya = (ya - y0) / hy + 1.0
yb = (yb - y0) / hy + 1.0
za = (za - z0) / hz + 1.0
zb = (zb - z0) / hz + 1.0

# start with fine cells around box, allow for small tolerance
i1 =     floor(Integer,xa) 
i2 = min(floor(Integer,xb), nx)
j1 =     floor(Integer,ya)
j2 = min(floor(Integer,yb), ny)
k1 =     floor(Integer,za)
k2 = min(floor(Integer,zb), nz)
t2 = sqrt(eps())
t1 = 1.0 - t2
if xa - i1 > t1 * i1
  i1 = i1 + 1
end
if xb - i2 < t2 * i2
  i2 = i2 - 1
end
if ya - j1 > t1 * j1
  j1 = j1 + 1
end
if yb - j2 < t2 * j2
  j2 = j2 - 1
end
if za - k1 > t1 * k1
  k1 = k1 + 1
end
if zb - k2 < t2 * k2
  k2 = k2 - 1
end
if i1 > i2
	i1,i2 = i2,i1
end
if j1 > j2
	j1,j2 = j2,j1
end
if k1 > k2
	k1,k2 = k2,k1
end

# make sure that indices live on grid with block size bsz
i1 -= mod(i1 - 1, bsz)
i2 -= mod(i2 - 1, bsz)
j1 -= mod(j1 - 1, bsz)
j2 -= mod(j2 - 1, bsz)
k1 -= mod(k1 - 1, bsz)
k2 -= mod(k2 - 1, bsz)

# initial number of cells per layer
nl = nf

# for all cell sizes, add padding cells until we fill the domain
for m = 1:nbsz
  
  # block size
  b = bsz * 2^(m-1)
  
  # pad
  if m < nbsz
    
    # pad by given number of cell layers
    # (i3 < i1 < i2 < i4, j3 < j1 < j2 < j4, k3 < k1 < k2 < k4)
    i3 = i1 - b * nl[1]
    i4 = i2 + b * nl[2]
    j3 = j1 - b * nl[3]
    j4 = j2 + b * nl[4]
    k3 = k1 - b * nl[5]
    k4 = k2 + b * nl[6]
    
    # make sure that box fills cells of size 2*h
    i3 = i3 - mod(i3 - 1,     b * 2)
    i4 = i4 + mod(i4 - 1 - b, b * 2)
    j3 = j3 - mod(j3 - 1,     b * 2)
    j4 = j4 + mod(j4 - 1 - b, b * 2)
    k3 = k3 - mod(k3 - 1,     b * 2)
    k4 = k4 + mod(k4 - 1 - b, b * 2)

    # domain limits
    i3 = max(1,          i3)
    i4 = min(nx + 1 - b, i4)
    j3 = max(1,          j3)
    j4 = min(ny + 1 - b, j4)
    k3 = max(1,          k3)
    k4 = min(nz + 1 - b, k4)
    
  else
    
    # fill whatever space is left
    i3 = 1
    i4 = nx + 1 - b
    j3 = 1
    j4 = ny + 1 - b
    k3 = 1
    k4 = nz + 1 - b
    
  end
  
  if m > 1
    
    # mask for inner region with finer cells, coordinate form
    u = [trues(round(Int64, (i1 - i3) / b)); falses(round(Int64, (i2 - i1) / b + 1)); trues(round(Int64, (i4 - i2) / b));]
    v = [trues(round(Int64, (j1 - j3) / b)); falses(round(Int64, (j2 - j1) / b + 1)); trues(round(Int64, (j4 - j2) / b));]
    w = [trues(round(Int64, (k1 - k3) / b)); falses(round(Int64, (k2 - k1) / b + 1)); trues(round(Int64, (k4 - k2) / b));]

    # we are done if mask is all false
    if !(any(u) || any(v) || any(w))
      break
    end
    
    # expand coordinate form to 3-D arrays
    u,v,w = ndgrid(u,v,w)
    
    # combine x/y/z masks
    q = u | v | w
    
  else
    
    # at the first, finest level, use all cells
    q = trues(round(Int64, (i4 - i3) / b + 1), round(Int64, (j4 - j3) / b + 1), round(Int64, (k4 - k3) / b + 1))
    
    # change to default number of cells
    nl = nc
    
  end
  
  # index vectors for current cell size
  i,j,k = ndgrid([i3:b:i4;], [j3:b:j4;], [k3:b:k4;])
  s = fill(b, size(i))
  
  # add only those cells that pad the inner, finer core region using mask
  S = [S; i[q] j[q] k[q] s[q];]
  
  # core region for next cell size, shift upper limit by current cell size
  # to position of next larger cell size
  i1 = i3
  i2 = i4 - b
  j1 = j3
  j2 = j4 - b
  k1 = k3
  k2 = k4 - b

end

# sanity check
if sum(S[:,4].^3) != nx * ny * nz
  error("Internal error: cell volumes don't fill the domain")
end

# convert to sparse3
 # S = sortrows(S)
 # S = sparse3(S[:,1], S[:,2], S[:,3], S[:,4], [nx, ny, nz])

S = sparse3(S[:,1], S[:,2], S[:,3], S[:,4], [nx; ny; nz;])

return S

end

createOcTreeFromBox(
  x0::AbstractFloat, y0::AbstractFloat, z0::AbstractFloat, 
  nx::Int, ny::Int, nz::Int,
  hx::AbstractFloat, hy::AbstractFloat, hz::AbstractFloat,
  xa::AbstractFloat, xb::AbstractFloat,
  ya::AbstractFloat, yb::AbstractFloat,
  za::AbstractFloat, zb::AbstractFloat,
  nf::Int,
  nc::Int,
  bsz::Int = 1) =
createOcTreeFromBox(
  x0, y0, z0,
  nx, ny, nz, 
  hx, hy, hz,
  xa, xb, ya, yb, za, zb,
  [nf; nf; nf; nf; nf; nf;],
  [nc; nc; nc; nc; nc; nc;],
  bsz)