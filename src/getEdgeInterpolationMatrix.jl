export getEdgeInterpolationMatrix

function getEdgeInterpolationMatrix(mesh::OcTreeMesh, x::Array{Float64,1}, y::Array{Float64,1}, z::Array{Float64,1})
	
	# map interpolation points to OcTree integer space
	x = (x - mesh.x0[1]) / mesh.h[1] + 1.0
	y = (y - mesh.x0[2]) / mesh.h[2] + 1.0
	z = (z - mesh.x0[3]) / mesh.h[3] + 1.0
	
	# interpolation point numbers
	n = length(x)
	P = repmat([1:n;],4)
	
	# get edge enumeration
	Ex,Ey,Ez = getEdgeNumbering(mesh.S)
	nx = nnz(Ex)
	ny = nnz(Ey)
	nz = nnz(Ez)
	
	# locate points within cells
	i,j,k,bsz = findBlocks(mesh.S, floor(Integer,x), floor(Integer,y), floor(Integer,z))
	
	# x-edge numbers 
	I  = [i, i, i, i;]
	J  = [j, j+bsz, j, j+bsz;]
	K  = [k, k, k+bsz, k+bsz;]
	Ex = vec(full(Ex.SV[sub2ind(size(Ex), I, J, K)]))
	
	# y-edge numbers 
	I  = [i, i+bsz, i, i+bsz;]
	J  = [j, j, j, j;]
	K  = [k, k, k+bsz, k+bsz;]
	Ey = vec(full(Ey.SV[sub2ind(size(Ey), I, J, K)]))
	
	# z-edge numbers 
	I  = [i, i+bsz, i, i+bsz;]
	J  = [j, j, j+bsz, j+bsz;]
	K  = [k, k, k, k;]
	Ez = vec(full(Ez.SV[sub2ind(size(Ez), I, J, K)]))
	
	# linear basis
	u2 = (x - i) ./ bsz
	u1 = 1.0 - u2
	v2 = (y - j) ./ bsz
	v1 = 1.0 - v2
	w2 = (z - k) ./ bsz
	w1 = 1.0 - w2
	
	# bilinear interpolation
	Qx = [v1 .* w1, v2 .* w1, v1 .* w2, v2 .* w2;] # constant along x
	Qy = [u1 .* w1, u2 .* w1, u1 .* w2, u2 .* w2;] # constant along y
	Qz = [u1 .* v1, u2 .* v1, u1 .* v2, u2 .* v2;] # constant along z
	
	# assemble sparse interpolation matrices
	Qx = sparse(P, Ex, Qx, n, nx)
	Qy = sparse(P, Ey, Qy, n, ny)
	Qz = sparse(P, Ez, Qz, n, nz)

	return Qx,Qy,Qz
	
end