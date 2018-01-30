export getFaceInterpolationMatrix

function getFaceInterpolationMatrix(mesh::OcTreeMesh, x::Array{Tf,1}, y::Array{Tf,1}, z::Array{Tf,1}) where Tf <: Real

	# map interpolation points to OcTree integer space
	x = (x - mesh.x0[1]) / mesh.h[1] + 1.0
	y = (y - mesh.x0[2]) / mesh.h[2] + 1.0
	z = (z - mesh.x0[3]) / mesh.h[3] + 1.0

	# interpolation point numbers
	Tn = typeof(mesh.nc)
	n = Tn(length(x))
	P = repmat([one(Tn):n;],2)

	# get face enumeration
	Fx,Fy,Fz = getFaceNumbering(mesh.S)
	nx = nnz(Fx)
	ny = nnz(Fy)
	nz = nnz(Fz)

	# locate points within cells
	i,j,k,bsz = findBlocks(mesh.S, floor.(Integer,x), floor.(Integer,y), floor.(Integer,z))

	# x-face numbers
	I  = [i, i+bsz;]
	J  = [j, j;]
	K  = [k, k;]
	Fx = vec(full(Fx.SV[sub2ind(size(Fx), I, J, K)]))

	# y-face numbers
	I  = [i, i;]
	J  = [j, j+bsz;]
	K  = [k, k;]
	Fy = vec(full(Fy.SV[sub2ind(size(Fy), I, J, K)]))

	# z-face numbers
	I  = [i, i;]
	J  = [j, j;]
	K  = [k, k+bsz;]
	Fz = vec(full(Fz.SV[sub2ind(size(Fz), I, J, K)]))

	# linear  basis
	u2 = (x - i) ./ bsz
	u1 = 1.0 - u2
	v2 = (y - j) ./ bsz
	v1 = 1.0 - v2
	w2 = (z - k) ./ bsz
	w1 = 1.0 - w2

	# linear interpolation
	Qx = [u1, u2;] # constant along y and z
	Qy = [v1, v2;] # constant along x and z
	Qz = [w1, w2;] # constant along x and y

	# assemble sparse interpolation matrices
	Qx = sparse(P, Fx, Qx, n, nx)
	Qy = sparse(P, Fy, Qy, n, ny)
	Qz = sparse(P, Fz, Qz, n, nz)

	return Qx,Qy,Qz

end
