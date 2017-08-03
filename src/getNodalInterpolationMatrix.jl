export getNodalInterpolationMatrix

function getNodalInterpolationMatrix(mesh::OcTreeMesh, x::Array{Float64,1}, y::Array{Float64,1}, z::Array{Float64,1})
	
	# map interpolation points to OcTree integer space
	x = (x - mesh.x0[1]) / mesh.h[1] + 1.0
	y = (y - mesh.x0[2]) / mesh.h[2] + 1.0
	z = (z - mesh.x0[3]) / mesh.h[3] + 1.0
	
	# interpolation point numbers
	n = length(x)
	P = repmat([1:n;],8)
	
	# get nodal enumeration
	N = getNodalNumbering(mesh)
	nn = nnz(N)
	
	# locate points within cells
	i,j,k,bsz = findBlocks(mesh.S, floor.(Integer,x), floor.(Integer,y), floor.(Integer,z))
	
	# node numbers 
	I = [i, i+bsz, i, i+bsz, i, i+bsz, i, i+bsz;]
	J = [j, j, j+bsz, j+bsz, j, j, j+bsz, j+bsz;]
	K = [k, k, k, k, k+bsz, k+bsz, k+bsz, k+bsz;]
	N = vec(full(N.SV[sub2ind(size(N), I,J,K)]));
	
	# trilinear interpolation
	xd = (x-i)./bsz
	yd = (y-j)./bsz
	zd = (z-k)./bsz

	u000 = (1.-xd).*(1.-yd).*(1.-zd)
	u100 = xd.*(1.-yd).*(1.-zd)
	u010 = (1.-xd).*yd.*(1.-zd)
	u110 = xd.*yd.*(1.-zd)
	u001 = (1.-xd).*(1.-yd).*zd
	u101 = xd.*(1.-yd).*zd
	u011 = (1.-xd).*yd.*zd
	u111 = xd.*yd.*zd

	Q = [u000, u100, u010, u110, u001, u101, u011, u111;]

	# assemble sparse interpolation matrix
	Q = sparse(P, N, Q, n, nn)

	return Q
	
end