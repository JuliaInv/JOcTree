function Base.display(Mesh::OcTreeMesh)
	
	n    = Mesh.n
	nc   = Mesh.nc
	nf   = Mesh.nf
	ne   = Mesh.ne
	nn   = Mesh.nn
	x0   = Mesh.x0
	h    = Mesh.h
	d    = h .* n
	bsz  = nonzeros(Mesh.S)
	hmin = (isempty(bsz) ? 0 : minimum(bsz)) .* h
	hmax = (isempty(bsz) ? 0 : maximum(bsz)) .* h
	
	println("OcTree mesh of size $(n[1]) x $(n[2]) x $(n[3])")
	println("Number of cells:   $nc")
	println("Number of faces:   $(nf[1]) + $(nf[2]) + $(nf[3]) = $(sum(nf))")
	println("Number of edges:   $(ne[1]) + $(ne[2]) + $(ne[3]) = $(sum(ne))")
	println("Number of nodes:   $nn")
	println("Coordinate origin: ($(x0[1])m, $(x0[2])m, $(x0[3])m)")
	println("Domain size:       $(d[1])m x $(d[2])m x $(d[3])m")
	println("Minimum cell size: $(hmin[1])m x $(hmin[2])m x $(hmin[3])m")
	println("Maximum cell size: $(hmax[1])m x $(hmax[2])m x $(hmax[3])m")
	
end