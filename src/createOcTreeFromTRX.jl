export getOcTreeFromTRX

function getOcTreeFromTRX(trx::Array{Transmitter,1},cellSize,padding,depthFine;numFineLayers=2,numCoarseLayers=1,outFile=[])
	# trx = readdata(data_locations.txt; only_loc = true)
	#trx = readdata(dataFile)

	# find bounding box for survey area
	xmin =  Inf
	xmax = -Inf
	ymin =  Inf
	ymax = -Inf
	zmin =  Inf
	zmax = -Inf
	for trxi in trx
		xmin = min(xmin, minimum(trxi.trxpts[1,:]), minimum(trxi.rcvpts[1,:]))
		xmax = max(xmax, maximum(trxi.trxpts[1,:]), maximum(trxi.rcvpts[1,:]))
		ymin = min(ymin, minimum(trxi.trxpts[2,:]), minimum(trxi.rcvpts[2,:]))
		ymax = max(ymax, maximum(trxi.trxpts[2,:]), maximum(trxi.rcvpts[2,:]))
		zmin = min(zmin, minimum(trxi.trxpts[3,:]), minimum(trxi.rcvpts[3,:]))
		zmax = max(zmax, maximum(trxi.trxpts[3,:]), maximum(trxi.rcvpts[3,:]))
	end

	if zmin != zmax
		error("Flat survey required")
	end

	hx = cellSize[1]
	hy = cellSize[2]
	hz = cellSize[3]

	# figure out domain size
	nx = nextpow2(ceil(Integer,((xmax - xmin) + 2.0 * padding) / hx))
	ny = nextpow2(ceil(Integer,((ymax - ymin) + 2.0 * padding) / hy))
	nz = nextpow2(ceil(Integer,((zmax - zmin) + 2.0 * padding) / hz))
	x0 = 0.5 * (xmin + xmax - nx * hx)
	y0 = 0.5 * (ymin + ymax - ny * hy)
	z0 = 0.5 * (zmin + zmax - nz * hz)

	# add depth of fine cells
	zmin -= depthFine

	# create Octree mesh from box
	S = createOcTreeFromBox(
		x0, y0, z0, nx, ny, nz, hx, hy, hz,
		xmin, xmax, ymin, ymax, zmin, zmax,
		numFineLayers, numCoarseLayers)

	M = getOcTreeMeshFV(S, [hx, hy, hz]; x0 = [x0, y0, z0])
	
	if !isempty(outFile)
		exportOcTreeRoman(outFile, M)	
	end
	return M

end

