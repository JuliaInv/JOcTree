export createOcTreeFromSrcRec

function createOcTreeFromSrcRec(x0,n,src,rec,h,w,sigback)
# 	createOcTreeFromSrcRec(src,rec,w,sigback)
	
	f = 2*pi*w 
	# compute skin depth
	sd = 500*sqrt(1/sigback/f)
	hRec = sd/5
	if minimum(h) > hRec
		warning("use smaller h, skin depth is",sd)
	end
	
	# define the box
	xmin = minimum([src[:,1];rec[:,1];rec[:,4]])
	xmax = maximum([src[:,1];rec[:,1];rec[:,4]])
	ymin = minimum([src[:,2];rec[:,2];rec[:,5]])
	ymax = maximum([src[:,2];rec[:,2];rec[:,5]])
	zmin = minimum([src[:,3];rec[:,3];rec[:,6]])
	zmax = maximum([src[:,3];rec[:,3];rec[:,6]])
						
	S = createOcTreeFromBox(
	  	x0[1], x0[2],x0[3],
		n[1],n[2],n[3],
		h[1],h[2],h[3],
	  	xmin, xmax,
	  	ymin, ymax,
	  	zmin, zmax,
	  	4,4)
		
	 S = regularizeOcTree(S)	
	 Msh  = getOcTreeMeshFV(S,h)	
	 	
	 EX,EY,EZ = getEdgeNumbering(MFor.S)
	 Smat     = getEdgeIntegralOfPolygonalChain(Msh,src,EX,EY,EZ)	
     Rmat     = getEdgeIntegralOfPolygonalChain(Msh,rec,EX,EY,EZ)		

	 return Msh, Smat, Rmat
end