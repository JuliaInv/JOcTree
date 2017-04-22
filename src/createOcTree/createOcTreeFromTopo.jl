using Mesh

function createOcTreeFromTopo(Z,n,h,x0,boxXY)
# Z a matrix of size n[1]*n[2]

A = zeros(UInt8,n[1],n[2],n[3])

i1 = int(round((boxXY[1,1]-x0[1])/h[1])) + 1
i2 = int(round((boxXY[1,2]-x0[1])/h[1])) + 1
j1 = int(round((boxXY[2,1]-x0[2])/h[2])) + 1
j2 = int(round((boxXY[2,2]-x0[2])/h[2])) + 1

	
for i=1:n[1]
	for j=1:n[2]	 		
		k = int(round((Z[i,j]-x0[3])/h[3])) + 1
		if k<1
			error("Topo under the mesh")
		end
		if k>n[3]
			error("Topo above the mesh")
		end
		if (i>i1) && (i<i2) && (j>j1) && (j<j2)
			A[i,j,1:k] = 1
		end	
		
	end
end
	
	S = regularizeOcTree(createOcTreeFromImage(A,.1))
	M = getOcTreeMeshFV(S,vec(h);x0=vec(x0))
	X = getCellCenteredGrid(M)
   Xc = X[:,1]; 
	Yc = X[:,2]; 
	Zc = X[:,3]	
	
	# move to indices
	Xc = (Xc - x0[1])/h[1]+1
   Yc = (Yc - x0[2])/h[2]+1 
		
	iact = zeros(UInt8,size(X,1))
	for k=1:length(iact)
		zk = Zc[k]
		i  = int(Xc[k]); j = int(Yc[k])
		height = Z[i,j]
		if height > zk
			iact[k] = 1
		end
	end
	
	return S, iact
end
		