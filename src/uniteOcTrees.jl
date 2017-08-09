export uniteOcTrees

function uniteOcTrees(S1::SparseArray3D, S2::SparseArray3D)
	
	if size(S1) != size(S2)
		error("Incompatible size")
	end
	SZ = S1.sz
	
	T1 = sparsevec(S1.SV.nzind, fill(1,length(S1.SV.nzval)), S1.SV.n)
	T2 = sparsevec(S2.SV.nzind, fill(2,length(S2.SV.nzval)), S2.SV.n)
	T  = T1+T2
	
	M = find(T)
	V = T.nzval
	
	j1 = 1
	j2 = 1
	for k=1:length(M)
		if V[k]==3
			V[k]= min(S1.SV.nzval[j1],S2.SV.nzval[j2])
			j2+=1
			j1+=1
		elseif V[k]==2
			V[k] = S2.SV.nzval[j2]
			j2+=1
		elseif V[k]==1
			V[k] = S1.SV.nzval[j1]
			j1+=1
		end
	end	# 
	I,J,K = ind2sub(size(S1), M)
	S = sparse3(I,J,K,V,SZ)
	return S
end

function uniteOcTrees(M1::OcTreeMeshFV, M2::OcTreeMeshFV)
    MU_S = uniteOcTrees(M1.S, M2.S)
    MU = getOcTreeMeshFV(MU_S, M1.h; x0=M1.x0);
    return MU
end

function uniteOcTrees(SL::Array{SparseArray3D,1})
  n = length(SL)
  if n > 2
    n1 = div(n,2)
    S1 = uniteOcTrees(SL[1:n1])
    S2 = uniteOcTrees(SL[n1+1:end])
    return uniteOcTrees(S1,S2)
  elseif n == 2
    return uniteOcTrees(SL[1],SL[2])
  else
    return SL[1]
  end
end