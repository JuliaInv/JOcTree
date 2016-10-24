export OcTreeMeshFEM, getOcTreeMeshFEM

type OcTreeMeshFEM <: OcTreeMesh
	S::SparseArray3D
	h::Vector{Float64}  # (3) cell size
	x0::Vector{Float64}
	n::Vector{Int64}
	nc::Int          # number of cells
	nf::Vector{Int}  # (3) number of faces
	ne::Vector{Int}  # (3) number of edges
	Af::SparseMatrixCSC
	Ae::SparseMatrixCSC
	V::SparseMatrixCSC
	N::SparseMatrixCSC # nullspace matrix
	FX::SparseArray3D  # X face size
	FY::SparseArray3D  # Y face size
	FZ::SparseArray3D  # Z face size
	EX::SparseArray3D  # X edge size
	EY::SparseArray3D  # Y edge size
	EZ::SparseArray3D  # Z edge size
	NFX::SparseArray3D
	NFY::SparseArray3D
	NFZ::SparseArray3D
	NEX::SparseArray3D
	NEY::SparseArray3D
	NEZ::SparseArray3D	
end  # type OcTreeMeshFEM

	
function getOcTreeMeshFEM(S,h;x0=zeros(3))
	
		empt  = spzeros(0,0)
		 
		# get number of cells
		i,   = find3(S)
		nc   = size(i,1)
		
		# get number of faces
		FX,FY,FZ = getFaceSize(S)
		iex, = find3(FX)
		iey, = find3(FY)
		iez, = find3(FZ)
		nf   = [size(iex,1); size(iey,1); size(iez,1)]
		
		# get number of edges
		EX,EY,EZ = getEdgeSize(S)
		iex, = find3(EX)
		iey, = find3(EY)
		iez, = find3(EZ)
		ne   = [size(iex,1); size(iey,1); size(iez,1)]
		
		empt3 = sparse3([size(S,1),size(S,2),size(S,3)])
		
		return OcTreeMeshFEM(S,h,x0,
                     		S.sz,nc,nf,ne,
                     		empt,empt, empt,empt,  # no Af,Ae, V,N
   								FX,FY,FZ, EX,EY,EZ,
   								empt3,empt3,empt3,  # no NFX,NFY,NFZ
   							   empt3,empt3,empt3)  # no NEX,NEY,NEZ
		
end  # function getOcTreeMeshFEM

import Base.clear!
function clear!(M::OcTreeMeshFEM)
	M.S    = clear!(M.S)
	M.Af   = clear!(M.Af  )
	M.Ae   = clear!(M.Ae  )
	M.V    = clear!(M.V   )
	M.FX   = clear!(M.FX )
	M.FY   = clear!(M.FY )
	M.FZ   = clear!(M.FZ )
	M.N    = clear!(M.N )
	M.EX   = clear!(M.EX )
	M.EY   = clear!(M.EY )
	M.EZ   = clear!(M.EZ )
	M.NFX  = clear!(M.NFX)
	M.NFY  = clear!(M.NFY)
	M.NFZ  = clear!(M.NFZ)
	M.NEX  = clear!(M.NEX)
	M.NEY  = clear!(M.NEY)
	M.NEZ  = clear!(M.NEZ)
	
end  # function clear



import Base.==
function ==(M1::OcTreeMeshFEM, M2::OcTreeMeshFEM)
	isEqual = fill(true,21)


	# check mandatory fields
	isEqual[1] =  M1.S==M2.S
	isEqual[2] =  M1.h==M2.h
	isEqual[3] =  (M1.x0    == M2.x0)
	isEqual[5] =  (M1.n     == M2.n)
	isEqual[6] =  (M1.nc    == M2.nc)
	isEqual[7] =  (M1.nf    == M2.nf)
	isEqual[8] =  (M1.ne    == M2.ne)
	
	# check fields that might be empty
	if !(isempty(M1.Af)) && !(isempty(M2.Af))
		isEqual[12] = (M1.Af == M2.Af)
	end
	if !(isempty(M1.Ae)) && !(isempty(M2.Ae))
		isEqual[13] = (M1.Ae == M2.Ae)
	end
	if !(isempty(M1.V)) && !(isempty(M2.V))
		isEqual[15] = (M1.V == M2.V)
	end
	if !(isempty(M1.N)) && !(isempty(M2.N))
		isEqual[21] = (M1.N == M2.N)
	end
		
	return all(isEqual)
end  # function ==

include("getMatricesFEM.jl")


