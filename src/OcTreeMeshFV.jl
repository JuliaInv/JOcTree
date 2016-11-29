#module JOcTree


export OcTreeMeshFV, getOcTreeMeshFV

type OcTreeMeshFV <: OcTreeMesh
	S::SparseArray3D    # i,j,k, bsz
	h::Vector{Float64}  # (3) underlying cell size
	x0::Vector{Float64}
	n::Vector{Int64}
	nc::Int          # number of cells
	nf::Vector{Int}  # (3) number of faces
	ne::Vector{Int}  # (3) number of edges
	Div::SparseMatrixCSC
	Grad::SparseMatrixCSC
	Curl::SparseMatrixCSC
	Pf::SparseMatrixCSC # FaceMassMatrixIntegrationMatrix
	Pe::SparseMatrixCSC # EdgeMassMatrixIntegrationMatrix
	Af::SparseMatrixCSC # Face to cell-centre matrix
	Ae::SparseMatrixCSC # Edge to cell-centre matrix
	V::SparseMatrixCSC # cell volume
	L::SparseMatrixCSC # edge lengths
	Ne::SparseMatrixCSC # Edge nullspace matrix
	Qe::SparseMatrixCSC # Edge projection matrix
	Nn::SparseMatrixCSC # Node nullspace matrix
	Qn::SparseMatrixCSC # Node projection matrix
	Nf::SparseMatrixCSC # Face nullspace matrix
	Qf::SparseMatrixCSC # Face projection matrix
	FX::SparseArray3D  # X face size
	FY::SparseArray3D  # Y face size
	FZ::SparseArray3D  # Z face size
	EX::SparseArray3D  # X edge size
	EY::SparseArray3D  # Y edge size
	EZ::SparseArray3D  # Z edge size
	NFX::SparseArray3D # X FaceNumbering
	NFY::SparseArray3D # Y FaceNumbering
	NFZ::SparseArray3D # Z FaceNumbering
	NEX::SparseArray3D # X EdgeNumbering
	NEY::SparseArray3D # Y EdgeNumbering
	NEZ::SparseArray3D # Z EdgeNumbering
end # type OcTreeMeshFV

	
function getOcTreeMeshFV(S,h;x0=zeros(3))
	
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
		return OcTreeMeshFV(S, h, x0,
                    		  S.sz,nc,nf,ne,
                    		  empt,empt,empt,       # no Div, Grad, Curl
                    		  empt,empt,empt,empt,empt,empt,empt,empt,  # no Pf,Pe,Af,Ae,V,L,Ne,Qe
                    		  empt,empt,empt,empt, #no Nn,Qn,Nf,Qf
   							  FX,FY,FZ, EX,EY,EZ,
   							  empt3,empt3,empt3,  # no NFX,NFY,NFZ
   							  empt3,empt3,empt3)  # no NEX,NEY,NEZ	
end  # function getOcTreeMeshFV



import Base.clear!
function clear!(M::OcTreeMeshFV)
	M.S    = clear!(M.S)
	M.Div  = clear!(M.Div )
	M.Grad = clear!(M.Grad)
	M.Curl = clear!(M.Curl)
	M.Pf   = clear!(M.Pf  )
	M.Pe   = clear!(M.Pe  )
	M.Ae   = clear!(M.Ae  )
	M.Af   = clear!(M.Af  )
	M.V    = clear!(M.V   )
	M.L    = clear!(M.L   )
	M.FX   = clear!(M.FX )
	M.FY   = clear!(M.FY )
	M.FZ   = clear!(M.FZ )
	M.Ne   = clear!(M.Ne )
	M.Qe   = clear!(M.Qe )
	M.Nn   = clear!(M.Nn )
	M.Qn   = clear!(M.Qn )
	M.Nf   = clear!(M.Nf )
	M.Qf   = clear!(M.Qf )
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
function ==(M1::OcTreeMeshFV,M2::OcTreeMeshFV)
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
	if !(isempty(M1.Div)) && !(isempty(M2.Div))
		isEqual[9] = (M1.Div == M2.Div)
	end
	if !(isempty(M1.Grad)) && !(isempty(M2.Grad))
		isEqual[10] = (M1.Grad == M2.Grad)
	end
	if !(isempty(M1.Curl)) && !(isempty(M2.Curl))
		isEqual[11] = (M1.Curl == M2.Curl)
	end
	if !(isempty(M1.Af)) && !(isempty(M2.Af))
		isEqual[12] = (M1.Af == M2.Af)
	end
	if !(isempty(M1.Ae)) && !(isempty(M2.Ae))
		isEqual[13] = (M1.Ae == M2.Ae)
	end
	if !(isempty(M1.V)) && !(isempty(M2.V))
		isEqual[15] = (M1.V == M2.V)
	end
	if !(isempty(M1.Ne)) && !(isempty(M2.Ne))
		isEqual[21] = (M1.Ne == M2.Ne)
	end
	if !(isempty(M1.Nn)) && !(isempty(M2.Nn))
		isEqual[21] = (M1.Nn == M2.Nn)
	end
		
	return all(isEqual)
end  # function ==

