export SparseArray3D, sparse3, find3, nonzeros, nnz, getindex, setindex!

# Extend polymorphic methods in module Base
import Base.nonzeros
import Base.nnz
import Base.getindex
import Base.setindex!
import Base.ndims 

type SparseArray3D
        SV::SparseVector{Int64}
        sz::Vector{Int64}    # size of fine mesh
end

Base.size(S::SparseArray3D) = (S.sz[1], S.sz[2], S.sz[3])
Base.size(S::SparseArray3D,dim::Int) = S.sz[dim]
Base.find(S::SparseArray3D) = find(S.SV)
Base.ndims(S::SparseArray3D) = 3

import Base.clear!
clear!(S::SparseArray3D) = sparse3(S.sz)

function sparse3(sz::Vector{Int})
        S = spzeros(Int,prod(sz),1)
        return SparseArray3D(S,sz)
end

import Base.isempty
isempty(S::SparseArray3D) = (nnz(S.SV)==0)

function sparse3(i::Vector{Int},j::Vector{Int},k::Vector{Int},v::Range,sz::Vector{Int})
        return sparse3(i,j,k,[v;],sz)
end
function sparse3(i::Vector{Int},j::Vector{Int},k::Vector{Int},v::Range,sz::Vector{Int},combine::Function)
        return sparse3(i,j,k,[v;],sz,combine)
end

function sparse3(i::Vector{Int},j::Vector{Int},k::Vector{Int},v::Vector{Int},sz::Vector{Int})
	IND = sub2ind(sz,i,j,k)
   
  # Note that the following line would sort IND, and v would be permuted.  For
  # duplicate IND values, a SMALLEST v would be used.
  IND, v = sortpermFast(IND, v)
  S = sparsevec(IND,v, prod(sz))
	S3 = SparseArray3D(S,sz)
	return S3
end


function sparse3(i::Vector{Int},j::Vector{Int},k::Vector{Int},v::Vector{Int},sz::Vector{Int},combine::Function)

        IND = sub2ind(sz,i,j,k)
     	  S = sparsevec(IND,v, prod(sz), combine)
        return SparseArray3D(S,sz)
end

function find3(S::SparseArray3D)
        IND = find(S.SV)
        i,j,k = ind2sub((S.sz[1],S.sz[2],S.sz[3]),IND)
        return i, j, k, nonzeros(S.SV)
end

function nonzeros(S::SparseArray3D)
        return nonzeros(S.SV)
end

function nnz(S::SparseArray3D)
    return nnz(S.SV)
end

function getindex(S::SparseArray3D,i::Int,j::Int,k::Int)
	return S.SV[sub2ind(S.sz,i,j,k)]
end

function getindex(S::SparseArray3D,i::Vector{Int},j::Vector{Int},k::Vector{Int})
	sz = (length(i), length(j), length(k))
	I,J,K = ndgrid(i,j,k)
	si = sub2ind(S.sz,vec(I),vec(J),vec(K))
	idx = reshape(full(S.SV[si]), sz)   # SLOW  !!!!
	return idx
end

function setindex!(S::SparseArray3D, v::Int, i::Int, j::Int, k::Int)
	S.SV[sub2ind(S.sz,i,j,k)] = v
end

function setindex!(S::SparseArray3D, V::Vector{Int}, I::Vector{Int}, J::Vector{Int}, K::Vector{Int})
	for ind=1:length(I)
		S.SV[sub2ind(S.sz,I[ind],J[ind],K[ind])] = V[ind]
	end
end

import Base.==
function ==(S1::SparseArray3D,S2::SparseArray3D)
	return (S1.SV==S2.SV) && (S1.sz == S2.sz)
end
