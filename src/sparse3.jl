export SparseArray3D, sparse3, find3, nonzeros, nnz, getindex, setindex!

# Extend polymorphic methods in module Base
import Base.nonzeros
import Base.nnz
import Base.getindex
import Base.setindex!
import Base.ndims
import Base.convert

mutable struct SparseArray3D{N<:Integer,N2<:Integer}
        SV::SparseVector{N,N2}
        sz::Tuple{N2,N2,N2}    # size of fine mesh
end

function SparseArray3D{N<:Integer,N2<:Integer}(SV::SparseVector{N,N2},sz::Tuple)
    sz = (convert(N2,sz[1]), convert(N2,sz[2]), convert(N2,sz[3]))
    SparseArray3D(SV,sz)
end

Base.size(S::SparseArray3D) = S.sz
Base.size(S::SparseArray3D,dim::Integer) = S.sz[dim]
Base.find(S::SparseArray3D) = find(S.SV)
Base.ndims(S::SparseArray3D) = 3

function convert(::Type{N}, ::Type{N2}, S::SparseArray3D) where {N <: Integer, N2 <: Integer}
    sz = (N2(S.sz[1]),N2(S.sz[2]),N2(S.sz[3]))
    return SparseArray3D(convert(SparseVector{N,N2}, S.SV), sz)
end

function convert(::Type{N}, S::SparseArray3D) where N <: Integer
    return convert(N, eltype(S.SV.nzind), S)
end

function sparse3(::Type{N}, sz::Tuple{N2,N2,N2}) where {N <: Integer, N2 <: Integer}
    S = spzeros(N,N2,prod(sz))
    return SparseArray3D(S,sz)
end

function sparse3(sz::Tuple{N2,N2,N2}) where N2 <: Integer
  S = spzeros(N2,N2,prod(sz))
  return SparseArray3D(S,sz)
end

import Base.isempty
isempty(S::SparseArray3D) = (nnz(S.SV)==0)

function sparse3(i::Vector,j::Vector,k::Vector,v::Range,sz::Tuple)
        return sparse3(i,j,k,[v;],sz)
end
function sparse3(i::Vector,j::Vector,k::Vector,v::Range,sz::Tuple,combine::Function)
        return sparse3(i,j,k,[v;],sz,combine)
end

function sparse3(i::Vector,j::Vector,k::Vector,v::Vector,sz::Tuple)
	IND = sub2ind(sz,i,j,k)

  # Note that the following line would sort IND, and v would be permuted.  For
  # duplicate IND values, a SMALLEST v would be used.
  IND, v = sortpermFast(IND, v)
  S = sparsevec(IND,v, prod(sz))
	S3 = SparseArray3D(S,sz)
	return S3
end


function sparse3(i::Vector,j::Vector,k::Vector,v::Vector,sz::Tuple,combine::Function)

        IND = sub2ind(sz,i,j,k)
     	  S = sparsevec(IND,v, prod(sz), combine)
        return SparseArray3D(S,sz)
end

function find3(S::SparseArray3D)
        IND = find(S.SV)
        i,j,k = ind2sub(S.sz,IND)
        # Tn = eltype(S.SV.nzval)
        # i = convert(Vector{Tn},i)
        # j = convert(Vector{Tn},j)
        # k = convert(Vector{Tn},k)
        return i, j, k, nonzeros(S.SV)
end

function nonzeros(S::SparseArray3D)
        return nonzeros(S.SV)
end

function nnz(S::SparseArray3D)
    return nnz(S.SV)
end

function getindex(S::SparseArray3D,i::Integer,j::Integer,k::Integer)
	return S.SV[sub2ind(S.sz,i,j,k)]
end

function getindex(S::SparseArray3D,i::Vector{Integer},j::Vector{Integer},
                  k::Vector{Integer})
	sz = (length(i), length(j), length(k))
	I,J,K = ndgrid(i,j,k)
	si = sub2ind(S.sz,vec(I),vec(J),vec(K))
	idx = reshape(full(S.SV[si]), sz)   # SLOW  !!!!
	return idx
end

function setindex!(S::SparseArray3D, v::Integer, i::Integer, j::Integer,
                   k::Integer)
	S.SV[sub2ind(S.sz,i,j,k)] = v
end

function setindex!(S::SparseArray3D, V::Vector{Integer}, I::Vector{Integer},
                   J::Vector{Integer}, K::Vector{Integer})
	for ind=1:length(I)
		S.SV[sub2ind(S.sz,I[ind],J[ind],K[ind])] = V[ind]
	end
end

import Base.==
function ==(S1::SparseArray3D,S2::SparseArray3D)
	return (S1.SV==S2.SV) && (S1.sz == S2.sz)
end

import Base.clear!
function clear!{N<:Integer,N2<:Integer}(S::SparseArray3D{N,N2})
  S.SV = spzeros(N,N2,prod(S.sz))
  return
end
