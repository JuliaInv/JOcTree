export DiagTimesM, MTimesDiag, DiagTimesMTimesDiag!

# Exported utilities

"""
A = DiagTimesM(d,A)

Computes D`*`A, the product of a diagonal matrix D and a general
sparse matrix A in an optimized way, with A being updated in place.

Input:

        d::Vector -- Diagonal matrix D, stored as a vector
                     of the diagonal elements. D = spdiagm(d).
        A::SparseMatrixCSC

Output:

        A::SparseMatrixCSC -- D*A, with input A overwritten
"""
function DiagTimesM( d::Vector, A::SparseMatrixCSC )
# Return diag{d} * A

   n = size(A,2)
   if (length(d) != size(A,1))
      error("length(d) != size(A,1)")
   end

   for ir = 1:n
      j1 = A.colptr[ir]
      j2 = A.colptr[ir+1] - 1
      for ic = j1:j2
         jcol = A.rowval[ic]
         A.nzval[ic] *= d[jcol]
      end
   end

   return A
end



"""
A = MTimesDiag(A,d)

Computes A`*`D, the product of a general
sparse matrix A and a diagonal matrix D in an optimized way,
with A being updated in place.

Input:

        A::SparseMatrixCSC
        d::Vector -- Diagonal matrix D, stored as a vector
                     of the diagonal elements. D = spdiagm(d).

Output:

        A::SparseMatrixCSC -- A*D, with input A overwritten
"""
function MTimesDiag( A::SparseMatrixCSC, d::Vector )
# Return A * diag{d}

   n = size(A,2)
   if (length(d) != n)
      error("length(d) != n")
   end

   for ir = 1:n
      j1 = A.colptr[ir]
      j2 = A.colptr[ir+1] - 1
      dgir = d[ir]
      for ic = j1:j2
         A.nzval[ic] *= dgir
      end
   end

   return A
end


"""
A = DiagTimesMTimesDiag!(d,A,e)

Equivalent output to spdiagm(d)`*`A`*`spdiagm(e), computed in an
optimized way, with A updated in-place.

Input:

       d::Vector -- Diagonal elements of a Diagonal matrix, stored
                    as a vector
       A::SparseMatrixCSC
       e::Vector -- Diagonal elements of a Diagonal matrix, stored
                    as a vector

Output:

       A::SparseMatrixCSC -- Result of spdiagm(d)\*A\*spdiagm(e),
                             overwrites input matrix A.

"""
function DiagTimesMTimesDiag!( d::Vector, A::SparseMatrixCSC, e::Vector )
# Return diag{d} * A * diag{e}

   n = size(A,2)
   if (length(d) != size(A,1) || length(e) != n)
    error("length(d) != size(A,1) || length(e) != n")
   end

   for ir = 1:n
      j1 = A.colptr[ir]
      j2 = A.colptr[ir+1] - 1
      dge = e[ir]
      for ic = j1:j2
         jcol = A.rowval[ic]
         A.nzval[ic] = d[jcol] * A.nzval[ic] * dge
      end
   end

   return A
end

#-----------------------------------------------------------------

# Internal utilities

function getNodesFromIndices(sv,mm::Tuple,i0::Vector,j0::Vector,k0::Vector)
    Ti = promote_type(eltype(mm),promote_type(eltype(i0),
                      promote_type(eltype(j0),eltype(k0))))
    mm = (Ti(mm[1]),Ti(mm[2]),Ti(mm[3]))
    i0 = convert(Vector{Ti},i0)
    j0 = convert(Vector{Ti},j0)
    k0 = convert(Vector{Ti},k0)
	jj = sub2ind(mm,i0,j0,k0)
    v  = Array{eltype(sv)}(length(jj))
    for i = 1:length(jj)
        v[i] = sv[jj[i]]
    end
	return v
end

"""
    merge!(a::Vector{Int64}, b::Vector{Int64})

Replace zero entries in `a` by values from `b`.
"""
function merge!(a::Vector{T}, b::Vector{T}) where T <: Number
  n = length(a)
  (length(b) == n) || throw(DimensionMismatch("length(a) != length(b)"))
  @inbounds begin
    for i = 1:n
      if a[i] == 0
        a[i] = b[i]
      end
    end
  end
end

import Base.speye
speye(Tt::Type{T}, Tn::Type{N}, m::Integer) where T where N = speye(Tt, Tn, m, m)
function speye(::Type{T}, ::Type{N}, m::Integer, n::Integer) where T where N
    ((m < 0) || (n < 0)) && throw(ArgumentError("invalid array dimensions"))
    nnz = min(m,n)
    colptr = Vector{N}(1+n)
    colptr[1:nnz+1] = 1:nnz+1
    colptr[nnz+2:end] = nnz+1
    SparseMatrixCSC(Int(m), Int(n), colptr, Vector{N}(1:nnz), ones(T,nnz))
end
