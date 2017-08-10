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

function getNodesFromIndices(sv,mm,i0::Vector{Int},j0::Vector{Int},k0::Vector{Int})

	jj = sub2ind(mm,i0,j0,k0)
  v  = Array{Int64}(length(jj))
  for i = 1:length(jj)
    v[i] = sv[jj[i]]
  end
	return v

end
