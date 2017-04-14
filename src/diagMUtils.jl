export DiagTimesM, MTimesDiag, DiagTimesMTimesDiag


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

#----------------------------------------------------

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

function DiagTimesMTimesDiag( d::Vector, A::SparseMatrixCSC, e::Vector )
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
