export getFaceMassMatrix #, getdFaceMassMatrix (derivative not yet implemented 
                                               #for nodal interpolation approach to mass matrices



function getFaceMassMatrix(M::OcTreeMeshFV,sigma::Vector)
    if isempty(M.Pf)
        M.Pf = getFaceMassMatrixIntegrationMatrix(M.S,M.h)
    end
    P = M.Pf    
    n = length(sigma)
    if n == M.nc
        M = P' * kron(speye(24),spdiagm(sigma)) * P
    elseif n == 3 * M.nc
        M = P' * kron(speye(8),spdiagm(sigma)) * P
    elseif n == 6 * M.nc
        S11 = spdiagm(sigma[:,1])
        S22 = spdiagm(sigma[:,2])
        S33 = spdiagm(sigma[:,3])
        S12 = spdiagm(sigma[:,4])
        S13 = spdiagm(sigma[:,5])
        S23 = spdiagm(sigma[:,6])    
        Sig = [S11 S12 S13;
               S12 S22 S23;
               S13 S23 S33]
        M = P'* kron(speye(8),Sig) * P
    else
        error("Invalid size of sigma")
    end
    return M
end

function getFaceMassMatrixIntegrationMatrix(S::SparseArray3D,h)
    #  M = getFaceMassMatrix(S,h,sigma)
    #
    
    n = S.sz;
    nex = n + [1, 0, 0]
    ney = n + [0, 1, 0]
    nez = n + [0, 0, 1]
    
    i,j,k,bsz = find3(S)
    ex,ey,ez = getFaceNumbering(S)
    
    nx = nnz(ex)
    ny = nnz(ey)
    nz = nnz(ez)
    nc = nnz(S)
    c  = collect(1:nc)
    uc = sqrt(bsz.^3 * prod(h)/8)
    
    Px1 = ex.SV[sub2ind(nex,i,j,k),1]
    Py1 = ey.SV[sub2ind(ney,i,j,k),1]
    Pz1 = ez.SV[sub2ind(nez,i,j,k),1]
    
    Px2 = ex.SV[sub2ind(nex,i+bsz,j,k),1]
    Py2 = ey.SV[sub2ind(ney,i,j,k),1]
    Pz2 = ez.SV[sub2ind(nez,i,j,k),1]
    
    Px3 = ex.SV[sub2ind(nex,i,j,k),1]
    Py3 = ey.SV[sub2ind(ney,i,j+bsz,k),1]
    Pz3 = ez.SV[sub2ind(nez,i,j,k),1]
    
    Px4 = ex.SV[sub2ind(nex,i+bsz,j,k),1]
    Py4 = ey.SV[sub2ind(ney,i,j+bsz,k),1]
    Pz4 = ez.SV[sub2ind(nez,i,j,k),1]
    
    Px5 = ex.SV[sub2ind(nex,i,j,k),1]
    Py5 = ey.SV[sub2ind(ney,i,j,k),1]
    Pz5 = ez.SV[sub2ind(nez,i,j,k+bsz),1]
    
    Px6 = ex.SV[sub2ind(nex,i+bsz,j,k),1]
    Py6 = ey.SV[sub2ind(ney,i,j,k),1]
    Pz6 = ez.SV[sub2ind(nez,i,j,k+bsz),1]
    
    Px7 = ex.SV[sub2ind(nex,i,j,k),1]
    Py7 = ey.SV[sub2ind(ney,i,j+bsz,k),1]
    Pz7 = ez.SV[sub2ind(nez,i,j,k+bsz),1]
    
    Px8 = ex.SV[sub2ind(nex,i+bsz,j,k),1]
    Py8 = ey.SV[sub2ind(ney,i,j+bsz,k),1]
    Pz8 = ez.SV[sub2ind(nez,i,j,k+bsz),1]
    
    
    sp1(Q) = sparse(1:Base.nnz(Q),Q.nzval,uc,nc,nx)
    sp2(Q) = sparse(1:Base.nnz(Q),Q.nzval,uc,nc,ny)
    sp3(Q) = sparse(1:Base.nnz(Q),Q.nzval,uc,nc,nz)
    # sp1(Q) = sparse(1:Base.nnz(Q),Base.nonzeros(Q),ones(Base.nnz(Q)),nc,nx)
    # sp2(Q) = sparse(1:Base.nnz(Q),Base.nonzeros(Q),ones(Base.nnz(Q)),nc,ny)
    # sp3(Q) = sparse(1:Base.nnz(Q),Base.nonzeros(Q),ones(Base.nnz(Q)),nc,nz)
    
    P1 = blkdiag(sp1(Px1),sp2(Py1),sp3(Pz1))
    P2 = blkdiag(sp1(Px2),sp2(Py2),sp3(Pz2))
    P3 = blkdiag(sp1(Px3),sp2(Py3),sp3(Pz3))
    P4 = blkdiag(sp1(Px4),sp2(Py4),sp3(Pz4))
    P5 = blkdiag(sp1(Px5),sp2(Py5),sp3(Pz5))
    P6 = blkdiag(sp1(Px6),sp2(Py6),sp3(Pz6))
    P7 = blkdiag(sp1(Px7),sp2(Py7),sp3(Pz7))
    P8 = blkdiag(sp1(Px8),sp2(Py8),sp3(Pz8))
    
    #v     = prod(h) * (bsz.^3)
    # Sigma = [spdiagm(sigma[:,1].*v)  spdiagm(sigma[:,4].*v)  spdiagm(sigma[:,5].*v);
    #          spdiagm(sigma[:,4].*v)  spdiagm(sigma[:,2].*v)  spdiagm(sigma[:,6].*v);
    #          spdiagm(sigma[:,5].*v)  spdiagm(sigma[:,6].*v)  spdiagm(sigma[:,3].*v)]
    # 
    # 
    # Massf = 0.125*(P1'*Sigma*P1 + P2'*Sigma*P2 + P3'*Sigma*P3 + P4'*Sigma*P4 +
    #                P5'*Sigma*P5 + P6'*Sigma*P6 + P7'*Sigma*P7 + P8'*Sigma*P8);
    
    #P = [P1;P2;P3;P4;P5;P6;P7;P8]
    # Following line faster than line above because hcat + transpose is
    # faster than vcat.
    P = [P1' P2' P3' P4' P5' P6' P7' P8']'
    return P
    
end  # function getFaceMassMatrix

