function createQuadTreeFromImage(A,tol,minsz)
# [S,A] = createOcTreeFromImage(A,tol);
# A - image
# tol - equal intensity color

m1,m2 = size(A)


maxbsz  = min((m1,m2))
S       = ones(Uint64,m1,m2)
Amin    = copy(A)
Amax    = copy(A)

cnt = 1

for bsz = 2.^[0:round(log2(maxbsz))-1]
    i,j = find(S.==bsz)
    
    I = find( (mod(i-1,2*bsz).==0) & (mod(j-1,2*bsz).==0)  )

    if ~isempty[I]
        I00 = sub2ind(size(S),i[I]     , j[I]    )
        I01 = sub2ind(size(S),i[I]+bsz , j[I]    )
        I10 = sub2ind(size(S),i[I]     , j[I]+bsz)
        I11 = sub2ind(size(S),i[I]+bsz , j[I]+bsz)

        Amin[I00] = min([ Amin[I00] , Amin[I01] , Amin[I10] , Amin[I11] ] ,[],2)
        Amax[I00] = max([ Amax[I00] , Amax[I01] , Amax[I10] , Amax[I11] ] ,[],2)
        
        Ic = find( (Amax[I00]-Amin[I00]) .<= tol );
        
        S(I00[Ic]) = 2*bsz;
        S(I01[Ic]) = 0;
        S(I10[Ic]) = 0;
        S(I11[Ic]) = 0;
    end;
    
end

S = sparse(S);
%% -------------------------------------------------------
%% -------------------------------------------------------
%% -------------------------------------------------------
function ret = isPowerOfTwo(n)
if 2^floor(log2(n)) == n
    ret = 1;
else 
    ret = 0;
end;