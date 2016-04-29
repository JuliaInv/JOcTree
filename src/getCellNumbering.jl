export getCellNumbering

function getCellNumbering(S::SparseArray3D)

m1,m2,m3  = S.sz
i,j,k     = find3(S)

return sparse3(i,j,k,1:length(i),[m1;m2;m3;]);

end
