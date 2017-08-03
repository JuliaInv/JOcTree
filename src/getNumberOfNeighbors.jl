export getNumberOfNeighbors

function  getNumberOfNeighbors(S::SparseArray3D)
# [upper,lower,left,right,front,back] = getNumberOfNeighbors(S)

m1,m2,m3  = S.sz
i,j,k,bsz = find3(S)

ns = length(i)
left  = zeros(Int64,ns)
right = zeros(Int64,ns)
upper = zeros(Int64,ns)
lower = zeros(Int64,ns)
front = zeros(Int64,ns)
back  = zeros(Int64,ns)

## UPPER ================================================
Iin   = find((i-bsz).>=1)
Itmp  = sub2ind(S.sz,i[Iin]-bsz[Iin],j[Iin],k[Iin])
upper[Iin] = 4
for kk = 1:length(Itmp)
   ik = Iin[kk]
   v = S.SV[Itmp[kk],1]
   if v == bsz[ik] || v == 0
      upper[ik] = 1
   end
end

##   CHECK FOR 4 NEIGHBORS IF i-sz/2 >=1
Iin   = find(((i-div.(bsz,2)) .>= 1) .& bsz.>1 )
Itmp  = sub2ind(S.sz,i[Iin]-div.(bsz[Iin],2),j[Iin],k[Iin]);
for kk = 1:length(Iin); if S.SV[Itmp[kk]] > 0; upper[Iin[kk]] = 4; end; end

## LOWER ==============================================
Iin  = find((i+bsz).<=m1)
Itmp = sub2ind(S.sz,i[Iin]+bsz[Iin],j[Iin],k[Iin])
lower[Iin] = 4
for kk = 1:length(Itmp)
   ik = Iin[kk]
   v = S.SV[Itmp[kk],1]
   bb = bsz[ik]
   if v == 2*bb || v == bb || v == 0
      lower[ik] = 1
   end
end

## LEFT ==================================
## IF j-sz >= 1
Iin   = find((j-bsz) .>= 1)
Itmp  = sub2ind(S.sz,i[Iin],j[Iin]-bsz[Iin],k[Iin])
left[Iin] = 4
##     CHECK FOR ONLY 1 NEIGHBOR
for kk = 1:length(Itmp)
   ik = Iin[kk]
   v = S.SV[Itmp[kk],1]
   if v == bsz[ik] || v == 0
      left[ik] = 1
   end
end
##   CHECK FOR 4 NEIGHBORS IF j-sz/2 >=1
Iin   = find(((j-div.(bsz,2)) .>= 1) .&  bsz.>1 )
Itmp  = sub2ind(S.sz,i[Iin],j[Iin]-div.(bsz[Iin],2),k[Iin])
for kk = 1:length(Itmp)
        if S.SV[Itmp[kk],1] > 0; left[Iin[kk]] = 4; end
end
left[Iin[find(S.SV[Itmp,1].>0)]] = 4

## RIGHT ===============================================================
Iin   = find(j+bsz.<=m2)
Itmp = sub2ind(S.sz,i[Iin],j[Iin]+bsz[Iin],k[Iin])
right[Iin] = 4
for kk = 1:length(Itmp)
   ik = Iin[kk]
   v = S.SV[Itmp[kk],1]
   bb = bsz[ik]
   if v == 2*bb || v == bb || v == 0
      right[ik] = 1
   end
end

## Front ==============================================================
Iin   = find((k-bsz) .>= 1)
Itmp  = sub2ind(S.sz,i[Iin],j[Iin],k[Iin]-bsz[Iin])
front[Iin] = 4
##     CHECK FOR ONLY 1 NEIGHBOR
for kk = 1:length(Itmp)
   ik = Iin[kk]
   v = S.SV[Itmp[kk],1]
   if v == bsz[ik] || v == 0
      front[ik] = 1
   end
end
##   CHECK FOR 4 NEIGHBORS IF j-sz/2 >=1
Iin   = find(((k-div.(bsz,2)) .>= 1)  .& bsz.>1 )
Itmp  = sub2ind(S.sz,i[Iin],j[Iin],k[Iin]-div.(bsz[Iin],2))
front[Iin[find(S.SV[Itmp,1].>0)]] = 4

## Back ==========================================================
Iin   = find(k+bsz.<=m3)
Itmp = sub2ind(S.sz,i[Iin],j[Iin],k[Iin]+bsz[Iin])
back[Iin] = 4
for kk = 1:length(Itmp)
   ik = Iin[kk]
   v = S.SV[Itmp[kk],1]
   bb = bsz[ik]
   if v == 2*bb || v == bb || v == 0
     back[ik] = 1
   end
end

return upper, lower, left, right, front, back

end
