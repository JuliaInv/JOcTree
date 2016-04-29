export getNumberOfNeighbors

function  getNumberOfNeighbors(S)
# [upper,lower,left,right,front,back] = getNumberOfNeighbors(S)

m1,m2,m3  = S.sz
i,j,k,bsz = find3(S)

left  = round(Int64,zeros(length(i)))
right = round(Int64,zeros(length(i)))
upper = round(Int64,zeros(length(i)))
lower = round(Int64,zeros(length(i)))
front = round(Int64,zeros(length(i)))
back  = round(Int64,zeros(length(i)))

## UPPER ================================================
Iin   = find((i-bsz).>=1)
Itmp  = sub2ind(S.sz,i[Iin]-bsz[Iin],j[Iin],k[Iin])
upper[Iin] = 4
for kk = 1:length(Itmp)
        if S.SV[Itmp[kk],1] == bsz[Iin[kk]]; upper[Iin[kk]] = 1; end
        if S.SV[Itmp[kk],1] == 0;upper[Iin[kk]] = 1; end
end

##   CHECK FOR 4 NEIGHBORS IF i-sz/2 >=1
Iin   = find(((i-bsz/2) .>= 1) .== true & (mod(bsz,2) .== 0) .== true)
Itmp  = sub2ind(S.sz,i[Iin]-round(Int64,bsz[Iin]/2),j[Iin],k[Iin]);
for kk = 1:length(Iin); if S.SV[Itmp[kk]] > 0; upper[Iin[kk]] = 4; end; end

## LOWER ==============================================
Iin  = find((i+bsz).<=m1)
Itmp = sub2ind(S.sz,i[Iin]+bsz[Iin],j[Iin],k[Iin])
lower[Iin] = 4
for kk = 1:length(Itmp)
        if S.SV[Itmp[kk],1] == 2*bsz[Iin[kk]]; lower[Iin[kk]] = 1; end
        if S.SV[Itmp[kk],1] == bsz[Iin[kk]]; lower[Iin[kk]] = 1; end
        if S.SV[Itmp[kk],1] == 0;lower[Iin[kk]] = 1; end
end

## LEFT ==================================
## IF j-sz >= 1
Iin   = find((j-bsz) .>= 1)
Itmp  = sub2ind(S.sz,i[Iin],j[Iin]-bsz[Iin],k[Iin])
left[Iin] = 4
##     CHECK FOR ONLY 1 NEIGHBOR
for kk = 1:length(Itmp)
        if S.SV[Itmp[kk],1] == bsz[Iin[kk]]; left[Iin[kk]] = 1; end
        if S.SV[Itmp[kk],1] == 0;left[Iin[kk]] = 1; end
end
##   CHECK FOR 4 NEIGHBORS IF j-sz/2 >=1
Iin   = find(((j-bsz/2) .>= 1) .== true &  (mod(bsz,2) .== 0) .== true)
Itmp  = sub2ind(S.sz,i[Iin],j[Iin]-round(Int64,bsz[Iin]/2),k[Iin])
for kk = 1:length(Itmp)
        if S.SV[Itmp[kk],1] > 0; left[Iin[kk]] = 4; end
end
left[Iin[find(S.SV[Itmp,1].>0)]] = 4

## RIGHT ===============================================================
Iin   = find((j+bsz.<=m2))
Itmp = sub2ind(S.sz,i[Iin],j[Iin]+bsz[Iin],k[Iin])
right[Iin] = 4
for kk = 1:length(Itmp)
        if S.SV[Itmp[kk],1] == 2*bsz[Iin[kk]]; right[Iin[kk]] = 1; end
        if S.SV[Itmp[kk],1] == bsz[Iin[kk]]; right[Iin[kk]] = 1; end
        if S.SV[Itmp[kk],1] == 0;right[Iin[kk]] = 1; end
end

## Front ==============================================================
Iin   = find((k-bsz) .>= 1)
Itmp  = sub2ind(S.sz,i[Iin],j[Iin],k[Iin]-bsz[Iin])
front[Iin] = 4
##     CHECK FOR ONLY 1 NEIGHBOR
for kk = 1:length(Itmp)
        if S.SV[Itmp[kk],1] == bsz[Iin[kk]]; front[Iin[kk]] = 1; end
        if S.SV[Itmp[kk],1] == 0;front[Iin[kk]] = 1; end
end
##   CHECK FOR 4 NEIGHBORS IF j-sz/2 >=1
Iin   = find(((k-bsz/2) .>= 1) .== true & (mod(bsz,2) .== 0) .== true)
Itmp  = sub2ind(S.sz,i[Iin],j[Iin],k[Iin]-round(Int64,bsz[Iin]/2))
front[Iin[find(S.SV[Itmp,1].>0)]] = 4

## Back ==========================================================
Iin   = find((k+bsz.<=m3))
Itmp = sub2ind(S.sz,i[Iin],j[Iin],k[Iin]+bsz[Iin])
back[Iin] = 4
for kk = 1:length(Itmp)
        if S.SV[Itmp[kk],1] == 2*bsz[Iin[kk]]; back[Iin[kk]] = 1; end
        if S.SV[Itmp[kk],1] == bsz[Iin[kk]]; back[Iin[kk]] = 1; end
        if S.SV[Itmp[kk],1] == 0;back[Iin[kk]] = 1; end
end

return upper, lower, left, right, front, back

end
