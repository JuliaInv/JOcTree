
export getSizeOfNeighbors

function getSizeOfNeighbors(S::SparseArray3D)
    m1, m2, m3 = S.sz
    i, j, k, bsz = find3(S)

    ns = length(i)
    left = zeros(Float64, ns)
    right = zeros(Float64, ns)
    upper = zeros(Float64, ns)
    lower = zeros(Float64, ns)
    front = zeros(Float64, ns)
    back = zeros(Float64, ns)

    ## UPPER
    Iin = find((i-bsz).>=1)
    Itmp = sub2ind(S.sz, i[Iin]-bsz[Iin] ,j[Iin],k[Iin])
    for kk = 1:length(Itmp)
        ik = Iin[kk]
        v = S.SV[Itmp[kk]]
        if v == bsz[ik]
            upper[ik] = 1.
        elseif v == 0
            upper[ik] = 2.
        else
            upper[ik] = 0.5
        end
    end

    ## LOWER
    Iin = find((i+bsz) .<= m1)
    Itmp = sub2ind(S.sz, i[Iin]+bsz[Iin], j[Iin], k[Iin])
    for kk = 1:length(Itmp)
        ik = Iin[kk]
        v = S.SV[Itmp[kk]]
        if v == bsz[ik]
            lower[ik] = 1.
        elseif v == 0 || v == 2*bsz[ik]
            lower[ik] = 2.
        else
            lower[ik] = 0.5
        end
    end

    ## LEFT
    Iin = find((j-bsz) .>= 1)
    Itmp = sub2ind(S.sz, i[Iin], j[Iin]-bsz[Iin], k[Iin])
    for kk = 1:length(Itmp)
        ik = Iin[kk]
        v = S.SV[Itmp[kk]]
        if v == bsz[ik]
            left[ik] = 1.
        elseif v == 0
            left[ik] = 2.
        else
            left[ik] = 0.5
        end
    end

    ## RIGHT
    Iin = find((j+bsz .<= m2))
    Itmp = sub2ind(S.sz, i[Iin], j[Iin]+bsz[Iin], k[Iin])
    for kk = 1:length(Itmp)
        ik = Iin[kk]
        v = S.SV[Itmp[kk]]
        if v == bsz[ik]
            right[ik] = 1.
        elseif v == 0 || v == 2*bsz[ik]
            right[ik] = 2.
        else
            right[ik] = 0.5
        end
    end

    ## Front
    Iin = find((k-bsz) .>= 1)
    Itmp = sub2ind(S.sz, i[Iin], j[Iin], k[Iin]-bsz[Iin])
    for kk = 1:length(Itmp)
        ik = Iin[kk]
        v = S.SV[Itmp[kk]]
        if v == bsz[ik]
            front[ik] = 1.
        elseif v == 0
            front[ik] = 2.
        else
            front[ik] = 0.5
        end
    end

    ## Back
    Iin = find((k+bsz .<= m3))
    Itmp = sub2ind(S.sz, i[Iin], j[Iin], k[Iin]+bsz[Iin])
    for kk = 1:length(Itmp)
        ik = Iin[kk]
        v = S.SV[Itmp[kk]]
        if v == bsz[ik]
            back[ik] = 1.
        elseif v == 0 || v == 2*bsz[ik]
            back[ik] = 2.
        else
            back[ik] = 0.5
        end
    end

    return upper, lower, left, right, front, back
end
