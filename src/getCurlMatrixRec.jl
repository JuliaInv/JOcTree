export getCurlMatrixRec, getCurlMatrix

function getCurlMatrix(M::OcTreeMeshFV)
# M.Curl = getCurlMatrix(M::OcTreeMesh)
     if isempty(M.Curl)
        M.Curl = getCurlMatrixRec(M)
     end
     return M.Curl
end


using MaxwellUtils.DiagTimesMTimesDiag

#function getCurlMatrixRec(S,h)
function getCurlMatrixRec(M)

# [CURL,T,ESZ,FSZ] = getCurlMatrix(S)

S = M.S; h = M.h

m1,m2,m3    = S.sz
FX,FY,FZ    = getFaceSize(M)
EX,EY,EZ    = getEdgeSize(M)
NFX,NFY,NFZ = getFaceNumbering(M)
NEX,NEY,NEZ = getEdgeNumbering(M)

######################################################################
######################################################################
##
## WORKING ON  X - F A C E S
##
######################################################################
######################################################################
# look for y edges next to x faces
i,j,k,fsz = find3(FX)
fn        = nonzeros(NFX)

##
# CREATE
#     DYZ : Y-EDGES --> X-FACES
# DYZ COMPUTES THE DERIVATIVE OF Y-EDGES INTO Z-DIRECTION AND
# DYZ : m1+1 x m2 x m3+1  -->  m1+1 x m2 x m3

front    = NEY.SV[sub2ind(NEY.sz,i,j,k),1];  front = vec(full(front))
back     = NEY.SV[sub2ind(NEY.sz,i,j,k+fsz),1];  back  = vec(full(back))
frontmid = zeros(Int64, length(front))
backmid  = zeros(Int64, length(back))

I = find( (j+div(fsz,2) .<= m2) & (fsz .>= 2) )
if ~isempty(I)
  frontmid[I] = NEY.SV[sub2ind(NEY.sz, i[I] , j[I]+div(fsz[I],2) , k[I]  ),1]
  backmid[I]  = NEY.SV[sub2ind(NEY.sz, i[I] , j[I]+div(fsz[I],2) , k[I]+fsz[I] ),1]
end
# front y

If1 = find(frontmid .== 0) # SINGLE FRONTSIDE EDGE PER FACE
If2 = find(frontmid .>  0) # 2 FRONT EDGES PER FACE
Ib1 = find(backmid .== 0)  # SINGLE BACKSIDE EDGE PER FACE
Ib2 = find(backmid .>  0)  # 2 BACKSIDE EDGES PER FACE

ii = [fn[If1];       fn[If2];           fn[If2];            fn[Ib1];        fn[Ib2];            fn[Ib2]          ]
jj = [front[If1];    front[If2];        frontmid[If2];      back[Ib1];      back[Ib2];          backmid[Ib2]     ]
vv = [-ones(length(If1));  -ones(length(If2));  -ones(length(If2)); ones(length(Ib1));    ones(length(Ib2));    ones(length(Ib2))]

DYZ = sparse(ii, jj, vv,nnz(NFX),nnz(NEY))


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# CREATE
#     DZY : Z-EDGES --> X-FACES
# DZY COMPUTES THE DERIVATIVE OF Z-EDGES INTO Y-DIRECTION AND
# THEREFOR DZY : m1+1 x m2+1 x m3  -->  m1+1 x m2 x m3
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

left     = NEZ.SV[sub2ind(NEZ.sz,i,j,k    ),1];  left  = vec(full(left))
right    = NEZ.SV[sub2ind(NEZ.sz,i,j+fsz,k),1];  right = vec(full(right))
leftmid  = zeros(Int64, length(left))
rightmid = zeros(Int64, length(right))

I = find((j+div(fsz,2) .<= m2) & (fsz .>= 2) )
if ~isempty(I)
    leftmid[I]  = NEZ.SV[sub2ind(NEZ.sz, i[I], j[I], k[I]+div(fsz[I],2) ),1]
    rightmid[I] = NEZ.SV[sub2ind(NEZ.sz, i[I], j[I]+fsz[I], k[I]+div(fsz[I],2)),1]
end
Il1 = find(leftmid .== 0)  # SINGLE LEFT EDGE PER FACE
Il2 = find(leftmid .>  0)  # 2 LEFT EDGES PER FACE
Ir1 = find(rightmid .== 0) # SINGLE RIGHT EDGE PER FACE
Ir2 = find(rightmid .>  0) # 2 RIGHT EDGES PER FACE

ii = [fn[Il1];  fn[Il2];  fn[Il2];  fn[Ir1];  fn[Ir2];  fn[Ir2]   ]
jj = [left[Il1];  left[Il2];   leftmid[Il2];  right[Ir1];  right[Ir2];  rightmid[Ir2] ]
vv = [-ones(length(Il1)); -ones(length(Il2));  -ones(length(Il2));  ones(length(Ir1)); ones(length(Ir2));  ones(length(Ir2))  ]

DZY = sparse(ii, jj, vv,nnz(NFX),nnz(NEZ))



######################################################################
######################################################################
##
## WORKING ON  Y - F A C E S
##
######################################################################
######################################################################
# look for y edges next to x faces
i,j,k,fsz = find3(FY)
fn        = nonzeros(NFY)




#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# CREATE
#     DXZ : X-EDGES --> Y-FACES
# DXZ COMPUTES THE DERIVATIVE OF X-EDGES INTO Z-DIRECTION AND
# DXZ : m1 x m2+1 x m3+1  -->  m1 x m2+1 x m3

front    = NEX.SV[sub2ind(NEX.sz,i,j,k    ),1];  front = vec(full(front))
back     = NEX.SV[sub2ind(NEX.sz,i,j,k+fsz),1];  back  = vec(full(back));
frontmid = zeros(Int64, length(front))
backmid  = zeros(Int64, length(back))

I = find((i+div(fsz,2) .<= m1) & (fsz .>= 2) )
if ~isempty(I)
    frontmid[I] = NEX.SV[sub2ind(NEX.sz, i[I]+div(fsz[I],2) , j[I] , k[I]),1]
    backmid[I]  = NEX.SV[sub2ind(NEX.sz, i[I]+div(fsz[I],2) , j[I] , k[I]+fsz[I]),1]
end

If1 = find(frontmid .== 0)  # SINGLE FRONT EDGE PER FACE
If2 = find(frontmid .>  0)  # 2 FRONT EDGES PER FACE
Ib1 = find(backmid .== 0)   # SINGLE BACK EDGE PER FACE
Ib2 = find(backmid .>  0)   # 2 BACK EDGES PER FACE

ii = [fn[If1];        fn[If2];           fn[If2];           fn[Ib1];        fn[Ib2];          fn[Ib2]          ]
jj = [front[If1];     front[If2];        frontmid[If2];     back[Ib1];      back[Ib2];        backmid[Ib2]     ]
vv = [-ones(length(If1));   -ones(length(If2));  -ones(length(If2));  ones(length(Ib1));    ones(length(Ib2));  ones(length(Ib2))];

DXZ = sparse(ii, jj, vv,nnz(NFY),nnz(NEX))


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# CREATE
#     DZX : Z-EDGES --> Y-FACES
# DZXCOMPUTES THE DERIVATIVE OF Z-EDGES INTO X-DIRECTION AND
# THEREFOR DZX : m1+1 x m2+1 x m3  -->  m1 x m2+1 x m3
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

upper    = NEZ.SV[sub2ind(NEZ.sz,i    ,j,k),1];  upper = vec(full(upper))
lower    = NEZ.SV[sub2ind(NEZ.sz,i+fsz,j,k),1];  lower = vec(full(lower))
uppermid = zeros(Int64, length(upper))
lowermid = zeros(Int64, length(lower))

I = find((i+div(fsz,2) .<= m1) & (fsz .>= 2) )
if ~isempty(I)
    uppermid[I] = NEZ.SV[sub2ind(NEZ.sz, i[I], j[I], k[I]+div(fsz[I],2) ),1]
    lowermid[I] = NEZ.SV[sub2ind(NEZ.sz, i[I]+fsz[I] , j[I] , k[I]+div(fsz[I],2) ),1]
end
Iu1 = find(uppermid .== 0)  # SINGLE UPPER EDGE PER FACE
Iu2 = find(uppermid .>  0)  # 2 UPPER EDGES PER FACE
Il1 = find(lowermid .== 0)  # SINGLE LOWER EDGE PER FACE
Il2 = find(lowermid .>  0)  # 2 LOWER EDGES PER FACE

ii = [fn[Iu1];        fn[Iu2];           fn[Iu2];           fn[Il1];        fn[Il2];          fn[Il2]          ]
jj = [upper[Iu1];     upper[Iu2];        uppermid[Iu2];     lower[Il1];     lower[Il2];       lowermid[Il2]    ]
vv = [-ones(length(Iu1));   -ones(length(Iu2));  -ones(length(Iu2));  ones(length(Il1)); ones(length(Il2)); ones(length(Il2))]

DZX = sparse(ii, jj, vv,nnz(NFY),nnz(NEZ))

######################################################################
######################################################################
#
# WORKING ON  Z - F A C E S
#
######################################################################
######################################################################
# look for y edges next to x faces
i,j,k,fsz = find3(FZ)
fn        = nonzeros(NFZ)


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# CREATE
#     DXY : X-EDGES --> Z-FACES
# DXY COMPUTES THE DERIVATIVE OF X-EDGES INTO Y-DIRECTION AND
# THEREFOR DXY : m1 x m2+1 x m3+1  -->  m1 x m2 x m3+1
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

left     = NEX.SV[sub2ind(NEX.sz,i,j,k    ),1];  left  = vec(full(left))
right    = NEX.SV[sub2ind(NEX.sz,i,j+fsz,k),1];  right = vec(full(right))
leftmid  = zeros(Int64, length(left))
rightmid = zeros(Int64, length(right))

I = find((i+div(fsz,2) .<= m1) & (fsz .>= 2) )

if ~isempty(I)
    leftmid[I]  = NEX.SV[sub2ind(NEX.sz, i[I]+div(fsz[I],2), j[I]        , k[I] ),1]
    rightmid[I] = NEX.SV[sub2ind(NEX.sz, i[I]+div(fsz[I],2), j[I]+fsz[I] , k[I] ),1]
end
Il1 = find(leftmid .== 0)  # SINGLE LEFT EDGE PER FACE
Il2 = find(leftmid .>  0)  # 2 LEFT EDGES PER FACE
Ir1 = find(rightmid .== 0) # SINGLE RIGHT EDGE PER FACE
Ir2 = find(rightmid .>  0) # 2 RIGHT EDGES PER FACE

ii = [fn[Il1];        fn[Il2];           fn[Il2];           fn[Ir1];        fn[Ir2];          fn[Ir2]          ];
jj = [left[Il1];      left[Il2];         leftmid[Il2];      right[Ir1];     right[Ir2];       rightmid[Ir2]    ];
vv = [-ones(length(Il1));   -ones(length(Il2));  -ones(length(Il2));  ones(length(Ir1));    ones(length(Ir2));  ones(length(Ir2))]

DXY = sparse(ii, jj, vv,nnz(NFZ),nnz(NEX))

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#% CREATE
#%     DYX : Y-EDGES --> Z-FACES
#% DYX COMPUTES THE DERIVATIVE OF Y-EDGES INTO X-DIRECTION AND
#% THEREFOR DYX : m1+1 x m2 x m3+1  -->  m1 x m2 x m3+1
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

upper    = NEY.SV[sub2ind(NEY.sz,i    ,j,k),1];  upper = vec(full(upper))
lower    = NEY.SV[sub2ind(NEY.sz,i+fsz,j,k),1];  lower = vec(full(lower))
uppermid = zeros(Int64, length(upper))
lowermid = zeros(Int64, length(lower))

I = find((j+div(fsz,2) .<= m2) & (fsz .>= 2) )
if ~isempty(I)
    uppermid[I] = NEY.SV[sub2ind(NEY.sz, i[I]        , j[I]+div(fsz[I],2) , k[I] ),1]
    lowermid[I] = NEY.SV[sub2ind(NEY.sz, i[I]+fsz[I] , j[I]+div(fsz[I],2) , k[I] ),1]
end
Iu1 = find(uppermid .== 0)  # SINGLE UPPER EDGE PER FACE
Iu2 = find(uppermid .>  0)  # 2 UPPER EDGES PER FACE
Il1 = find(lowermid .== 0)  # SINGLE LOWER EDGE PER FACE
Il2 = find(lowermid .>  0)  # 2 LOWER EDGES PER FACE

ii = [fn[Iu1];        fn[Iu2];           fn[Iu2];           fn[Il1];        fn[Il2];          fn[Il2]          ];
jj = [upper[Iu1];     upper[Iu2];        uppermid[Iu2];     lower[Il1];     lower[Il2];       lowermid[Il2]    ];
vv = [-ones(length(Iu1));   -ones(length(Iu2));  -ones(length(Iu2));  ones(length(Il1));   ones(length(Il2));  ones(length(Il2))];


DYX = sparse(ii, jj, vv,nnz(NFZ),nnz(NEY))

nedges = nnz(EX) + nnz(EY) + nnz(EZ)
nfaces = nnz(FX) + nnz(FY) + nnz(FZ)

#ESZ  = sdiag([nonzeros(EX)*h[1];nonzeros(EY)*h[2];nonzeros(EZ)*h[3]])
ESZ  = vcat(nonzeros(EX)*h[1], nonzeros(EY)*h[2], nonzeros(EZ)*h[3])

#FSZi = sdiag(1./[nonzeros(FX).^2*h[2]*h[3];
#                 nonzeros(FY).^2*h[1]*h[3];
#                 nonzeros(FZ).^2*h[1]*h[2]])
FSZi = 1.0 ./ vcat(nonzeros(FX).^2*(h[2]*h[3]),
                   nonzeros(FY).^2*(h[1]*h[3]),
                   nonzeros(FZ).^2*(h[1]*h[2]) )


ZEROEX = spzeros(nnz(NFX), nnz(NEX))
ZEROEY = spzeros(nnz(NFY), nnz(NEY))
ZEROEZ = spzeros(nnz(NFZ), nnz(NEZ))

T = [ ZEROEX   -DYZ       DZY;
         DXZ    ZEROEY   -DZX;
        -DXY    DYX     ZEROEZ]

#CURL = FSZi * T * ESZ
CURL = DiagTimesMTimesDiag(FSZi, T, ESZ)

return CURL
end  # function getCurlMatrixRec
