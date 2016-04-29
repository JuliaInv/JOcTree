function getLineIntegrationMatrix(M::OcTreeMesh, X1, X2)
#  Q = getLineIntegrationMatrix(mesh)
#  Integrate edge function along polygons defined in mesh.receivers.
#  Requires that polgon vertices are OcTree vertices and that they are
#  connected by edges of fine mesh cells.


S  = M.S.SV;
x0 = M.x0;
y0 = M.y0;
z0 = M.z0;
hx = M.h[1];
hy = M.h[2];
hz = M.h[3];

# get OcTree edge locations times two (so we have integer values)
EX, EY, EZ = getEdgeGrids(S);

# find the cell that corresponds to X1
tmp1 = findmin(abs(EX[:,1] - X1[1])); tmp1 = tmp1[2]
tmp2 = findmin(abs(EY[:,1] - X1[2])); tmp2 = tmp2[2]
tmp3 = findmin(abs(EZ[:,1] - X1[3])); tmp3 = tmp3[2] 



L = L(j); % edge length

q = D .* L;

Q = sparse(i,j,q,m,n);

Linv = sdiag(1 ./ sum(sparse(i, j, L), 2));

% for normalized measurement compute Linv * Q