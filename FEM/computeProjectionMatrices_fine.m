function fmdl = computeProjectionMatrices_fine(fmdl,meshpar,trunc)
E = eigenbasisFEM(meshpar,trunc);


% Plane matrix
fmdl.G = computePlanes(meshpar);
[A,B,C] = computeProjMatrix(meshpar);
fmdl.U_proj1 = A*E;
fmdl.U_proj2 = B*E;
fmdl.U_proj3 = C*E;

fmdl.E = E;