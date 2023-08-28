meshpar = mesh_comp(0.01);
meshpar_coarse = mesh_comp(0.0175);
N = 15;
trunc = (2*N+1)*N;
E = eigenbasisFEM(meshpar,trunc);
p1 = meshpar_coarse.p(1,:);
p2 = meshpar_coarse.p(2,:);
E_coarse = zeros(length(p1),size(E,2));

for i = 1:size(E,2)
    vq = interpolateMesh(E(:,i),p1',p2',meshpar);
    E_coarse(:,i) = vq;
    i
end