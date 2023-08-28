function fmdl = computeProjectionMatrices_coarse(fmdl,meshpar,priorpar,trunc)

s = load('Data/noise_model/eigenv/E_coarse_0.01_0.0175.mat');
E_coarse = s.E_coarse;
E_coarse = E_coarse(:,1:trunc);

if ~(strcmpi(priorpar.type,'starDG') || strcmpi(priorpar.type,'level'))
    fmdl.U_proj_coarse = fmdl.Carea*E_coarse;
else
    % Plane matrix
    fmdl.G = computePlanes(meshpar);
    [A,B,C] = computeProjMatrix(meshpar);
    fmdl.U_proj_coarse1 = A*E_coarse;
    fmdl.U_proj_coarse2 = B*E_coarse;
    fmdl.U_proj_coarse3 = C*E_coarse;

end