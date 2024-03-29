function fmdl = computeProjectionMatrices_coarse(fmdl,meshpar,meshpar_fine,priorpar,trunc)

%E_coarse = eigenbasisFEM(meshpar,trunc);
E_coarse = eigenbasisFEM_coarse(meshpar,meshpar_fine,trunc);
%s = load(strcat('Data/noise_model/eigenv/E_',num2str(meshpar.hmax),'.mat'));
%E_coarse = s.E;
E_coarse = E_coarse(:,1:trunc);

%if ~(strcmpi(priorpar.type,'starDG') || strcmpi(priorpar.type,'level'))
if ~strcmpi(priorpar.type,'starDG')
    fmdl.U_proj_coarse = fmdl.Carea*E_coarse;
else
    % Plane matrix
    fmdl.G = computePlanes(meshpar);
    [A,B,C] = computeProjMatrix(meshpar);
    fmdl.U_proj_coarse1 = A*E_coarse;
    fmdl.U_proj_coarse2 = B*E_coarse;
    fmdl.U_proj_coarse3 = C*E_coarse;

end