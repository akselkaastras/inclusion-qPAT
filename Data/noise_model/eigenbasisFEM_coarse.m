function E_coarse = eigenbasisFEM_coarse(meshpar,meshpar_fine,trunc)
% Interpolates eigenbasis to coarse mesh.
% It was no good just recomputing the eigenvalues on the coarse mesh.
% Indeed, eigenfunctions were computed in a different order.
% A possible solution might be simple to reorder and sign normalize.

filename_coarse = strcat('Data/noise_model/eigenv/E_coarse_',num2str(meshpar.hmax),'.mat');

if isfile(filename_coarse)
    s = load(filename_coarse);
    E_coarse = s.E_coarse;
else
    filename = strcat('Data/noise_model/eigenv/E_',num2str(meshpar_fine.hmax),'.mat');
    
    if isfile(filename)
        s = load(filename);
        E = s.E;
    else
        disp('Error occured: no eigenbasis for the fine mesh found')
        return;
    end
    
    E = E(:,1:trunc);
    E_coarse = zeros(length(meshpar.p),trunc);

    disp('Interpolating eigenbasis from fine to coarse mesh');
    for i = 1:trunc
        disp(['Interpolating ',num2str(i),'/',num2str(trunc)])
        E_coarse(:,i) = interpolateMesh(E(:,i),meshpar.p(1,:)',meshpar.p(2,:)',meshpar_fine);
    end
    
    % Save
    save(filename_coarse,'E_coarse');
end


