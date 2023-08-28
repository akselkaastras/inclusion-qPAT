function E = eigenbasisFEM(meshpar,trunc)

% Enters Dirichlet Laplacian eigenbasis into the columns of a matrix E
% We stop when we have trunc functions.


filename = strcat('Data/noise_model/eigenv/eigen_',num2str(meshpar.hmax),'_1.mat');
filename2 = strcat('Data/noise_model/eigenv/E_',num2str(meshpar.hmax),'.mat');

if ~isfile(filename) && ~isfile(filename2)
    findEigenLaplacian(meshpar);
end

if ~isfile(filename2)
    filename = strcat('Data/noise_model/eigenv/eigen_',num2str(meshpar.hmax));
    E = zeros(length(meshpar.p),trunc);
    k = 1;
    flag = 0;
    
    for i = 1:15
        str = strcat(filename,'_',num2str(i),'.mat');
        s = load(str);
        M = size(s.x,2);
        for j = 1:M
            E(:,k) = s.x(:,j);
            
            if k == trunc
                flag = 1;
                break;
            end
            k = k+1;
        end
        if flag == 1
            disp(' - - Truncation reached');
            break;
        end
        disp(strcat(' - - Iteration -- ', ' ', num2str(i), '-- done.'))
    
    end
    
    if flag == 0
        disp(' - - Error: Truncation was not reached')
        return;
    end

    filename = strcat('Data/noise_model/eigenv/E_',num2str(meshpar.hmax),'.mat');
    save(filename,'E');
end

filename = strcat('Data/noise_model/eigenv/eigen_',num2str(meshpar.hmax));
for i = 1:15
    str = strcat(filename,'_',num2str(i),'.mat');
    if isfile(str)
        delete(str);
    end
end

s = load(filename2);
E = s.E;



