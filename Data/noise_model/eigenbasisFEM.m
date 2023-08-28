function E = eigenbasisFEM(meshpar,trunc)

hmax = meshpar.hmax;
filename = strcat('Data/noise_model/eigenv/eigen_',num2str(hmax));




% Reading again and adding to noise expansion
E = zeros(length(meshpar.p),trunc);
k = 1;
flag = 0;
for i = 1:30
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


