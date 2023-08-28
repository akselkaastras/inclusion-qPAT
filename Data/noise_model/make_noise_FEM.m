function [xi, norm_noise, norm_noise_fem, D] = make_noise_FEM(meshpar,fmdl,trunc)

hmax = meshpar.hmax;
filename = strcat('Data/noise_model/eigenv/eigen_',num2str(hmax));
m = 0;
for i = 1:55
    str = strcat(filename,'_',num2str(i),'.mat');
    s = load(str);
    m = m + size(s.x,2);
end
disp(strcat('Found +','  ', num2str(m),' eigenvectors for hmax = ', num2str(hmax)));


% Noise expansion coefficients
e = randn(1,trunc);
rng(0);

% Reading again and adding to noise expansion
xi = zeros(length(meshpar.p),1);
D = zeros(m,1);
k = 1;
flag = 0;
for i = 1:55
    str = strcat(filename,'_',num2str(i),'.mat');
    s = load(str);
    M = size(s.x,2);
    for j = 1:M
        ee = randn;
        e(k) = ee;
        xi = xi + ee*s.x(:,j);
        D(j) = s.y(j);
        k = k+1;
        if k == trunc
            flag = 1;
            break;
        end
    end
    if flag == 1
        disp(' - - Truncation reached');
        break;
    end
    disp(strcat(' - - Iteration -- ', ' ', num2str(i), '-- done.'))
end



norm_noise = sum(e.^2);
norm_noise_fem = xi'*fmdl.Carea*xi;
