priortype = 'level';
noiseseed = [1,2,3,4,5];
noiselevel = [0.005,0.01,0.02,0.04];
datapar.meshpar_fine.hmax = 0.01;
datapar.meshpar.hmax = 0.02;

for i = 1:length(noiseseed)
    for j = 1:length(noiselevel)
        noise = noiselevel(j);
        seed = noiseseed(i);
        filename = strcat(num2str(datapar.meshpar_fine.hmax),'_',num2str(datapar.meshpar.hmax),'_',priortype,'_',num2str(noise), '.mat');
        lol = strcat('Data/sigma2/sigma2_',filename);
        load(lol)
        sigma.m = mean(sigma.V,2);
        C = cov(sigma.V');
        sigma.sigmasq_m = trace(C)/length(sigma.m);
        save(lol,'sigma');
        filename = strcat(num2str(datapar.meshpar_fine.hmax),'_',num2str(datapar.meshpar.hmax),'_',priortype,'_',num2str(noise), '.mat');
        lol = strcat('Data/data/data_','noiseseed_',num2str(seed),'_',filename);
        load(lol)
        datapar.m = m;
        datapar.sigmasq_m = sigma.sigmasq_m;
        datapar.epssq_approx = datapar.epssq + datapar.sigmasq_m;
        save(lol,'datapar');
    end
end