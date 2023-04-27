%% Run
clc; clear; close all;
addpath(genpath(pwd));
%%
priortype = 'star';
%noiseseed = [1,2,3,4,5];
%noiselevel = [0.005,0.01,0.02,0.04];

iter = 1000;
stepsize = 0.01;
noiselevel = 0.01;
x0seed = 1;
noiseseed = 1;

for i=1
    for j=1
        iterstr = sprintf('%0.5g',iter);
        stepstr = sprintf('%0.5g',stepsize);
        noisestr = sprintf('%0.5g',noiselevel);
        x0seedstr = sprintf('%0.5g',x0seed);
        noiseseedstr = sprintf('%0.5g',noiseseed);
        
        matlab_call = 'matlab -nodisplay -r ''';
        call2 = 'addpath(genpath(pwd));';
        func_call = ['starPCN(',iterstr,',',stepstr,',',noisestr,',',x0seedstr,',',noiseseedstr,')'];
        end_call = ';exit''';
        
        cmd = [matlab_call, call2, func_call, end_call];
        
        submitToDTUCluster('1',cmd)
    end
end