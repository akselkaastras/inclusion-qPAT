%% Run
clc; clear; close all;
addpath(genpath(pwd));
%%
priortype = 'star';
%noiseseed = [1,2,3,4,5];
%noiselevel = [0.005,0.01,0.02,0.04];

iter = 1*1e6;
stepsize = [0.004,0.006,0.008,0.01,0.03];
noiselevel = [0.01,0.02,0.04,0.08,0.16];
x0seed = 1;
noiseseed = 1;

for i=1:5
    iterstr = sprintf('%0.5g',iter);
    stepstr = sprintf('%0.5g',stepsize(i));
    noisestr = sprintf('%0.5g',noiselevel(i));
    x0seedstr = sprintf('%0.5g',x0seed);
    noiseseedstr = sprintf('%0.5g',noiseseed);
    
    matlab_call = 'matlab -nodisplay -r ''';
    call2 = 'addpath(genpath(pwd));';
    func_call = ['starPCN(',iterstr,',',stepstr,',',noisestr,',',x0seedstr,',',noiseseedstr,')'];
    end_call = ';exit''';
    
    cmd = [matlab_call, call2, func_call, end_call];
    
    submitToDTUCluster(num2str(i),cmd)
end