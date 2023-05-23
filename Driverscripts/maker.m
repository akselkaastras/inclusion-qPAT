%% Run
clc; clear; close all;
addpath(genpath(pwd));
%% starDG
priortype = 'starDG';

iter = 1e6;
stepsize = [0.02,0.025,0.05,0.1,0.2];
noiselevel = [0.01,0.02,0.04,0.08,0.16];
x0seed = 1;
noiseseed = 1;
std = [0.15,0.2,0.25,0.3,0.35];

for i=1:5
    iterstr = sprintf('%0.5g',iter);
    stepstr = sprintf('%0.5g',stepsize(i));
    noisestr = sprintf('%0.5g',noiselevel(i));
    x0seedstr = sprintf('%0.5g',x0seed);
    noiseseedstr = sprintf('%0.5g',noiseseed);
    stdstr = sprintf('%0.5g',std(i));
    
    matlab_call = 'matlab -nodisplay -r ''';
    call2 = 'addpath(genpath(pwd));';
    func_call = ['starDGPCN(',iterstr,',',stepstr,',',noisestr,',',x0seedstr,',',noiseseedstr,',',stdstr,')'];
    end_call = ';exit''';
    
    cmd = [matlab_call, call2, func_call, end_call];
    
    submitToDTUCluster(['starDG_', num2str(i)],cmd)
end

%% level

priortype = 'level';

iter = 1.5e6;
stepsize = [0.005, 0.006, 0.01, 0.05, 0.1];
noiselevel = [0.01,0.02,0.04,0.08,0.16];
x0seed = 1;
noiseseed = 1;
std = [0.95,1,1.05,1.1,1.15];

for i=1:5
    iterstr = sprintf('%0.5g',iter);
    stepstr = sprintf('%0.5g',stepsize(i));
    noisestr = sprintf('%0.5g',noiselevel(i));
    x0seedstr = sprintf('%0.5g',x0seed);
    noiseseedstr = sprintf('%0.5g',noiseseed);
    
    matlab_call = 'matlab -nodisplay -r ''';
    call2 = 'addpath(genpath(pwd));';
    func_call = ['levelPCN(',iterstr,',',stepstr,',',noisestr,',',x0seedstr,',',noiseseedstr,',',stdstr,')'];
    end_call = ';exit''';
    
    cmd = [matlab_call, call2, func_call, end_call];
    
    submitToDTUCluster(['level_', num2str(i)],cmd)
end