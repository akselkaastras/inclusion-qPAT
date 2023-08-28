%% Run
clc; clear; close all;
addpath(genpath(pwd));
%% starDG
priortype = 'starDG';

iter = 1.5e6;
stepsize = [0.015,0.025,0.035,0.045,0.1];
noiselevel = [0.01,0.02,0.04,0.08,0.16];
sampleseed = [1,2,3,4,5];
samplemat = [0;10;20;30;40];
sampleseed = samplemat + sampleseed;
x0seed = 1;
noiseseed = [1,2,3,4,5];
std = [0.15,0.2,0.25,0.3,0.35];

for j=1:5
    for i=1:5
        iterstr = sprintf('%0.5g',iter);
        stepstr = sprintf('%0.5g',stepsize(i));
        noisestr = sprintf('%0.5g',noiselevel(i));
        x0seedstr = sprintf('%0.5g',x0seed);
        noiseseedstr = sprintf('%0.5g',noiseseed(j));
        sampleseedstr = sprintf('%0.5g',sampleseed(j,i));
        stdstr = sprintf('%0.5g',std(i));
        
        matlab_call = 'matlab -nodisplay -r ''';
        call2 = 'addpath(genpath(pwd));';
        func_call = ['starDGPCN(',iterstr,',',stepstr,',',noisestr,',',x0seedstr,',',noiseseedstr,',',sampleseedstr,',',stdstr,',13)'];
        end_call = ';exit''';
        
        cmd = [matlab_call, call2, func_call, end_call];
        
        submitToDTUCluster(['starDG_', num2str(i)],cmd)
    end
end

%% level

priortype = 'level';

iter = 2.2e6;
stepsize = [0.002, 0.003, 0.006, 0.01, 0.05];
noiselevel = [0.01,0.02,0.04,0.08,0.16];
sampleseed = [1,2,3,4,5];
samplemat = [10;20;30;40];
sampleseed = samplemat + sampleseed;
deltavec = [0.005,0.01,0.02,0.04,0.1];
x0seed = 1;
noiseseed = [2,3,4,5];

std = [0.75,1,1.5,2,2.5];
for j = 1:4
    for i = 1:5
        iterstr = sprintf('%0.5g',iter);
        stepstr = sprintf('%0.5g',stepsize(i));
        noisestr = sprintf('%0.5g',noiselevel(i));
        x0seedstr = sprintf('%0.5g',x0seed);
        noiseseedstr = sprintf('%0.5g',noiseseed(j));
        sampleseedstr = sprintf('%0.5g',sampleseed(j,i));
        stdstr = sprintf('%0.5g',std(i));
        deltastr = sprintf('%0.5g',deltavec(i));
        
        matlab_call = 'matlab -nodisplay -r ''';
        call2 = 'addpath(genpath(pwd));';
        func_call = ['levelPCN(',iterstr,',',stepstr,',',noisestr,',',x0seedstr,',',noiseseedstr,',',sampleseedstr,',',stdstr,',',deltastr,',13)'];
        end_call = ';exit''';
        
        cmd = [matlab_call, call2, func_call, end_call];
        
        submitToDTUCluster(['level_', num2str(i)],cmd)
    end
end