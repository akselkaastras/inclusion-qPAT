%% Run
clc; clear; close all;
addpath(genpath(pwd));
%%
priortype = 'level';
noiseseed = [2,3,4,5];
noiselevel = [0.01,0.02,0.04,0.08,0.16];


for i=1:4
    for j=1:5
        seedstr = sprintf('%0.5g',noiseseed(i));
        noisestr = sprintf('%0.5g',noiselevel(j));
        
        %matlab_call = 'matlab -nodisplay -r ''cd("..");';
        matlab_call = 'matlab -nodisplay -r ''';
        call2 = 'addpath(genpath(pwd));';
        func_call = ['dataSave(',seedstr,',',noisestr,',"',priortype,'")'];
        end_call = ';exit''';
        
        cmd = [matlab_call, call2, func_call, end_call];
        
        submitToDTUCluster([priortype,num2str(i),num2str(j)],cmd)
    end
end