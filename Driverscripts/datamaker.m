%% Run
clc; clear; close all;
addpath(genpath(pwd));
%%
priortype = 'star';
noiseseed = [1,2,3,4,5];
noiselevel = [0.005,0.01,0.02,0.04];


for i=1:length(noiseseed)
    for j=1:length(noiselevele)
        seedstr = sprintf('%0.5g',noiseseed);
        noisestr = sprintf('%0.5g',noiselevel);
        priorstr = sprintf('''%s''',priortype);
        
        matlab_call = 'matlab -nodisplay -r ''cd(''Data/dataSave'');';
        func_call = ['dataSave(',seedstr,',',noisestr,',',priortype,')'];
        end_call = ';exit''';
        
        cmd = [matlab_call, func_call, end_call];
        
        submitToDTUCluster([num2str(i),num2str(j)],cmd)
    end
end