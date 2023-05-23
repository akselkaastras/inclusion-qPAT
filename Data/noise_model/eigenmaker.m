%% Run
clc; clear; close all;
addpath(genpath(pwd));
%%



    
matlab_call = 'matlab -nodisplay -r ''';
call2 = 'addpath(genpath(pwd)); meshpar = mesh_comp(0.01);';
func_call = 'findEigenLaplacian(meshpar)';
end_call = ';exit''';

cmd = [matlab_call, call2, func_call, end_call];

submitToDTUCluster(num2str(1),cmd)
