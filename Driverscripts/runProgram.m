%% Run
clc; clear; close all;
addpath(genpath(pwd));
%%
funcCall = 'cd(); functionname(argument1, argument2, argumentN);exit'
'matlab -nodisplay -r "cd('FEM/auxfun'); distToLine([2,3]', [3,4]', [0,1]);exit"'
cmd = 'matlab -nodisplay -r curve_project';
submitToDTUCluster('1',cmd)