%% Run
clc; clear; close all;
addpath(genpath(pwd));
%%
priortype = 'star';
dataSave(noiseseed, noiselevel, priortype)
funcCall = 'matlab -nodisplay -r ''cd(''FEM/auxfun''); functionname(argument1, argument2, argumentN), [0,1]);exit';
cmd = 'matlab -nodisplay -r curve_project';
submitToDTUCluster('1',cmd)