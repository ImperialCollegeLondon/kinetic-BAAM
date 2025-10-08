%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Imperial College London, United Kingdom
% Multiphase Systems Laboratory
% Year:     2025
% MATLAB:   R2024a
% Authors:  Hassan Azzan (HA)
%
% Purpose:
% A routine for optimizing a PVSA process by minimizing over two conflicting
% objectives for different systems
%
% Last modified:
% - 2025-09-21, HA: Initial creation
%
% Input arguments:
%
% Output arguments:
%
% Dependencies:
%   - run_NSGA.m
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

addpath(genpath(pwd))

clear; close all; clc;

%%
load('Z13X_AW_2022.mat')
parameters.rp = 1e-3;
parameters.p_L = 0.02e5;
parameters.V_column = 0.0661;
parameters.outputType = "opt";
parameters.pressType = "FP";

%%
parameters.adsorbentName = "Z13X";
parameters.modelType = "nonisothermal";
parameters.OptType = "Const";
parameters.processType = "PVSA";
run_NSGA(parameters);

%%
parameters.adsorbentName = "Z13X";
parameters.modelType = "nonisothermal";
parameters.OptType = "Const";
parameters.processType = "VSA";
run_NSGA(parameters);

%%
parameters.adsorbentName = "Hypo";
parameters.modelType = "nonisothermal";
parameters.OptType = "Unc";
parameters.processType = "AdsorbentVSA";
run_NSGA(parameters);

%%
parameters.adsorbentName = "Hypo";
parameters.modelType = "nonisothermal";
parameters.OptType = "Unc";
parameters.processType = "AdsorbentPVSA";
run_NSGA(parameters);
