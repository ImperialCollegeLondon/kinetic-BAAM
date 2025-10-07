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
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


addpath(genpath(pwd))

clear; close all; clc;

% load('13X_AKR_16.mat')
load('13X_AW_22_2.mat')

parameters.rp = 1e-3;
parameters.p_L = 0.02e5;
parameters.V_column = 0.0661;
parameters.e_bed = 0.37;
parameters.t_press = 20;
parameters.outputType = "opt";
parameters.pressType = "LPP";
%%
parameters.adsorbentName = "Z13X";
parameters.modelType = "nonisothermal";
parameters.OptType = "Const";
parameters.processType = "PVSA";
run_NSGA(parameters);
% run_sobol(parameters);c

%%
parameters.adsorbentName = "Z13X";
parameters.modelType = "isothermal";
parameters.OptType = "Const";
parameters.processType = "PVSA";
run_NSGA(parameters);
% run_sobol(parameters);

%%
parameters.adsorbentName = "Hypo";
parameters.modelType = "nonisothermal";
parameters.OptType = "Unc";
parameters.processType = "AdsorbentPVSA";
parameters.y1_in = 0.01;

run_NSGA(parameters);
% run_sobol(parameters);

%%
parameters.adsorbentName = "Hypo";
parameters.modelType = "isothermal";
parameters.OptType = "Unc";
parameters.processType = "AdsorbentVSA";

run_NSGA(parameters);
% run_sobol(parameters);
