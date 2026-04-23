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
addpath(genpath('../PVSA_opt-main'))

clear; close all; clc;

addpath(genpath(pwd))
addpath(genpath('../PVSA_opt-main'))

adsFiles = ["UTSA-16_VB_2019";
    "Mg-MOF-74_VB_2019";
    "Z13X_AW_2022";
    "CSAC_RH_2014"];

adsFiles = ["MgMOF74";
    "Z13X";
    "CSAC";
    "CALF20"];


adsFiles = ["13XH";
    "ZNaY";
    "ZHY";
    "ACRB3";
    "CALF20"];

for jj = 1:length(adsFiles)
    load(adsFiles(jj))
    %%
    parameters.rp = 1e-3;
    parameters.outputType = "opt";
    parameters.plot0D = 0;
    parameters.plotVideo = 0;
    parameters.layered =  0;
    parameters.Lbyr = 7;
    parameters.pressureDrop = 1;
    parameters.equilibrium = 0;
    parameters.cCSTR = 0;
    parameters.testBT = 0;
    parameters.testEvac =  0;
    parameters.normPlot = 0;
    parameters.plot0D = 0;
    parameters.rigid = 1;
    parameters.plotVideo = 0;
    parameters.layered =  0;
    parameters.adsorbentName = adsFiles(jj);
    parameters.modelType = "nonisothermal";


    parameters.pressType = "FP";
    
    parameters.processType = "PVSA";
    parameters.OptType = "Const"; 
    run_NSGA(parameters);
    parameters.OptType = "Unc";
    run_NSGA(parameters);


    parameters.processType = "VSA";

    parameters.OptType = "Const"; 
    run_NSGA(parameters);
    parameters.OptType = "Unc";
    run_NSGA(parameters);




    parameters.pressType = "LPP";
    
    parameters.processType = "PVSA";
    parameters.OptType = "Const"; 
    run_NSGA(parameters);
    parameters.OptType = "Unc";
    run_NSGA(parameters);


    parameters.processType = "VSA";

    parameters.OptType = "Const"; 
    run_NSGA(parameters);
    parameters.OptType = "Unc";
    run_NSGA(parameters);
end

% %%
% parameters.adsorbentName = "Z13X";
% parameters.modelType = "nonisothermal";
% parameters.OptType = "Const";
% parameters.processType = "PVSA";
% run_NSGA(parameters);
%
% %%
% parameters.adsorbentName = "Z13X";
% parameters.modelType = "nonisothermal";
% parameters.OptType = "Const";
% parameters.processType = "VSA";
% run_NSGA(parameters);
%
% %%
% parameters.adsorbentName = "Z13X"
% parameters.modelType = "nonisothermal";
% parameters.OptType = "Unc";
% parameters.processType = "VSA";
% run_NSGA(parameters);
%
% %%
% parameters.adsorbentName = "Z13X";
% parameters.modelType = "nonisothermal";
% parameters.OptType = "Unc";
% parameters.processType = "PVSA";
% run_NSGA(parameters);
% %%
% parameters.adsorbentName = "Hypo";
% parameters.modelType = "nonisothermal";
% parameters.OptType = "Unc";
% parameters.processType = "AdsorbentVSA";
% run_NSGA(parameters);
%
% %%
% parameters.adsorbentName = "Hypo";
% parameters.modelType = "nonisothermal";
% parameters.OptType = "Unc";
% parameters.processType = "AdsorbentPVSA";
% run_NSGA(parameters);

%%
% parameters.adsorbentName = "Resin_TLS";
% parameters.modelType = "nonisothermal";
% parameters.OptType = "Const";
% parameters.processType = "Resin";
% parameters.outputType = "opt";
% parameters.fixResins = 0;
% parameters.pressureDrop = 1;
%     parameters.equilibrium = 0;
%     parameters.cCSTR = 0;
%     parameters.testBT = 0;
%     parameters.testEvac =  0;
%     parameters.normPlot = 0;
% run_NSGA(parameters);
%%
% parameters.adsorbentName = "Resin_TLS";
% parameters.modelType = "nonisothermal";
% parameters.OptType = "Const";
% parameters.processType = "Resin";
% parameters.outputType = "opt";
% parameters.fixResins = 1;
% parameters.pressureDrop = 1;
%     parameters.equilibrium = 0;
%     parameters.cCSTR = 0;
%     parameters.testBT = 0;
%     parameters.testEvac =  0;
%     parameters.normPlot = 0;
%     parameters.plot0D = 0;
%     parameters.rigid = 1;
%     parameters.plotVideo = 0;
%     parameters.layered =  0;
% run_NSGA(parameters);
% % %%
% % parameters.adsorbentName = "Resin";
% % parameters.modelType = "nonisothermal";
% % parameters.OptType = "Unc";
% % parameters.processType = "Resin";
% % parameters.outputType = "opt";
% %
% % run_NSGA(parameters);
