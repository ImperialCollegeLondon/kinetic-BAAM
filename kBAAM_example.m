clear all
clc
close all

%% Run to CSS and plot cycle corresponding to theta
% run createParameters.m first to generate the parameter struct
createParameters
load('13XH_T.mat')

parameters.outputType = "plot";
parameters.processType = "PVSA";
parameters.OptType = "Unc";
parameters.Lbyr = 7;
parameters.equilibrium = 0;
parameters.cCSTR = 0 ;
parameters.testBT = 0;
parameters.testEvac = 0;
parameters.normPlot = 0;
parameters.forwardEvac = 0;
parameters.pressureDrop = 1;
parameters.plotVideo = 0;
parameters.equilibrium = 0;
theta = [0.1 0.2e5  85 100 170 0.02e5 4e5]; % [v_in, p_I, t_ads, t_blo, t_evac, p_H]

tic
KPIs2 = kBAAM_Outputs_nonIsothermal_dP(parameters,theta); % KPIS = [-Recovery, -Purity];
kbaamDP = toc

