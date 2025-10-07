%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Imperial College London, United Kingdom
% Multiphase Systems Laboratory  
% Year:     2025
% MATLAB:   R2024a
% Authors:  Hassan Azzan (HA), Ayca Yilmaz (AY)
%
% Purpose:
% Run an example for kinetic batch adsorber analogue model (k-BAAM) 
%
% Last modified:
% - 2025-09-17, HA: Initial creation
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% clc; clear all; 
% close all;
addpath(genpath(pwd))
load('13X_AKR_16.mat') 
load('Zeolite13Xparams.mat') 
% load('CSAC_AKR_16.mat');
%% Parameters structure
% pre-exponent for site b, b0 [m3/mol]
% pre-exponent for site d, d0 [m3/mol]
% saturation capacity for site b, qsb [mol/kg]
% saturation capacity for site d, qsd [mol/kg]
% internal energy change for adsorption for site b, delUb [J/mol]
% internal energy change for adsorption for site d, delUd [J/mol]
% High pressure in the adsorption step, p_H [Pa]
% Intermediate pressure in the blowdown step, p_I [Pa]
% Low pressure in the evacuation step, p_L [Pa]
% Step durations, t_ads t_blo t_evac t_press [s]
% Temperature of gas, T_amb T_feed [K]
% Inlet molar flowrate for adsorption step, F_in [mol/s] 
% Inlet molar composition, y1_in [-] 
% Column volume, V_column [m3] 
% Bed voidage, e_bed [-] 
% Adsorbent particle density, Rho_s [-] 
% Time constant for vacuum pump, lambda [1/s]
% Adsorbent particle radius, rp [m]
% Type of pressurisation ("LPP" or "FP"), pressType [-]

%%
plotFlag = 1; % to plot profiles at CSS
theta = [1,   3e5,    300,  100,  500,    10e5]; % vector of decision variables, [F_in, P_I, t_ads, t_blo, t_evac]
% parameters.t_ads = 300;
parameters.t_press = 10;
% parameters.pressType = "FP";
parameters.heating = 0;
% parameters.y1_in = 0.0004;
parameters.p_L = 0.02e5;
% parameters.p_H = 8e5;
% parameters.Theat = 393;
parameters.outputType = "plot"; 
% parameters.outputType = "opt"; 
parameters.OptType = "Const";
% parameters.OptType = "Unc";
parameters.outputType = "plot";
parameters.processType = "PVSA";
parameters.adsorbentName = "Z13X";
parameters.pressType = "FP";
parameters.modelType = "nonisothermal";

%% Run to CSS and plot profiles
% KPIs = kBAAM_Outputs(parameters,plotFlag) % KPIS = [Recovery, Purity]

%% Run to CSS and output KPIs corresponding to theta
tic
KPIs = kBAAM_Outputs_isothermal(parameters,theta) % KPIS = [Recovery, Purity]
toc
tic
KPIs = kBAAM_Outputs_nonIsothermal(parameters,theta) % KPIS = [Recovery, Purity]
toc
