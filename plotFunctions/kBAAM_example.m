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
clc; clear all; 
close all;
addpath(genpath(pwd))
% run createParameters.m first to generate the parameter struct
load('Z13X_AW_2022.mat') 

%% Run to CSS and plot profiles
tic
KPIs = kBAAM_Outputs_nonIsothermal(parameters) % KPIS = [-Recovery, -Purity]
toc
