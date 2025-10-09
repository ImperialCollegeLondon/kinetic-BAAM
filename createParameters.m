%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Imperial College London, United Kingdom
% Multiphase Systems Laboratory
% Year:     2025
% MATLAB:   R2024a
% Authors:  Hassan Azzan (HA)
%
% Purpose:
% Creates and saves parameters structure for a given adsorbent that can be
% used for design and optimization
%
% Last modified:
% - 2025-10-09, HA: Add properties required for wall energy balance
% - 2025-10-08, HA: Initial creation
%
% Input arguments:
%
% Output arguments: 
% - parameters: structure containing information described below
%
% Dependencies:
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function parameters = createParameters
%% %%%  %%% Parameters structure  %%%  %%% 
%% Adsorbent, Adsorbate and bed properties 
% Name of the adsorbentt, adsorbentName [-]
% saturation capacity for site b, qsb [mol/kg]
% saturation capacity for site d, qsd [mol/kg]
% pre-exponent for site b, b0 [m3/mol]
% pre-exponent for site d, d0 [m3/mol]
% internal energy change for adsorption for site b, delUb [J/mol]
% internal energy change for adsorption for site d, delUd [J/mol]
% Gas specific heat capacity, cp_g [J/molK]
% Adsorbed phase specific heat capacity, cp_a [J/molK]
% Specific heat capacity of the column wall, cp_w [J/kgK]
% Adsorbent specific heat capacity, cp_s [J/kgK]
% Adsorbent particle density,rho_s [kg/m3] 
% Column wall density,rho_s [kg/m3] 
% Adsorbent particle radius, rp [m]
% Column volume, V_column [m3] 
% Column wall thickness, t_wall [m] 
% Bed voidage, e_bed [-] 
% Inside heat transfer coefficient for the bed (Haghpanah 2014) [W/m2K]
% Outside heat transfer coefficient for the bed (Haghpanah 2014) [W/m2K]
% Particle voidage, epsilon_p [-]
% Molecular diffusivity of gas, Dm [m2/s]
% Tortuosity, tau [-]
parameters.adsorbentName = "Z13X_AW_2022";
parameters.qsb_1 = 3.09;
parameters.qsd_1 = 2.54;
parameters.bo_1 = 8.65e-7;
parameters.do_1 = 2.63e-8;
parameters.delUb_1 = -3.66e4;
parameters.delUd_1 = -3.57e4;
parameters.qsb_2 = 5.84;
parameters.qsd_2 = 0;
parameters.bo_2 = 2.50e-7;
parameters.do_2 = 0;
parameters.delUb_2 = -1.58e4;
parameters.delUd_2 = 0;
parameters.cp_g = 30.7;
parameters.cp_a = parameters.cp_g;
parameters.cp_w = 502;
parameters.cp_s = 1070;
parameters.rho_s = 1130;
parameters.rho_w = 7800;
parameters.rp = 1e-3;
parameters.V_column = 0.066;
parameters.t_wall = 0.003;
parameters.e_bed = 0.37;
parameters.h_in = 8.6; 
parameters.h_out = 2.5; 
parameters.epsilon_p = 0.35; 
parameters.Dm = 1.6e-5;
parameters.tau = 3; 

%% Cycle Properties
% Flow velocity for adsorption step, v_in [m/s] 
% Inlet molar composition, y1_in [-] 
% Temperature of gas, T_feed [K]
% High pressure in the adsorption step, p_H [Pa]
% Intermediate pressure in the blowdown step, p_I [Pa]
% Low pressure in the evacuation step, p_L [Pa]
% Step durations, t_ads t_blo t_evac t_press [s]
% Heating in evacuation step (1 or 0) if temperature swing, heating [-]
% Type of pressurisation ("LPP" or "FP"), pressType [-]
% Type of cycle ("VSA" or "PSA"), processType [-]
% Time constant for vacuum pump, lambda [1/s]
parameters.v_in = 0.5;
parameters.y1_in = 0.15;
parameters.T_feed = 298;
parameters.p_H = 1.00e5;
parameters.p_I = 0.30e5;
parameters.p_L = 0.01e5;
parameters.t_ads = 300;
parameters.t_blo = 100;
parameters.t_evac = 300;
parameters.t_press = 20;
parameters.heating = 0;
parameters.pressType = "FP";
parameters.processType = "PVSA";
parameters.lambda = 0.5; 

%% Model Properties
% Type of model ("isothermal" or "nonisothermal"), modelType [-]
% Type of output if  ("p[ot" or "opt"), outputType [-]
% Type of optimization if outputType == "opt" ("isothermal" or "nonisothermal"), OptType [-]
parameters.modelType = "nonisothermal";
parameters.OptType = "Unc";
parameters.outputType = "plot";

% Save parameters structure in AdsorbentFfiles folder
save(['AdsorbentFiles/',convertStringsToChars(strcat(parameters.adsorbentName)),'.mat'])

end