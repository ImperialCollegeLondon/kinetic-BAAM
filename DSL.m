%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Imperial College London, United Kingdom
% Multiphase Systems Laboratory  
% Year:     2025
% MATLAB:   R2024a
% Authors:  Hassan Azzan (HA), Ayca Yilmaz (AY)
%
% Purpose:
% Function that takes pressure, mole fraction of component 1, temperature,
% and isotherm parameters. Outputs corresponding equilibrium adsorbed 
% amounts for both components.
%
% Last modified:
% - 2025-09-17, HA: Initial creation
%
% Input arguments: 
%   - P: pressure [Pa]
%   - y1: mole fraction of component 1 [-]
%   - T: temperature [K]
%   - DSL isotherm parameters
%
% Output arguments:
%   - q1_star: equilibrium adsorbed amount for component 1 [mol/kg]
%   - q2_star: equilibrium adsorbed amount for component 2 [mol/kg]
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [q1_star, q2_star] = DSL(P, y1, T, qsb_1, qsd_1, qsb_2, qsd_2, bo_1, do_1, bo_2, do_2, delUb_1, delUd_1, delUb_2, delUd_2)
R = 8.314;
b1 = bo_1 .* exp(-delUb_1 ./ (R .* T));
d1 = do_1 .* exp(-delUd_1 ./ (R .* T));
b2 = bo_2 .* exp(-delUb_2 ./ (R .* T));
d2 = do_2 .* exp(-delUd_2 ./ (R .* T));

q1_star = (qsb_1 .* b1 .* P .* y1 ./(R .* T) ) ./   (1 + b1 .* P .* y1 ./(R .* T) + b2 .* P .* (1 - y1)./(R .* T) ) + qsd_1 .* d1 .* P .* y1 ./(R .* T) ./ (1 + (d1  .* P .* y1 ./(R .* T) ) + d2  .* P .* (1 - y1) ./(R .* T) ) ;%mol/kg
q2_star = (qsb_2 .* b2 .* P .* (1 - y1) ./(R .* T) ) ./   (1 + b2 .* P .* (1 - y1) ./(R .* T) + b1 .* P .* y1./(R .* T) ) + qsd_2 .* d2 .* P .* (1- y1) ./(R .* T) ./ (1 + (d1  .* P .* y1 ./(R .* T) ) + d2  .* P .* (1 - y1) ./(R .* T) ) ;%mol/kg


end


