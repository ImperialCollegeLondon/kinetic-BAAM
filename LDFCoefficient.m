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
% corresponding equilibrium adsorbed amounts, and parameters. Outputs
% concentration dependent mass transfer coefficients
%
% Last modified:
% - 2025-09-17, HA: Initial creation
%
% Input arguments:
%   - P: pressure [Pa]
%   - y1: mole fraction of component 1 [-]
%   - T: temperature [K]
%   - q1_star: equilibrium adsorbed amount of component 1 at P, y1, T
%   - q2_star: equilibrium adsorbed amount of component 2 at P, y1, T
%   - parameters: contains adsorbent properties and process parameters
%
% Output arguments:
%   - k1: concentration dependent mass transfer coefficient for component 1 [1/s]
%   - k2: concentration dependent mass transfer coefficient for component 2 [1/s]
%
% Dependencies:
%  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [k1, k2] = LDFCoefficient(P,y1,T,q1_star,q2_star,parameters)
R = 8.314;

rp = parameters.rp; % Particle radius [m]
Dp = parameters.Dm/parameters.tau; % Effective pore diffusivity [m2/s]

k01 = (15*parameters.epsilon_p*Dp/(rp^2)); % Adsorption rate constant for component 1 [/s]
k02 = (15*parameters.epsilon_p*Dp/(rp^2)); % Adsorption rate constant for component 1 [/s]

if y1 > 0 && y1 < 1
    ratio1 = ((P.*y1./(R.*T))./q1_star)./parameters.rho_s; % c1/q1_Star*rho_s (mol/m^3)
    ratio2 = ((P.*(1-y1)./(R.*T))./q2_star)./parameters.rho_s; % c2/q2_Star*rho_s (mol/m^3)
elseif y1 <= 0
    ratio1 = ((P.*(1-y1)./(R.*T))./q2_star)./parameters.rho_s;
    ratio2 = 1; % c2/q2_Star*rho_s (mol/m^3)
elseif y1 >= 1
    ratio1 = ((P.*y1./(R.*T))./q1_star)./parameters.rho_s; % c1/q1_Star*rho_s (mol/m^3)
    ratio2 = 1; % to avoid division by 0
else
    ratio1 = 1; % c1/q1_Star*rho_s (mol/m^3)
    ratio2 = 1; % to avoid division by 0
end
if isnan(ratio1)
    ratio1 = 1;
end
if isnan(ratio2)
    ratio2 = 1;
end
k1 = k01.*ratio1; % mass transfer coefficient of CO2
k2 = k02.*ratio2; % mass transfer coefficient of N2
end