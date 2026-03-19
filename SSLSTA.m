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
function [q1_star, q2_star] = SSLSTA(P, y1, T, parameters)

isoparams1 = parameters.SSLSTA1;
isoparams2 = parameters.SSLSTA2;


qsNPA   =   isoparams1(1,1);
b01NPA    = isoparams1(2,1);
delU1NPA  = isoparams1(3,1);
qsLPA   =   isoparams1(4,1);
b01LPA    = isoparams1(5,1);
delU1LPA  = isoparams1(6,1);
kgate  = isoparams1(7,1);
cgate  = isoparams1(8,1);
sval  = isoparams1(9,1);


qsNPB   =   isoparams2(1,1);
b01NPB    = isoparams2(2,1);
delU1NPB  = isoparams2(3,1);
qsLPB   =   isoparams2(4,1);
b01LPB    = isoparams2(5,1);
delU1LPB  = isoparams2(6,1);

[q1_star, q2_star] = computeSSLSTALoadingBinary(P./(1e5),T,y1,qsNPA,qsLPA,b01NPA,b01LPA,delU1NPA,delU1LPA,qsNPB,qsLPB,b01NPB,b01LPB,delU1NPB,delU1LPB,kgate,cgate,sval);
end


