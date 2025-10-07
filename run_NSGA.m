%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Imperial College London, United Kingdom
% Multiphase Systems Laboratory
% Year:     2025
% MATLAB:   R2024a
% Authors:  Hassan Azzan (HA)
%
% Purpose:
% Function that takes parameters as inputs and carries out multiobjective
% process optimization for the system defined in parameters.
%
% Last modified:
% - 2025-09-21, HA: Initial creation
%
% Input arguments:
%   - parameters: contains adsorbent properties and process parameters
%
% Output arguments:
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function run_NSGA(parameters)
% Load PVSA model onto path
addpath(genpath(pwd))

if parameters.processType == "PVSA"
    lb = [0.1, log10(0.021e5),    10,    10,    10,     log10(1e5)];   % Lower bounds [v_in, p_I, t_ads, t_blo, t_evac, V_column, p_H]
    ub = [2,    log10(3e5),     200,  200,  200,     log10(10e5)];  % Upper bounds [v_in, p_I, t_ads, t_blo, t_evac, V_column, p_H]
    parameters.xRef = ub;
    % Linear inequality constraints: P_blo < P_ads
    A = [0, 1, 0, 0, 0,  -1];
    b = 0;
elseif parameters.processType == "VSA"
    lb = [0.1, log10(0.02e5),    5,   5,   5,  0.1];   % Lower bounds [v_in, p_I, t_ads, t_blo, t_evac, V_column, p_H]
    ub = [2,   log10(0.5e5),    600,  200,  600,  2.0];  % Upper bounds [v_in, p_I, t_ads, t_blo, t_evac, V_column, p_H]
    parameters.xRef = ub;
    % Linear inequality constraints: none
    A = [0, 0, 0, 0, 0 0];
    b = 0;
elseif parameters.processType == "AdsorbentVSA"
    lb = [log10(0.02e5), 100,  40,  0.5, 1e-6, 10e3, 10e3];  % Lower bounds [p_I, t_ads, t_blo, qsb, b0, delU1b, delU2b]
    ub = [log10(0.9e5),  3000, 300, 8,   1e-3, 45e3, 45e3];  % Upper bounds [p_I, t_ads, t_blo, qsb, b0, delU1b, delU2b]]
    parameters.xRef = ub;
    % Linear inequality constraints: |delU2| < |delU1|
    A = [0, 0, 0, 0, 0, -1, 1];
    b = 0;
elseif parameters.processType == "AdsorbentPVSA"
    lb = [log10(0.02e5), 100,  40,  0.5, 1e-6, 10e3, 10e3, log10(1e5)];   % Lower bounds [p_I, t_ads, t_blo, qsb, b0, delU1b, delU2b, p_H]
    ub = [log10(10e5),    3000, 300, 8,   1e-3, 45e3, 45e3, log10(10e5)];  % Upper bounds [p_I, t_ads, t_blo, qsb, b0, delU1b, delU2b, p_H]]
    parameters.xRef = ub;
    % Linear inequality constraints: |delU2| < |delU1| & P_blo < P_ads
    A = [1, 0, 0, 0, 0, -1, 1,-1];
    b = 0;
end

%% Solve optimisation problem using genetic algorithm
% Genetic algorithm settings
nVars = length(lb);
ngens = 90;     % Number of generations
pop_size = 300; % Number of members in population of each generation [-]
% pop_size = 30; % Number of members in population of each generation [-]

rng default
parameters.fileName = convertStringsToChars(strcat(parameters.adsorbentName,"_",parameters.processType,"_",parameters.OptType,"_",parameters.modelType,"_",datestr(now,'ddmmyyhhMM')));

p = sobolset(pop_size,'Skip',1e2,'Leap',1e1);
p = scramble(p,'MatousekAffineOwen');
X0 = net(p,nVars+2); X0 = X0';
X0 = X0(:,2:end-1);

X0 = lhsdesign(pop_size,nVars);
initPop = lb+X0.*ub;
parameters.outputType = "opt";
options = optimoptions(@gamultiobj, 'Display', 'iter', 'Generations', ngens, 'PopulationSize', pop_size,'UseParallel',true,'InitialPopulationMatrix',initPop,'CrossoverFraction',0.85,'ParetoFraction',0.35,'PlotFcn','gaplotpareto','MigrationInterval',5); % ,'PlotFcn','gaplotpareto' GA options
parameters.fileName = convertStringsToChars(strcat(parameters.adsorbentName,"_",parameters.processType,"_",parameters.OptType,"_",parameters.modelType,"_",datestr(now,'ddmmyyhhMM')));


[x, fval] = gamultiobj(@(theta) kBAAM_Outputs_nonIsothermal(parameters,theta), nVars, A, b, [], [], lb./parameters.xRef, ub./parameters.xRef, options); % Optimize objectives using GA

save(['matFiles/',parameters.fileName,'.mat'])
end