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
% - 2025-10-08, HA: Add reverse engineering method
% - 2025-09-21, HA: Initial creation
%
% Input arguments:
%   - parameters: contains adsorbent properties and process parameters
%
% Output arguments:
%
% Dependencies:
%   - kBAAM_Outputs_nonIsothermal.m
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function run_NSGA(parameters)

% Load model onto path
addpath(genpath(pwd))

if parameters.processType == "PVSA"  % Optimization for fixed adsorbent defined in parameters using PVSA
    lb = [0.3.*0.37,  (0.13e5),    40,    30,    30,   (1e5)];   % Lower bounds [v_in, p_I, t_ads, t_blo, t_evac]
    ub = [3.*0.37,    (5e5),     300,  300,  300,  (10e5)];  % Upper bounds [v_in, p_I, t_ads, t_blo, t_evac]

    % lb = [0.1.*0.37,  (0.13e5),    20,    30,    30,  (0.02e5),  (1e5)];   % Lower bounds [v_in, p_I, t_ads, t_blo, t_evac]
    % ub = [3.*0.37,    (3e5),     300,  300,  300,  (0.5e5),  (10e5)];  % Upper bounds [v_in, p_I, t_ads, t_blo, t_evac]


    parameters.xRef = ones(1,length(ub));
    % Linear inequality constraints: P_blo < P_ads    
    %                                P_evac < P_blo
    A = [0, -1, 0, 0, 0, 0 ;
        0, 1, 0, 0, 0,  -1 ];
    b = [0;0];
elseif parameters.processType == "VSA"  % Optimization for fixed adsorbent defined in parameters using VSA

    if parameters.amine
        lb = [0.3.*0.37,  (0.021e5),    100,    20,    30,  (0.02e5)];   % Lower bounds [v_in, p_I, t_ads, t_blo, t_evac]
        ub = [3.*0.37,    (0.9e5),     3e4,  1200,  3.5e4,  (0.5e5)];  % Upper bounds [v_in, p_I, t_ads, t_blo, t_evac]
    else
        lb = [0.3.*0.37,  (0.021e5),    30,    30,    30,  (0.02e5)];   % Lower bounds [v_in, p_I, t_ads, t_blo, t_evac]
        ub = [3.*0.37,    (0.9e5),     200,  300,  300,  (0.5e5)];  % Upper bounds [v_in, p_I, t_ads, t_blo, t_evac]
    end
    parameters.xRef = ub;
    % Linear inequality constraints: none
    A = [0, -1, 0, 0, 0,1];
    b = 0;
elseif parameters.processType == "AdsorbentVSA" % Reverse engineer adsorbent assuming equal pre-exponent for each gas using VSA
    lb = [log10(0.02e5), 100,  40,  0.5, 1e-6, 10e3, 10e3];  % Lower bounds [p_I, t_ads, t_blo, qsb, b0, delU1b, delU2b]
    ub = [log10(0.9e5),  3000, 300, 8,   1e-3, 45e3, 45e3];  % Upper bounds [p_I, t_ads, t_blo, qsb, b0, delU1b, delU2b]]
    parameters.xRef = ub;
    % Linear inequality constraints: |delU2| < |delU1|
    A = [0, 0, 0, 0, 0, -1, 1];
    b = 0;
elseif parameters.processType == "AdsorbentPVSA" % Reverse engineer adsorbent assuming equal pre-exponent for each gas using PVSA
    lb = [log10(0.02e5), 100,  40,  0.5, 1e-6, 10e3, 10e3, log10(1e5)];   % Lower bounds [p_I, t_ads, t_blo, qsb, b0, delU1b, delU2b, p_H]
    ub = [log10(10e5),    3000, 300, 8,   1e-3, 45e3, 45e3, log10(10e5)];  % Upper bounds [p_I, t_ads, t_blo, qsb, b0, delU1b, delU2b, p_H]]
    parameters.xRef = ub;
    % Linear inequality constraints: |delU2| < |delU1| & P_blo < P_ads
    A = [1, 0, 0, 0, 0, -1, 1,-1];
    b = 0;
elseif parameters.processType == "AdsorbentVSAb0" % Reverse engineer adsorbent assuming equal deltaU for each gas using VSA
    lb = [log10(0.02e5), 100,  40,  0.5, 1e-6, 10e3, 10e3];  % Lower bounds [p_I, t_ads, t_blo, qsb, b0, delU1b, delU2b]
    ub = [log10(0.9e5),  3000, 300, 8,   1e-3, 45e3, 45e3];  % Upper bounds [p_I, t_ads, t_blo, qsb, b0, delU1b, delU2b]]
    parameters.xRef = ub;
    % Linear inequality constraints: b02 < b01
    A = [1, 0, 0, 0, -1, 1, 0];
    b = 0;
elseif parameters.processType == "AdsorbentPVSAb0" % Reverse engineer adsorbent assuming equal deltaU for each gas using PVSA
    lb = [log10(0.02e5), 100,  40,  0.5, 1e-6, 1e-6, 10e3, log10(1e5)];   % Lower bounds [p_I, t_ads, t_blo, qsb, b0, delU1b, delU2b, p_H]
    ub = [log10(10e5),    3000, 300, 8,   1e-3, 1e-3, 45e3, log10(10e5)];  % Upper bounds [p_I, t_ads, t_blo, qsb, b0, delU1b, delU2b, p_H]]
    parameters.xRef = ub;
    % Linear inequality constraints: b02 < b01 & P_blo < P_ads
    A = [1, 0, 0, 0, -1, 1, 0,-1];
    b = 0;
    % elseif parameters.processType == "Resin" % Reverse engineer adsorbent assuming equal deltaU for each gas using PVSA
    %     lb = [log10(0.051e5), 1000,  20,  0.5, 400, log10(0.05e5), 0.5];   % Lower bounds [p_I, t_ads, t_blo, qsb, t_evac, p_L, v_in]
    %     ub = [log10(0.975e5),    30000, 400, 2.6, 20000, log10(0.5e5), 1];  % Upper bounds [p_I, t_ads, t_blo, qsb, t_evac, p_L, v_in]]
    %     parameters.xRef = ub;
    %     % Linear inequality constraints: b02 < b01 & P_blo < P_ads
    %     A = [-1, 0, 0, 0, 0, 1,0];
    %     b = 0;
elseif parameters.processType == "Resin" % Reverse engineer adsorbent assuming equal deltaU for each gas using PVSA
    lb = [log10(0.021e5),    1000,   20,   500,   log10(0.05e5), 0.005];   % Lower bounds [p_I, t_ads, t_blo, qsb, t_evac, p_L, v_in]
    ub = [log10(0.975e5),    30000, 1000, 45000, log10(0.2e5), 1];     % Upper bounds [p_I, t_ads, t_blo, qsb, t_evac, p_L, v_in]]
    parameters.xRef = ub;
    % Linear inequality constraints: b02 < b01 & P_blo < P_ads
    A = [-1, 0, 0, 0, 1,0];
    b = 0;
    % elseif parameters.processType == "Resin" % Reverse engineer adsorbent assuming equal deltaU for each gas using PVSA
    %     lb = [log10(0.051e5), 1000,  20,  0.5, 400, log10(0.05e5)];   % Lower bounds [p_I, t_ads, t_blo, qsb, t_evac, p_]
    %     ub = [log10(0.975e5),    30000, 400, 2.5, 20000, log10(0.5e5)];  % Upper bounds [p_I, t_ads, t_blo, qsb, t_evac, p_L]]
    %     parameters.xRef = ub;
    %     % Linear inequality constraints: b02 < b01 & P_blo < P_ads
    %     A = [-1, 0, 0, 0, 0, 1];
    %     b = 0;
end

parameters.fileName = convertStringsToChars(strcat(parameters.adsorbentName,"_",parameters.processType,"_",parameters.pressType,"_",parameters.OptType,"_",parameters.modelType,"_",datestr(now,'ddmmyyhhMM')));

if parameters.processType == "Resin"
    if ~parameters.fixResins
        qsbvals = [0.4:0.2:3.4];
        for kk = 1:length(qsbvals)
            parameters.qsb_1 = qsbvals(kk);
            %% Solve optimisation problem using genetic algorithm
            % Genetic algorithm settings
            nVars = length(lb);
            ngens = 50;     % Maximum number of generations
            pop_size = 120; % Number of members in population of each generation [-]

            rng default
            X0 = lhsdesign(pop_size,nVars); % Latin hypercube sampling to generate initial population matrix
            initPop = lb+X0.*(ub-lb); % Initial population matrix

            parameters.outputType = "opt"; % Output type needs to be "opt"
            options = optimoptions(@gamultiobj, 'Display', 'iter', 'Generations', ngens, 'PopulationSize', pop_size,'UseParallel',true,'InitialPopulationMatrix',initPop,'CrossoverFraction',0.85,'ParetoFraction',0.35,'PlotFcn','gaplotpareto','MigrationInterval',5); % ,'PlotFcn','gaplotpareto' GA options


            [x, fval] = gamultiobj(@(theta) kBAAM_Outputs_nonIsothermal_dP(parameters,theta), nVars, A, b, [], [], lb./parameters.xRef, ub./parameters.xRef, options); % Optimize objectives using GA

            save(['matFiles/',parameters.fileName,'.mat'])
        end
    else
        resinVals = [0.454357067	4.06469E-17	0.0399;
            0.51399776	6.11738E-17	0.034815;
            0.984989599	1.58041E-16	0.019913;
            1.336757211	1.65654E-16	0.01002;
            0.781131901	8.65318E-17	0.0247;
            1.174854923	1.15484E-16	0.0199;
            2.37	2.16419E-16	0.0007];

        for kk = 1:length(resinVals(:,1))
            parameters.qsb_1 = resinVals(kk,1);
            parameters.bo_1 = resinVals(kk,2);
            parameters.LDF = resinVals(kk,3);
            %% Solve optimisation problem using genetic algorithm
            % Genetic algorithm settings
            nVars = length(lb);
            ngens = 20;     % Maximum number of generations
            pop_size = 300; % Number of members in population of each generation [-]

            rng default
            X0 = lhsdesign(pop_size,nVars); % Latin hypercube sampling to generate initial population matrix
            initPop = lb+X0.*(ub-lb); % Initial population matrix

            parameters.outputType = "opt"; % Output type needs to be "opt"
            options = optimoptions(@gamultiobj, 'Display', 'iter', 'Generations', ngens, 'PopulationSize', pop_size,'UseParallel',true,'InitialPopulationMatrix',initPop,'CrossoverFraction',0.85,'ParetoFraction',0.35,'PlotFcn','gaplotpareto','MigrationInterval',5); % ,'PlotFcn','gaplotpareto' GA options


            [x, fval] = gamultiobj(@(theta) kBAAM_Outputs_nonIsothermal_dP(parameters,theta), nVars, A, b, [], [], lb./parameters.xRef, ub./parameters.xRef, options); % Optimize objectives using GA

            save(['matFiles/',parameters.fileName,'.mat'])
        end
    end
else

    %% Solve optimisation problem using genetic algorithm
    % Genetic algorithm settings
    nVars = length(lb);
    ngens = 70;     % Maximum number of generations
    pop_size = 120; % Number of members in population of each generation [-]

    rng default
    X0 = lhsdesign(pop_size,nVars); % Latin hypercube sampling to generate initial population matrix
    initPop = lb+X0.*(ub-lb); % Initial population matrix

    parameters.outputType = "opt"; % Output type needs to be "opt"
    % options = optimoptions(@gamultiobj, 'Display', 'iter', 'Generations', ngens, 'PopulationSize', pop_size,'UseParallel',true,'InitialPopulationMatrix',initPop,'CrossoverFraction',0.85,'ParetoFraction',0.35,'PlotFcn','gaplotpareto','MigrationInterval',5); % ,'PlotFcn','gaplotpareto' GA options
    options = optimoptions(@gamultiobj, 'Display', 'iter', 'Generations', ngens, 'PopulationSize', pop_size,'UseParallel',true,'CrossoverFraction',0.85,'ParetoFraction',0.35,'PlotFcn','gaplotpareto','MigrationInterval',5); % ,'PlotFcn','gaplotpareto' GA options


    [x, fval] = gamultiobj(@(theta) kBAAM_Outputs_nonIsothermal_dP(parameters,theta), nVars, A, b, [], [], lb./parameters.xRef, ub./parameters.xRef, options); % Optimize objectives using GA

    save(['matFiles/',parameters.fileName,'.mat'])
end
end