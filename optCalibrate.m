%% Run to CSS and output KPIs corresponding to theta
% run createParameters.m first to generate the parameter struct
load('Z13X_AW_2022.mat')
load('funVals1D4.mat')

% load('CSAC_RH_2014.mat')
% load('funVals1D_AC.mat')

parameters.adsorbentName = "13XCalibrate_fin";
% parameters.adsorbentName = "CSACCalibrate_fin";

parameters.outputType = "Opt";
parameters.processType = "VSA";
parameters.modelType = "nonisothermal";
parameters.OptType = "sampling";
parameters.fileName = convertStringsToChars(strcat(parameters.adsorbentName,"_",parameters.processType,"_",parameters.OptType,"_",parameters.modelType,"_",datestr(now,'ddmmyyhhMM')));

lb = [1,     1,     0.01,    0.01];   % Lower bounds [h_in, h_out, lambdab, lambdae, Lbyr]
ub = [30,    20,     5,  5];  % Upper bounds [h_in, h_out, lambdab, lambdae, Lbyr]
% Linear inequality constraints: none
A = [0, 0, 0, 0];
b = 0;

nVars = length(lb);
ngens = 50;     % Maximum number of generations
pop_size = 30; % Number of members in population of each generation [-]

rng default
X0 = lhsdesign(pop_size,nVars); % Latin hypercube sampling to generate initial population matrix

initPop = lb+X0.*(ub-lb); % Initial population matrix
% 
% options = optimoptions(@gamultiobj, 'Display', 'iter', 'Generations', ngens, 'PopulationSize', pop_size,'UseParallel',true,'InitialPopulationMatrix',initPop,'CrossoverFraction',0.85,'ParetoFraction',0.35,'PlotFcn','gaplotpareto','MigrationInterval',5); % ,'PlotFcn','gaplotpareto' GA options
% 
% 
% [x, fval] = gamultiobj(@(theta) objFuncCalib(theta,parameters,funVals1D), nVars, A, b, [], [], lb, ub, options); % Optimize objectives using GA


options = optimoptions(@ga, 'Display', 'iter', 'Generations', ngens, 'PopulationSize', pop_size,'UseParallel',true,'InitialPopulationMatrix',initPop,'CrossoverFraction',0.85,'PlotFcn','gaplotbestf','MigrationInterval',5); % ,'PlotFcn','gaplotpareto' GA options


[x, fval] = ga(@(theta) objFuncCalib(theta,parameters,funVals1D), nVars, A, b, [], [], lb, ub, [], options) % Optimize objectives using GA


parameters.h_in = x(1);
parameters.h_out = x(2);
parameters.lambdab = x(3);
parameters.lambdae = x(4);
% parameters.Lbyr = x(5);

% save(['AdsorbentFiles/',convertStringsToChars(strcat(parameters.adsorbentName)),'.mat'],'parameters')
