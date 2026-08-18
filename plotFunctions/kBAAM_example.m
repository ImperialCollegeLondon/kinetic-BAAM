%% Run to CSS and output KPIs corresponding to theta
createParameters

load("13XH_T.mat")
parameters.outputType = "plot";
parameters.processType = "PVSA";
parameters.OptType = "Unc";
parameters.Lbyr = 7;
parameters.equilibrium = 0;
parameters.cCSTR = 0 ;
parameters.testEvac = 0;
parameters.normPlot = 0;
% parameters.modelType = "isothermal";
parameters.forwardEvac = 0;
parameters.pressureDrop = 1;
parameters.plotVideo = 0;
parameters.equilibrium = 0;
parameters.plot0D = 1;
parameters.pressType = "FP";
theta = [0.3.*0.37 0.4e5 110 40 100 0.02e5 3e5];% Upper bounds [v_in, p_I, t_ads, t_blo, t_evac]

tic
KPIs2 = kBAAM_Outputs_nonIsothermal_dP(parameters,theta) % KPIS = [-Recovery, -Purity];
kbaamDP = toc
