% clc; clear all;
% close all;
addpath(genpath(pwd))
load('13X_AKR_16.mat')
load('13X_AW_22_2.mat')
% load('CSAC_AKR_16.mat');
% parameters.lambda = 0.5;
maxTarget = "Productivity";

%
% fileNames =    {'13X_AKR_VSA_kBAAM_2209250845.txt',...
%     'CSAC_AKR_VSA_kBAAM_2209250850.txt'};

% fileNames =    {'CSAC_AKR_VSA_kBAAM_NI_2309251301.txt'};
% fileNames =    {'13X_AKR_VSA_kBAAM_2309251228.txt'};
fileNames =    {'Z13X_PVSA_Const_nonisothermal_0310251317.txt'};

for jj = 1:length(fileNames)

    fileName = fileNames{jj};
    data = load(fileName);

    data = data(find(data(:,8)>95),:);
    data = data(find(data(:,9)>90),:);

    try
        if maxTarget == "Purity"
            Obj = data(:,8);
        elseif maxTarget == "Recovery"
            Obj = data(:,9);
        elseif maxTarget == "Productivity"
            Obj = data(:,10);
        elseif maxTarget == "both"
            Obj = sqrt(data(:,8).^2+data(:,9).^2);
        else
        end
    catch
        Obj = data(:,8);
    end

    thetaVals = data(find(Obj==max(Obj)),:);
    thetaVals = thetaVals(1,:);
    theta1 = [thetaVals(end) thetaVals(2) thetaVals(5) thetaVals(6) thetaVals(7) thetaVals(1)]; % vector of decision variables, [F_in, P_I, t_ads, t_blo, t_evac]

    parameters.rp = 5e-4;
    parameters.p_L = thetaVals(3);
    parameters.V_column = thetaVals(12);
    parameters.outputType = "opt";
    parameters.processType = "PVSA";
    parameters.pressType = "FP";
    parameters.outputType = "plot";
    parameters.OptType = "unc";
    % parameters.y1_in = 0.01;
    % parameters.OptType = "Const";
    parameters.fileName = 'testsads';
    parameters.modelType = "isothermal";
    parameters.t_press = 20;
    % parameters.modelType = "isothermal";
    %% Run to CSS and output KPIs corresponding to theta
    % KPIs = kBAAM_Outputs_isothermal(parameters,theta1) % KPIS = [Recovery, Purity];
    KPIs = kBAAM_Outputs_nonIsothermal(parameters,theta1) % KPIS = [Recovery, Purity];

end