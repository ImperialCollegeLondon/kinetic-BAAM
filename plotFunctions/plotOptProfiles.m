%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Imperial College London, United Kingdom
% Multiphase Systems Laboratory
% Year:     2025
% MATLAB:   R2024a
% Authors:  Hassan Azzan (HA)
%
% Purpose:
% Plot cyclic steady state profiles for the optimal design based on
% specific targets
%
% Last modified:
% - 2025-09-21, HA: Initial creation
%
% Input arguments:
%
% Output arguments:
%
% Dependencies:
%   - kBAAM_Outputs_nonIsothermal.m
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clc; clear all;
close all;
addpath(genpath(pwd))
maxTarget = "bothcon"; % Purity, Recovery, Productivity, both (Pu and Rec)
fileNames =    {'Z13X_PVSA_Const_nonisothermal_0810251603'};

for jj = 1:length(fileNames)

    fileName = fileNames{jj};
    data = load([fileName, '.txt']);
    load([fileName, '.mat'])
    data = data(find(data(:,8)>95),:);
    data = data(find(data(:,9)>90),:);
    load("Z13X_AW_2022_4.mat")

    try
        if maxTarget == "Purity"
            Obj = data(:,8);
        elseif maxTarget == "Recovery"
            Obj = data(:,9);
        elseif maxTarget == "Productivity"
            Obj = data(:,10);
        elseif maxTarget == "Energy"
            Obj = -data(:,11);
        elseif maxTarget == "both"
            Obj = sqrt(data(:,8).^2+data(:,9).^2);
        elseif maxTarget == "bothcon"
            Obj = sqrt((data(:,10)./max(data(:,10))).^2+((1./data(:,11)./max(1./data(:,11)))).^2);
        else
        end
    catch
        Obj = data(:,8);
    end

    thetaVals = data(find(Obj==max(Obj)),:);
    thetaVals = thetaVals(1,:);
    theta = [thetaVals(end) thetaVals(2) thetaVals(5) thetaVals(6) thetaVals(7) thetaVals(1)]; % vector of decision variables, [F_in, P_I, t_ads, t_blo, t_evac, p_H]
    parameters.outputType = "plot";
    parameters.OptType = "Unc";

    %% Run to CSS and output KPIs corresponding to theta
    tic
    KPIs = kBAAM_Outputs_nonIsothermal(parameters,theta) % KPIS = [-Recovery, -Purity];
    kbaam = toc

    model1D = 14.079827;
    
    x = ["0-D" "1-D"];
    y = [kbaam model1D];
    figure
    % bar(y)

    x = ["0D - this work" "1D"];
y = [kbaam model1D];
bar(x,y)
ylabel('Simulation time [s]')
set(gca,'FontSize',30)
box on;grid off;set(gca,'YScale','linear','XScale','linear','LineWidth',2)
set(gcf,'Color','white')

    % parameters.h_out = 1e5;
    % tic
    % KPIs = kBAAM_Outputs_nonIsothermal(parameters) % KPIS = [-Recovery, -Purity];
    % toc
end