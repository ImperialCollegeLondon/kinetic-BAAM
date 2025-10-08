%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Imperial College London, United Kingdom
% Multiphase Systems Laboratory
% Year:     2025
% MATLAB:   R2024a
% Authors:  Hassan Azzan (HA)
%
% Purpose:
% Plot optimal isotherms from reverse engineering based on
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
%   - DSL.m
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clc; clear all;
close all;
addpath(genpath(pwd))
maxTarget = "Productivity"; % Purity, Recovery, Productivity, both (Pu and Rec)addpath(genpath(pwd))

Pvals = linspace(0,1e5,1000);
T = 298;

fileNames =    {'Hypo_AdsorbentPVSA_Unc_nonisothermal_0310251213'};

for kk = 1:length(fileNames)
    fileName = fileNames{kk};
    data = readmatrix([fileName, '.txt']);
    load([fileName, '.mat']);
    figure
    try
        if maxTarget == "Purity"
            Obj = data(:,8);
        elseif maxTarget == "Recovery"
            Obj = data(:,9);
        elseif maxTarget == "Productivity"
            Obj = data(:,10);
        elseif maxTarget == "Energy"
            Obj = -data(:,11);
        elseif maxTarget == "EnPu"
            Obj = sqrt((data(:,8)./100).^2+((1./data(:,11)./max(1./data(:,11)))).^2);
        elseif maxTarget == "EnPuRe"
            Obj = sqrt((data(:,8)./100).^2+((1./data(:,11)./max(1./data(:,11)))).^2+(data(:,9)./100).^2);
        elseif maxTarget == "both"
            Obj = sqrt(data(:,8).^2+data(:,9).^2);
        else
        end
    catch
        Obj = data(:,8);computeStatSTALoading2
    end

    ObjTemp = Obj(1);

    for jj = 1:length(data(:,1))
        if data(jj,8) > 95 && data(jj,9) > 90
            thetaVals = data(jj,:);
            thetaVals = thetaVals(1,:);
            theta = thetaVals(end-6:end); % vector of decision variables, [F_in, P_I, t_ads, t_blo, t_evac]

            parameters.qsb_1 = theta(4);
            parameters.qsb_2 = theta(4);
            parameters.qsd_1 = 0;
            parameters.qsd_2 = 0;
            parameters.bo_1 = theta(5);
            parameters.bo_2 = theta(5);
            parameters.do_1 = 0;
            parameters.do_2 = 0;
            parameters.delUb_1 = -theta(6);
            parameters.delUb_2 = -theta(7);
            parameters.delUd_1 = 0;
            parameters.delUd_2 = 0;

            [q1, x] = DSL(Pvals, 1, T, parameters.qsb_1, parameters.qsd_1, parameters.qsb_2, parameters.qsd_2, parameters.bo_1, parameters.do_1, parameters.bo_2, parameters.do_2, parameters.delUb_1, parameters.delUd_1, parameters.delUb_2, parameters.delUd_2);
            [x, q2] = DSL(Pvals, 0, T, parameters.qsb_1, parameters.qsd_1, parameters.qsb_2, parameters.qsd_2, parameters.bo_1, parameters.do_1, parameters.bo_2, parameters.do_2, parameters.delUb_1, parameters.delUd_1, parameters.delUb_2, parameters.delUd_2);


            subplot(1,2,1)
            hold on
            plot(Pvals./1e5,q1, 'Color',[0 0 0 0.05],'LineWidth', 0.1,'LineStyle','-'); xlabel('P [bar]'); ylabel('Adsorbed amount of CO_2 [mol/kg]'); hold on;%check unit
            % xline(parameters.y1_in.*max(Pvals),  'k--', 'LineWidth', 1, 'HandleVisibility','off');
            box on; grid off; set(gca,"LineWidth",2,"FontName","CMU Serif","FontSize",15)
            subplot(1,2,2)
            hold on
            plot(Pvals./1e5,q2, 'Color',[0 0 0 0.05],'LineWidth', 0.1,'LineStyle','-'); xlabel('P [bar]'); ylabel('Adsorbed amount of N_2 [mol/kg]'); hold on;%check unit
            % xline(1-parameters.y1_in.*max(Pvals),  'k--', 'LineWidth', 1, 'HandleVisibility','off');
            box on; grid off; set(gca,"LineWidth",2,"FontName","CMU Serif","FontSize",15)

            if Obj(jj) > ObjTemp
                ObjTemp = Obj(jj);
                thetaVals = data(jj,:);
            end

        end
    end


    theta = thetaVals(end-6:end); % vector of decision variables, [F_in, P_I, t_ads, t_blo, t_evac]

    parameters.qsb_1 = theta(4);
    parameters.qsb_2 = theta(4);
    parameters.qsd_1 = 0;
    parameters.qsd_2 = 0;
    parameters.bo_1 = theta(5);
    parameters.bo_2 = theta(5);
    parameters.do_1 = 0;
    parameters.do_2 = 0;
    parameters.delUb_1 = -theta(6);
    parameters.delUb_2 = -theta(7);
    parameters.delUd_1 = 0;
    parameters.delUd_2 = 0;

    [q1, x] = DSL(Pvals, 1, T, parameters.qsb_1, parameters.qsd_1, parameters.qsb_2, parameters.qsd_2, parameters.bo_1, parameters.do_1, parameters.bo_2, parameters.do_2, parameters.delUb_1, parameters.delUd_1, parameters.delUb_2, parameters.delUd_2);
    [x, q2] = DSL(Pvals, 0, T, parameters.qsb_1, parameters.qsd_1, parameters.qsb_2, parameters.qsd_2, parameters.bo_1, parameters.do_1, parameters.bo_2, parameters.do_2, parameters.delUb_1, parameters.delUd_1, parameters.delUb_2, parameters.delUd_2);

    subplot(1,2,1)
    hold on
    plot(Pvals./1e5,q1, 'Color','r','LineWidth', 3); xlabel('P [bar]'); ylabel('Adsorbed amount of CO_2 [mol/kg]'); hold on;%check unit
    % xline(parameters.y1_in.max(Pvals),  'k--', 'LineWidth', 1, 'HandleVisibility','off');
    box on; grid off; set(gca,"LineWidth",2,"FontName","CMU Serif","FontSize",15)
    subplot(1,2,2)
    hold on
    plot(Pvals./1e5,q2, 'Color','b','LineWidth', 3); xlabel('P [bar]'); ylabel('Adsorbed amount of N_2 [mol/kg]'); hold on;%check unit
    % xline(parameters.y1_in.max(Pvals),  'k--', 'LineWidth', 1, 'HandleVisibility','off');
    box on; grid off; set(gca,"LineWidth",2,"FontName","CMU Serif","FontSize",15)
end