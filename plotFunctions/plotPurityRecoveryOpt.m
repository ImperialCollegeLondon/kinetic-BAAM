%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Imperial College London, United Kingdom
% Multiphase Systems Laboratory
% Year:     2025
% MATLAB:   R2024a
% Authors:  Hassan Azzan (HA)
%
% Purpose:
% Plot results from unconstrained optimization Recovery vs Purity
%
% Last modified:
% - 2025-09-21, HA: Initial creation
%
% Input arguments:
%
% Output arguments:
%
% Dependencies:
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
addpath(genpath(pwd))
close all
clear

saveFlag = 0;
ProductivityFlag = 0;
colorVals = ["b","k","none","none"];
EdgecolorVals = ["b","k","b","k"];
markerVals = ["square", "square","o","o"];

matNames = ["Zeolite 13X", "Activated Carbon"];

fileNames =    {'13X_AKR_VSA_kBAAM_NI_2309251238.txt',...
    'CSAC_AKR_VSA_kBAAM_NI_2309251301.txt'};


figure
tiledlayout(1,1,"Padding","tight","TileSpacing","compact")
set(gcf,'units','inch','position',[5,5,6,4],'color','w')
for mm = 1:length(fileNames)
    Targetflag = 0;
    data = readmatrix(fileNames{mm});

    input = [];
    data(:,10) = data(:,10) .*0.044.*3600.*24;
    for kk = 1:length(data(:,1))
        data(kk,8) = max(0,data(kk,8));
        data(kk,9) = max(0,data(kk,9));

        if data(kk,8) > 100 || data(kk,9) > 100
            data(kk,8) = 0;
            data(kk,9) = 0;
        end
        input = [input; 1./data(kk,8),1./data(kk,9)];
    end
    hold on
    p = patch('vertices', [95, 90; 100, 90; 100, 100; 95 100], ...
        'faces', [1, 2, 3, 4], ...
        'FaceColor', 'b', ...
        'FaceAlpha', 0.05,'HandleVisibility','off');

    [flag, ParetoPoints]=find_pareto_frontier(input);
    hold on
    outPareto = [1./ParetoPoints(:,1),1./ParetoPoints(:,2)];
    p1 = scatter(data(flag,8),data(flag,9),50,markerVals(mm),'MarkerEdgeColor',EdgecolorVals(mm),'LineWidth',0.8,'MarkerFaceColor',colorVals(mm),'MarkerFaceAlpha',1);
    box on; grid on;
    set(gca,'YScale','linear','XScale','linear','FontSize',12,'LineWidth',1.5,'GridLineWidth',1)
    xlabel('CO_{2} Purity [%]');
    ylabel('CO_{2} Recovery [%]')
    ylim([60 100]);
    xlim([80 100]);
end
try
    load("litOpt.mat")
    scatter(litAC(:,1),litAC(:,2),50,'k',Marker='diamond',MarkerEdgeColor='g',MarkerFaceColor='k',MarkerFaceAlpha=0.5);
    scatter(lit13X(:,1),lit13X(:,2),50,'b',Marker='diamond',MarkerEdgeColor='g',MarkerFaceColor='b',MarkerFaceAlpha=0.5);
    scatter(LitAdam(:,1),LitAdam(:,2),50,'k',Marker='diamond',MarkerEdgeColor='m',MarkerFaceColor='m',MarkerFaceAlpha=1);
    legend('13X', 'AC', 'Balashankar et al. 2019 (1D)','Balashankar et al.  2019 (1D)', 'Ward & Pini 2022 (1D)','Location','southwest')
catch
    disp('fileNotFound')
end
