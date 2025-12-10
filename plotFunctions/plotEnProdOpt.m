
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Imperial College London, United Kingdom
% Multiphase Systems Laboratory
% Year:     2025
% MATLAB:   R2024a
% Authors:  Hassan Azzan (HA)
%
% Purpose:
% Plot results from constrained optimization Energy vs Productivity
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
% close all
clear

saveFlag = 0;
ProductivityFlag = 0;
colorVals = ["b","r","m","none"];
EdgecolorVals = ["b","r","m","k"];
markerVals = ["square", "square","o","o"];

load('Z13X_AW_2022.mat')
fsz = 15;

matNames = ["Zeolite 13X", "Activated Carbon"];

fileNames =    {'Z13X_PVSA_Const_nonisothermal_0910250741',...
    'Z13X_PVSA_Const_nonisothermal_1010250858'};
fileNames =    {['Resin_Resin_Const_nonisothermal_0512251300']};
fileNames =    {['Resin_Resin_Const_nonisothermal_0512251656']};

figure
tiledlayout(1,1,"Padding","tight","TileSpacing","compact")
set(gcf,'units','inch','position',[5,5,7,5],'color','w')
% title('Resin Reverse Engineering')

dataTemp = readmatrix([fileNames{1}, '.txt']);
qsbvals = unique(dataTemp(:,17));
input = [];
conInd = [];
abTemp = [];


ColorsForPlotAll = brewermap(length(qsbvals),'RdYlGn');
ColorsForPlot = ColorsForPlotAll;
C = [175,1,1;
    215,41,0;
    255,82,0;
    255,169,3;
    254,255,7;
    127,253,3;
    0,250,0;
    0,156,7;
    0,63,13]./255;
ColorsForPlot = C;
ColorsForPlot = createcolormap(length(qsbvals),C(1,:),C(2,:),C(3,:),C(4,:),C(5,:),C(6,:),C(7,:),C(8,:),C(9,:));

for mm = 2:length(qsbvals)
    data = dataTemp(dataTemp(:,17)==qsbvals(mm),:);
    input = [];
    conInd = [];
    ab = [];
    for kk = 1:length(data(:,1))
        if data(kk,8) > 95 && data(kk,9) > 0
            input = [input; 1./data(kk,10),data(kk,11)];
            ab = [ab;data(kk,:)];
            conInd = [conInd;kk];
        else
            conInd = [conInd;1];
        end
    end
    % data = sortrows(ab,17);
    [flag, ParetoPoints]=find_pareto_frontier(input);
    outPareto = [1./ParetoPoints(:,1),ParetoPoints(:,2)];

    hold on
    kvals = 0.132.*exp(-2.076.*qsbvals);
    scatter(ab(flag,10).*0.044.*24.*3600,ab(flag,11)*1e-6,50,repmat(qsbvals(mm),1,length(ab(flag,10))),'filled',MarkerFaceAlpha=1);
    % scatter(ab(flag,10).*0.044.*24.*3600,ab(flag,11)*1e-6,50,repmat(kvals(mm),1,length(ab(flag,10))),'filled',MarkerFaceAlpha=1);

    clim([min(qsbvals) max(qsbvals)]);
    colormap(ColorsForPlot)
    hcb=colorbar('ver'); % colorbar handle
    hcb.FontSize = fsz;
    hcb.Title.String = "q_{sat}";
    hcb.Title.String = "k_{LDF}";
    hcb.Title.FontSize = fsz;
    ylim([0 12])
    xlim([0 120])
    hcb.Limits = [min(qsbvals) max(qsbvals)];
    % hcb.Limits = [[min(kvals) max(kvals)]];
    box on; grid off;
    yline(1.66,'LineWidth',1,'Color','k','LineStyle','--')
    set(gca,'YScale','linear','XScale','linear','FontSize',fsz,'LineWidth',2,'GridLineWidth',1,'Color',[0.9 0.9 0.9])
    ylabel('Specific Energy Consumption [MJ/kg_{CO_{2}}]');
    xlabel('CO_{2} Productivity [kg_{CO_{2}}/m_{bed}^{3}/day]')
end
set(gcf,'units','inch','position',[5,5,10,6],'color','w')

patch([17.1 17.1 40.7 40.7],[3.3 5.1 5.1 3.3],'cyan','FaceAlpha',0.01,'EdgeColor','b','LineWidth',1.5);