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
close all
clear

saveFlag = 0;
ProductivityFlag = 0;
colorVals = ["b","r","m","none"];
EdgecolorVals = ["b","r","m","k"];
markerVals = ["square", "square","o","o"];


matNames = ["Zeolite 13X", "Activated Carbon"];

fileNames =    {'Z13X_PVSA_Const_nonisothermal_0810251603',...
    'Z13X_VSA_Const_nonisothermal_0810251609'};

figure
tiledlayout(1,1,"Padding","tight","TileSpacing","compact")
set(gcf,'units','inch','position',[5,5,6,4],'color','w')


for mm = 1:length(fileNames)
    fileName = fileNames{mm};
    data = readmatrix([fileName, '.txt']);

    input = [];
    conInd = [];
    ab = [];
    for kk = 1:length(data(:,1))
        if data(kk,8) > 95 && data(kk,9) > 90
            input = [input; 1./data(kk,10),data(kk,11)];
            ab = [ab;data(kk,:)];
            conInd = [conInd;kk];
        else
            conInd = [conInd;1];
        end
    end
    hold on

    scatter(ab(:,10),ab(:,11).*2.77778e-7.*1e3,1,'MarkerEdgeAlpha',0.05,'MarkerEdgeColor',colorVals(mm),'MarkerFaceColor','none','MarkerFaceAlpha',0.1,'HandleVisibility','off');
    [flag, ParetoPoints]=find_pareto_frontier(input);
    hold on
    outPareto = [1./ParetoPoints(:,1),ParetoPoints(:,2)];
    scatter(outPareto(:,1),outPareto(:,2).*2.77778e-7.*1e3,50,markerVals(mm),'MarkerEdgeColor','none','LineWidth',0.8,'MarkerFaceColor',colorVals(mm),'MarkerFaceAlpha',1);
    box on; grid on;
    set(gca,'YScale','linear','XScale','linear','FontSize',12,'LineWidth',1.5,'GridLineWidth',1)
    ylabel('E_{T} [kWh/tonne]');
    xlabel('CO_{2} Productivity [mol/m3/s]')
end

try
    load("AdamConst002.mat")
catch
    AdamData = [  0.713159491000000	433.823529400000;
        1.15762078000000	490.196078400000;
        1.50975338600000	568.627451000000;
        2.15933899300000	573.529411800000;
        2.97417626800000	593.137254900000;
        3.47003234300000	642.156862700000;
        4.21290681200000	671.568627500000;
        4.30609460300000	691.176470600000;
        4.76076409900000	742.647058800000;
        4.92657166000000	784.313725500000;
        4.93784111600000	830.882352900000;
        5.20664038800000	867.647058800000;
        5.34066100700000	867.647058800000;
        5.46447341800000	872.549019600000;
        5.78421265400000	879.901960800000;
        5.92914897900000	909.313725500000;
        6.00197089100000	941.176470600000;
        6.18824540100000	975.490196100000;
        6.30179907000000	982.843137300000;
        6.38457651100000	997.549019600000;
        6.47761269500000	1009.80392200000;
        6.51895087900000	1014.70588200000;
        6.60162724900000	1024.50980400000;
        6.68440469000000	1039.21568600000;
        6.79841318000000	1068.62745100000;
        6.71907216500000	1220.58823500000];
end

scatter(AdamData(:,1),AdamData(:,2),50,'o',Marker='diamond',MarkerEdgeColor='k',MarkerFaceColor='g',MarkerFaceAlpha=1);
legend('PVSA', 'VSA','Adam (TSEMO)')
