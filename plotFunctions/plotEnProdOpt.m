addpath(genpath(pwd))
close all
clear

saveFlag = 0;
ProductivityFlag = 0;
colorVals = ["b","r","m","none"];
EdgecolorVals = ["b","r","m","k"];
markerVals = ["square", "square","o","o"];

% colorVals = ["k","none","none","none"];
% EdgecolorVals = ["k","k","b","k"];

matNames = ["Zeolite 13X", "Activated Carbon"];


fileNames =    {'Z13X_PVSA_Const_isothermal_0510251103.txt',...
    'Z13X_PVSA_Const_nonisothermal_0510251049.txt'};
% fileNames =    {'Hypo_AdsorbentVSA_Const_isothermal_0110251600.txt',...
%     'Hypo_AdsorbentVSA_Const_nonisothermal_0110251550.txt'};


% fileNames =    {'Z13X_PVSA_Const_isothermal_0210251241.txt'};
% fileNames =    {'Z13X_PVSA_Const_isothermal_0210251241.txt'};

load("litOpt.mat")
figure
tiledlayout(1,1,"Padding","tight","TileSpacing","compact")
set(gcf,'units','inch','position',[5,5,6,4],'color','w')


for mm = 1:length(fileNames)
    Targetflag = 0;
    data = readmatrix(fileNames{mm});
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

    % scatter(data(:,10),data(:,11).*2.77778e-7.*1e3,9,'MarkerEdgeAlpha',0.1,'MarkerEdgeColor',colorVals(mm),'MarkerFaceColor','none','MarkerFaceAlpha',0.1,'HandleVisibility','off');
    % scatter(ab(:,10),ab(:,11).*2.77778e-7.*1e3,1,'MarkerEdgeAlpha',0.05,'MarkerEdgeColor',colorVals(mm),'MarkerFaceColor','none','MarkerFaceAlpha',0.1,'HandleVisibility','off');
    [flag, ParetoPoints]=find_pareto_frontier(input);
    hold on
    outPareto = [1./ParetoPoints(:,1),ParetoPoints(:,2)];
    scatter(outPareto(:,1),outPareto(:,2).*2.77778e-7.*1e3,50,markerVals(mm),'MarkerEdgeColor','none','LineWidth',0.8,'MarkerFaceColor',colorVals(mm),'MarkerFaceAlpha',1);
    box on; grid on;
    set(gca,'YScale','linear','XScale','linear','FontSize',12,'LineWidth',1.5,'GridLineWidth',1)
    ylabel('E_{T} [kWh/tonne]');
    xlabel('CO_{2} Productivity [mol/m3/s]')
    % ylim([0 1200])
    % xlim([0.2 8])
end
load("AdamConst002.mat")
scatter(AdamData(:,1),AdamData(:,2),50,'o',Marker='diamond',MarkerEdgeColor='k',MarkerFaceColor='g',MarkerFaceAlpha=0.2);


legend('Isothermal', 'Nonisothermal','Adam (TSEMO)')

% figure
% scatter(data(:,10),data(:,5)+data(:,6)+data(:,7)+30)
% figure
% scatter(data(:,10),data(:,4).*1e5./data(:,1))