addpath(genpath(pwd))
close all
clear

saveFlag = 0;
ProductivityFlag = 0;
colorVals = ["b","k","none","none"];
EdgecolorVals = ["b","k","b","k"];
markerVals = ["square", "square","o","o"];

colorVals = ["k","none","none","none"];
EdgecolorVals = ["k","k","b","k"];

% colorVals = ["b","none","none","none"];
% EdgecolorVals = ["b","b","b","k"];


matNames = ["Zeolite 13X", "Activated Carbon"];


fileNames =    {'13X_AKR_VSA_kBAAM_2209251156.txt',...
                'CSAC_AKR_VSA_kBAAM_2209251153.txt',...
                '13X_AKR_VSA_kBAAM_2209250845.txt',...
                'CSAC_AKR_VSA_kBAAM_2209250850.txt'};

fileNames =    {'13X_AKR_VSA_kBAAM_2209251156.txt',...
                '13X_AKR_VSA_kBAAM_2209250845.txt'};


fileNames =    {'13X_AKR_VSA_kBAAM_2309251228.txt',...
                'CSAC_AKR_VSA_kBAAM_2309251257.txt',...
                '13X_AKR_VSA_kBAAM_NI_2309251238.txt',...
                'CSAC_AKR_VSA_kBAAM_NI_2309251301.txt'};


fileNames =    {'13X_AKR_VSA_kBAAM_NI_2309251238.txt',...
                '13X_AKR_VSA_kBAAM_2209250845.txt'};


fileNames =    {'CSAC_AKR_VSA_kBAAM_NI_2309251301.txt',...
                'CSAC_AKR_VSA_kBAAM_2309251257.txt'};

% 
% fileNames =    {'13X_AKR_VSA_kBAAM_2309251228.txt',...
%                 '13X_AKR_VSA_kBAAM_NI_2309251238.txt'};

load("litOpt.mat")
figure
tiledlayout(1,1,"Padding","tight","TileSpacing","compact")
set(gcf,'units','inch','position',[5,5,6,4],'color','w')

for mm = 1:length(fileNames)
    Targetflag = 0;
    data = load(fileNames{mm});
    
    if ProductivityFlag
        colY = 10;
    else
        colY = 9;
    end
    input = [];
    data(:,10) = data(:,10) .*0.044.*3600.*24;
    for kk = 1:length(data(:,1))
        data(kk,8) = max(0,data(kk,8));
        data(kk,colY) = max(0,data(kk,colY));
        
        if data(kk,8) > 100 || data(kk,9) > 100
            data(kk,8) = 0;
            data(kk,9) = 0;
        end
        input = [input; 1./data(kk,8),1./data(kk,colY)];
    end
    % nexttile
    hold on
    if ~ProductivityFlag && mm == 1
    p = patch('vertices', [95, 90; 100, 90; 100, 100; 95 100], ...
        'faces', [1, 2, 3, 4], ...
        'FaceColor', 'b', ...
        'FaceAlpha', 0.1,'HandleVisibility','off');
    elseif mm == 1

    p = patch('vertices', [95, 0; 2e5, 0; 2e5, 2e5; 95 2e5], ...
        'faces', [1, 2, 3, 4], ...
        'FaceColor', 'b', ...
        'FaceAlpha', 0.1,'HandleVisibility','off');
    end
    % scatter(1*data(:,8),1*data(:,colY),7,'MarkerEdgeAlpha',1,'MarkerEdgeColor',colorVals(mm),'MarkerFaceColor',colorVals(mm),'MarkerFaceAlpha',0.1);
    [flag, ParetoPoints]=find_pareto_frontier(input);
    hold on
    outPareto = [1./ParetoPoints(:,1),1./ParetoPoints(:,2)];
    p1 = scatter(data(flag,8),data(flag,colY),50,markerVals(mm),'MarkerEdgeColor',EdgecolorVals(mm),'LineWidth',0.8,'MarkerFaceColor',colorVals(mm),'MarkerFaceAlpha',1);
    box on; grid on; 
    % axis square
    set(gca,'YScale','linear','XScale','linear','FontSize',12,'LineWidth',1.5,'GridLineWidth',1)

    xlabel('CO_{2} Purity [%]');
    
    if ~ProductivityFlag
    ylabel('CO_{2} Recovery [%]')
    ylim([60 100]);
    else
    ylabel('CO_{2} Productivity [kg/m3/day]')
    ylim([0 1600]);
    end


    xlim([80 100]);
    % title(matNames(mm))
    % adsname = fileNames{mm};
    % adsname = strsplit((adsname),'_');
    % FileNameX = ['OptPops/',adsname(2),'_OptPops_',adsname(5),'_','0pt',adsname(6),'bar','.xlsx'];
    % FileNameX = [FileNameX{:}];
    % OptimizationPop = table([1:length(data(:,1))]', data(:,1)./1e5, data(:,2)./1e5, data(:,3)./1e5, data(:,4), data(:,5), data(:,6), data(:,7), data(:,9), data(:,11),data(:,12),data(:,13), ...
    %     'VariableNames',["Simulation number [-]", "Adsorption pressure [bar]", "Blowdown pressure [bar]", "Evacuation pressure [bar]", "Feed temperature [K]", "Desorption temperature [K]",...
    %     "Adsorption time [s]", "Blowdown time [s]", "Evacuation and Heating time[s]","CO_2 purity [%]","CO_2 Recovery [%]","CO_2 productivity [mol/m3/s]"]);
    % writetable(OptimizationPop,FileNameX,'Sheet','Optimization Points','WriteVariableNames',true);
end

scatter(litAC(:,1),litAC(:,2),50,'k',Marker='diamond',MarkerEdgeColor='g',MarkerFaceColor='k',MarkerFaceAlpha=0.5);
% scatter(lit13X(:,1),lit13X(:,2),50,'b',Marker='diamond',MarkerEdgeColor='g',MarkerFaceColor='b',MarkerFaceAlpha=0.5); scatter(LitAdam(:,1),LitAdam(:,2),50,'k',Marker='diamond',MarkerEdgeColor='m',MarkerFaceColor='m',MarkerFaceAlpha=1);

legend('Nonisothermal', 'Isothermal', 'Balashankar 2019', 'Adam 2022')
%
% if saveFlag
%     FileNameF = ['Figures/',adsname{2},'_PurRec'];
%
%     saveas(gcf,[FileNameF,'.pdf'])
% end
