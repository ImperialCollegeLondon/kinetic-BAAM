%%
addpath(genpath('/Users/ha3215/Box Sync/Hassan Azzan Supervision Folder/Process Models/DBAAM_FP/'))
addpath(genpath(pwd))
load('Z13X_AW_2022.mat')
% load('CSAC_RH_2014.mat')
load('funVals1D_13X_new.mat')
load('funvals1D_13X_all.mat')
funVals1D = sortrows(funVals1D,10,"descend");
funVals1D(funVals1D(:,9)==0,:) = [];
% funVals1D(funVals1D(:,8)<10,:) = [];
% load('funVals1D_AC.mat')
parameters.adsorbentName = "13XCalibrate_fin";
parameters.outputType = "Opt";
parameters.processType = "VSA";
parameters.modelType = "nonisothermal";
parameters.OptType = "sampling";

addpath(genpath('/Users/ha3215/Box Sync/Hassan Azzan Supervision Folder/Process Models/DBAAM_FP/'))

plot0D = 1;
plot0Dp = 0;


% % % funVals1D(:,9) = funVals1D(:,9).*funVals1D(:,10).*(funVals1D(:,3)+funVals1D(:,4)+funVals1D(:,5)+20).*parameters.V_column.*3.600.*0.044;
% % Fvals0D = [];
% % Fvals0Dp = [];
% % Fvals0De = [];
% % Fvals1D = [];
% % timekBAAM = [];
% % timeeBAAM = [];
% % time1D = [];
% % 
% % parfor ii = 1:length(funVals1D(:,1))
% % % parfor ii = 1:5
% %     % for ii = 1:5
% %     parametersT = parameters;
% %     theta = funVals1D(ii,1:6);
% %     % theta(1) = theta(1);
% %     parametersT.pressureDrop = 1;
% %     parametersT.equilibrium = 0;
% % 
% %     KPIs = kBAAM_Outputs_nonIsothermal_dP(parametersT,theta); % KPIS = [-Recovery, -Purity];
% %     kBAAMt = KPIs(end-1);
% %     % KPIs(3)=KPIs(3).*KPIs(4).*(funVals1D(ii,3)+funVals1D(ii,4)+funVals1D(ii,5)+20).*parameters.V_column.*3.600.*0.044;
% % 
% %     % parametersT.pressureDrop = 0;
% %     % KPIs2 = kBAAM_Outputs_nonIsothermal_dP(parametersT,theta); % KPIS = [-Recovery, -Purity];
% %     % KPIs2(3)=KPIs2(3).*KPIs2(4).*(funVals1D(ii,3)+funVals1D(ii,4)+funVals1D(ii,5)+20).*parameters.V_column.*3.600.*0.044;
% % 
% %     tic
% %     KPIs3 = DBAAM_Outputs(parametersT,theta);
% %     eBAAMt = toc;
% %     % KPIs3(3)=KPIs3(3).*KPIs3(4).*(funVals1D(ii,3)+funVals1D(ii,4)+funVals1D(ii,5)+20).*parameters.V_column.*3.600.*0.044;
% % 
% %     Fvals0D = [Fvals0D;KPIs];
% %     Fvals0Dp = [Fvals0Dp;KPIs];
% %     Fvals0De = [Fvals0De;KPIs3];
% %     Fvals1D = [Fvals1D;funVals1D(ii,7:12)];
% %     timekBAAM = [timekBAAM;kBAAMt];
% %     timeeBAAM = [timeeBAAM;eBAAMt];
% %     time1D = [time1D;funVals1D(ii,11)];
% % end
% % z13X.Fvals1D = Fvals1D;
% % z13X.Fvals0Dp = Fvals0Dp;
% % z13X.Fvals0De = Fvals0De;
% % z13X.Fvals0D = Fvals0D;
% % z13X.timekBAAM = timekBAAM;
% % z13X.timeeBAAM = timeeBAAM;
% % z13X.time1D = time1D;
% % save(convertStringsToChars(strcat("calibrationResults_z13X","_",datestr(now,'ddmmyyhhMM'),".mat")),'z13X')
load('calibrationResults_z13X_3001261725.mat')

figure(999)
subplot(1,4,1)
hold on
[P2,S] = plot95confregion(z13X.Fvals1D(:,1),z13X.Fvals0De(:,1),'k')
scatter(z13X.Fvals1D(:,1),z13X.Fvals0De(:,1),25,'ok'); z13X.R20DePu = S.rsquared;
if plot0Dp
[P2,S] = plot95confregion(z13X.Fvals1D(:,1),z13X.Fvals0Dp(:,1),'r')
scatter(z13X.Fvals1D(:,1),z13X.Fvals0Dp(:,1),25,'or')
end
if plot0D
[P2,S] = plot95confregion(z13X.Fvals1D(:,1),z13X.Fvals0D(:,1),'b')
scatter(z13X.Fvals1D(:,1),z13X.Fvals0D(:,1),30,'ob','filled'); z13X.R20DkPu = S.rsquared;
end
plot([0 1e2],[0 1e2],LineStyle="--",Color='r',LineWidth=2)
xlim([75 100])
ylim([75 100])
xlabel('Pu_{CO_{2},1D} [%]'); ylabel('Pu_{CO_{2},0D} [%]')
box on; grid off; set(gca,"LineWidth",2,"FontName","times new roman","FontSize",15)

subplot(1,4,2)
hold on
[P2,S] = plot95confregion(z13X.Fvals1D(:,2),z13X.Fvals0De(:,2),'k')
scatter(z13X.Fvals1D(:,2),z13X.Fvals0De(:,2),25,'ok'); z13X.R20DeRe = S.rsquared;
if plot0Dp
[P2,S] = plot95confregion(z13X.Fvals1D(:,2),z13X.Fvals0Dp(:,2),'r')
scatter(z13X.Fvals1D(:,2),z13X.Fvals0Dp(:,2),25,'or')
end
if plot0D
[P2,S] = plot95confregion(z13X.Fvals1D(:,2),z13X.Fvals0D(:,2),'b')
scatter(z13X.Fvals1D(:,2),z13X.Fvals0D(:,2),30,'ob','filled'); z13X.R20DkRe = S.rsquared;
end
plot([0 1e2],[0 1e2],LineStyle="--",Color='r',LineWidth=2)
xlim([0 100])
ylim([0 100])
xlabel('Re_{CO_{2},1D} [%]'); ylabel('Re_{CO_{2},0D} [%]')
box on; grid off; set(gca,"LineWidth",2,"FontName","times new roman","FontSize",15)

subplot(1,4,3)
hold on
[P2,S] = plot95confregion(z13X.Fvals1D(:,3),z13X.Fvals0De(:,3),'k')
scatter(z13X.Fvals1D(:,3),z13X.Fvals0De(:,3),25,'ok'); z13X.R20DeWC = S.rsquared;
if plot0Dp
[P2,S] = plot95confregion(z13X.Fvals1D(:,3),z13X.Fvals0Dp(:,3),'r')
scatter(z13X.Fvals1D(:,3),z13X.Fvals0Dp(:,3),25,'or')
end
if plot0D
[P2,S] = plot95confregion(z13X.Fvals1D(:,3),z13X.Fvals0D(:,3),'b')
scatter(z13X.Fvals1D(:,3),z13X.Fvals0D(:,3),30,'ob','filled'); z13X.R20DkWC = S.rsquared;
end
plot([0 8e3],[0 8e3],LineStyle="--",Color='r',LineWidth=2)
xlim([0 1e3])
ylim([0 1e3])
xlabel('W_{eq,1D} [kWh/tonne]'); ylabel('W_{eq,0D} [kWh/tonne]')
box on; grid off; set(gca,"LineWidth",2,"FontName","times new roman","FontSize",15)

subplot(1,4,4)
hold on
[P2,S] = plot95confregion(z13X.Fvals1D(:,4),z13X.Fvals0De(:,4),'k')
scatter(z13X.Fvals1D(:,4),z13X.Fvals0De(:,4),25,'ok'); z13X.R20DePr = S.rsquared;
if plot0Dp
[P2,S] = plot95confregion(z13X.Fvals1D(:,4),z13X.Fvals0Dp(:,4),'r')
scatter(z13X.Fvals1D(:,4),z13X.Fvals0Dp(:,4),25,'or')
end
if plot0D
[P2,S] = plot95confregion(z13X.Fvals1D(:,4),z13X.Fvals0D(:,4),'b')
scatter(z13X.Fvals1D(:,4),z13X.Fvals0D(:,4),30,'ob','filled'); z13X.R20DkPr = S.rsquared;
end
plot([0 1e1],[0 1e1],LineStyle="--",Color='r',LineWidth=2)
xlim([0 2])
ylim([0 6])
xlabel('Pr_{CO_{2},1D} [mol/m^{3}s]'); ylabel('Pr_{CO_{2},0D} [mol/m^{3}s]')
box on; grid off; set(gca,"LineWidth",2,"FontName","times new roman","FontSize",15)



set(gcf,"Position",[100,1200,1100 300],"Units","pixels")

%%
figure(666)
subplot(1,4,1)
hold on
histogram(z13X.Fvals1D(:,1),20,'FaceColor','r','FaceAlpha',0.31)
xlim([75 100])
% ylim([60 100])
xlabel('Pu_{CO_{2},1D} [%]'); ylabel('Count')
box on; grid off; set(gca,"LineWidth",2,"FontName","times new roman","FontSize",15)

subplot(1,4,2)
hold on
histogram(z13X.Fvals1D(:,2),20,'FaceColor','r','FaceAlpha',0.31)
xlim([00 100])
% ylim([60 100])
xlabel('Re_{CO_{2},1D} [%]'); ylabel('Count')
box on; grid off; set(gca,"LineWidth",2,"FontName","times new roman","FontSize",15)

subplot(1,4,3)
hold on
histogram(z13X.Fvals1D(:,3),20,'FaceColor','r','FaceAlpha',0.31)
xlim([0 1e3])
% ylim([60 100])
xlabel('W_{eq,1D} [kWh/tonne]'); ylabel('Count')
box on; grid off; set(gca,"LineWidth",2,"FontName","times new roman","FontSize",15)

subplot(1,4,4)
hold on
histogram(z13X.Fvals1D(:,4),20,'FaceColor','r','FaceAlpha',0.31)
xlim([0 2])
% ylim([60 100])
xlabel('Pr_{CO_{2},1D} [mol/m^{3}s]'); ylabel('Count')
box on; grid off; set(gca,"LineWidth",2,"FontName","times new roman","FontSize",15)

set(gcf,"Position",[100,100,1100 300],"Units","pixels")

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
% load('Z13X_AW_2022.mat')
load('CSAC_RH_2014.mat')
% load('funVals1D_13X.mat')
load('funVals1D_AC_new.mat')
load('funvals1D_AC_all.mat')
parameters.adsorbentName = "13XCalibrate_fin";
parameters.outputType = "Opt";
parameters.processType = "VSA";
parameters.modelType = "nonisothermal";
parameters.OptType = "sampling";
funVals1D = sortrows(funVals1D,10,"ascend");
funVals1D(funVals1D(:,9)==0,:) = [];
% funVals1D(funVals1D(:,8)<10,:) = [];


% % funVals1D(:,9) = funVals1D(:,9).*funVals1D(:,10).*(funVals1D(:,3)+funVals1D(:,4)+funVals1D(:,5)+20).*parameters.V_column.*3.600.*0.044;
% Fvals0D = [];
% Fvals0Dp = [];
% Fvals0De = [];
% Fvals1D = [];
% timekBAAM = [];
% timeeBAAM = [];
% time1D = [];
% 
% parfor ii = 1:length(funVals1D(:,1))
% % parfor ii = 1:5
%     % for ii = 1:5
%     parametersT = parameters;
%     theta = funVals1D(ii,1:6);
%     % theta(1) = theta(1);
%     parametersT.pressureDrop = 1;
%     parametersT.equilibrium = 0;
% 
%     KPIs = kBAAM_Outputs_nonIsothermal_dP(parametersT,theta) % KPIS = [-Recovery, -Purity];
%     kBAAMt = KPIs(end-1);
%     % KPIs(3)=KPIs(3).*KPIs(4).*(funVals1D(ii,3)+funVals1D(ii,4)+funVals1D(ii,5)+20).*parameters.V_column.*3.600.*0.044;
% 
%     % parametersT.pressureDrop = 0;
%     % KPIs2 = kBAAM_Outputs_nonIsothermal_dP(parametersT,theta); % KPIS = [-Recovery, -Purity];
%     % KPIs2(3)=KPIs2(3).*KPIs2(4).*(funVals1D(ii,3)+funVals1D(ii,4)+funVals1D(ii,5)+20).*parameters.V_column.*3.600.*0.044;
% 
%     tic
%     KPIs3 = DBAAM_Outputs(parametersT,theta);
%     eBAAMt = toc;
%     % KPIs3(3)=KPIs3(3).*KPIs3(4).*(funVals1D(ii,3)+funVals1D(ii,4)+funVals1D(ii,5)+20).*parameters.V_column.*3.600.*0.044;
% 
%     Fvals0D = [Fvals0D;KPIs];
%     Fvals0Dp = [Fvals0Dp;KPIs];
%     Fvals0De = [Fvals0De;KPIs3];
%     Fvals1D = [Fvals1D;funVals1D(ii,7:12)];
%     timekBAAM = [timekBAAM;kBAAMt];
%     timeeBAAM = [timeeBAAM;eBAAMt];
%     time1D = [time1D;funVals1D(ii,11)];
% end
% ACRB3.Fvals1D = Fvals1D;
% ACRB3.Fvals0Dp = Fvals0Dp;
% ACRB3.Fvals0De = Fvals0De;
% ACRB3.Fvals0D = Fvals0D;
% ACRB3.timekBAAM = timekBAAM;
% ACRB3.timeeBAAM = timeeBAAM;
% ACRB3.time1D = time1D;
% save(convertStringsToChars(strcat("calibrationResults_ACRB3","_",datestr(now,'ddmmyyhhMM'),".mat")),'ACRB3')
load('calibrationResults_ACRB3_3001261725.mat')


% % % % % % 
% % % % % %     Fvals0D = [ACRB3.Fvals0D ACRB3.timekBAAM 12.*ones(length(ACRB3.timekBAAM),1)]
% % % % % %     Fvals0Dp = [ACRB3.Fvals0Dp ACRB3.timekBAAM 12.*ones(length(ACRB3.timekBAAM),1)];
% % % % % %     Fvals0De = [ACRB3.Fvals0De];
% % % % % %     Fvals1D = [ACRB3.Fvals1D ACRB3.time1D 12.*ones(length(ACRB3.timekBAAM),1)];
% % % % % %     timekBAAM = ACRB3.timekBAAM;
% % % % % %     timeeBAAM = ACRB3.timeeBAAM;
% % % % % %     time1D = ACRB3.time1D;
% % % % % % 
% % % % % % load('calibrationResults_ACRB3_3001261522.mat')
% % % % % % 
% % % % % % ACRB3.Fvals1D = [ACRB3.Fvals1D; Fvals1D];
% % % % % % ACRB3.Fvals0Dp = [ACRB3.Fvals0Dp;  Fvals0Dp];
% % % % % % ACRB3.Fvals0De = [ACRB3.Fvals0De;  Fvals0De];
% % % % % % ACRB3.Fvals0D = [ACRB3.Fvals0D;  Fvals0D];
% % % % % % ACRB3.timekBAAM = [ACRB3.timekBAAM;  timekBAAM];
% % % % % % ACRB3.timeeBAAM = [ACRB3.timeeBAAM;  timeeBAAM];
% % % % % % ACRB3.time1D = [ACRB3.time1D;  time1D];

figure(888)
subplot(1,4,1)
hold on
[P2,S] = plot95confregion(ACRB3.Fvals1D(:,1),ACRB3.Fvals0De(:,1),'k')
scatter(ACRB3.Fvals1D(:,1),ACRB3.Fvals0De(:,1),25,'ok'); ACRB3.R20DePu = S.rsquared;
if plot0Dp
[P2,S] = plot95confregion(ACRB3.Fvals1D(:,1),ACRB3.Fvals0Dp(:,1),'r')
scatter(ACRB3.Fvals1D(:,1),ACRB3.Fvals0Dp(:,1),25,'or')
end
if plot0D
[P2,S] = plot95confregion(ACRB3.Fvals1D(:,1),ACRB3.Fvals0D(:,1),'b')
scatter(ACRB3.Fvals1D(:,1),ACRB3.Fvals0D(:,1),30,'ob','filled'); ACRB3.R20DkPu = S.rsquared;
end
plot([0 1e2],[0 1e2],LineStyle="--",Color='r',LineWidth=2)
xlim([60 100])
ylim([60 100])
xlabel('Pu_{CO_{2},1D} [%]'); ylabel('Pu_{CO_{2},0D} [%]')
box on; grid off; set(gca,"LineWidth",2,"FontName","times new roman","FontSize",15)

subplot(1,4,2)
hold on
[P2,S] = plot95confregion(ACRB3.Fvals1D(:,2),ACRB3.Fvals0De(:,2),'k')
scatter(ACRB3.Fvals1D(:,2),ACRB3.Fvals0De(:,2),25,'ok'); ACRB3.R20DeRe = S.rsquared;
if plot0Dp
[P2,S] = plot95confregion(ACRB3.Fvals1D(:,2),ACRB3.Fvals0Dp(:,2),'r')
scatter(ACRB3.Fvals1D(:,2),ACRB3.Fvals0Dp(:,2),25,'or')
end
if plot0D
[P2,S] = plot95confregion(ACRB3.Fvals1D(:,2),ACRB3.Fvals0D(:,2),'b')
scatter(ACRB3.Fvals1D(:,2),ACRB3.Fvals0D(:,2),30,'ob','filled'); ACRB3.R20DkRe = S.rsquared;
end
plot([0 1e2],[0 1e2],LineStyle="--",Color='r',LineWidth=2)
xlim([0 100])
ylim([0 100])
xlabel('Re_{CO_{2},1D} [%]'); ylabel('Re_{CO_{2},0D} [%]')
box on; grid off; set(gca,"LineWidth",2,"FontName","times new roman","FontSize",15)

subplot(1,4,3)
hold on
[P2,S] = plot95confregion(ACRB3.Fvals1D(:,3),ACRB3.Fvals0De(:,3),'k')
scatter(ACRB3.Fvals1D(:,3),ACRB3.Fvals0De(:,3),25,'ok'); ACRB3.R20DeWC = S.rsquared;
if plot0Dp
[P2,S] = plot95confregion(ACRB3.Fvals1D(:,3),ACRB3.Fvals0Dp(:,3),'r')
scatter(ACRB3.Fvals1D(:,3),ACRB3.Fvals0Dp(:,3),25,'or')
end
if plot0D
[P2,S] = plot95confregion(ACRB3.Fvals1D(:,3),ACRB3.Fvals0D(:,3),'b')
scatter(ACRB3.Fvals1D(:,3),ACRB3.Fvals0D(:,3),30,'ob','filled'); ACRB3.R20DkWC = S.rsquared;
end
plot([0 8e3],[0 8e3],LineStyle="--",Color='r',LineWidth=2)
xlim([0 1e3])
ylim([0 1e3])
xlabel('W_{eq,1D} [kWh/tonne]'); ylabel('W_{eq,0D} [kWh/tonne]')
box on; grid off; set(gca,"LineWidth",2,"FontName","times new roman","FontSize",15)

subplot(1,4,4)
hold on
[P2,S] = plot95confregion(ACRB3.Fvals1D(:,4),ACRB3.Fvals0De(:,4),'k')
scatter(ACRB3.Fvals1D(:,4),ACRB3.Fvals0De(:,4),25,'ok'); ACRB3.R20DePr = S.rsquared;
if plot0Dp
[P2,S] = plot95confregion(ACRB3.Fvals1D(:,4),ACRB3.Fvals0Dp(:,4),'r')
scatter(ACRB3.Fvals1D(:,4),ACRB3.Fvals0Dp(:,4),25,'or')
end
if plot0D
[P2,S] = plot95confregion(ACRB3.Fvals1D(:,4),ACRB3.Fvals0D(:,4),'b')
scatter(ACRB3.Fvals1D(:,4),ACRB3.Fvals0D(:,4),30,'ob','filled'); ACRB3.R20DkPr = S.rsquared;
end
plot([0 1e1],[0 1e1],LineStyle="--",Color='r',LineWidth=2)
xlim([0 2])
ylim([0 2])
xlabel('Pr_{CO_{2},1D} [mol/m^{3}s]'); ylabel('Pr_{CO_{2},0D} [mol/m^{3}s]')
box on; grid off; set(gca,"LineWidth",2,"FontName","times new roman","FontSize",15)


set(gcf,"Position",[100,100,1100 300],"Units","pixels")


%%
figure(777)
subplot(1,4,1)
hold on
histogram(ACRB3.Fvals1D(:,1),20,'FaceColor','r','FaceAlpha',0.31)
xlim([60 100])
% ylim([60 100])
xlabel('Pu_{CO_{2},1D} [%]'); ylabel('Count')
box on; grid off; set(gca,"LineWidth",2,"FontName","times new roman","FontSize",15)

subplot(1,4,2)
hold on
histogram(ACRB3.Fvals1D(:,2),20,'FaceColor','r','FaceAlpha',0.31)
xlim([00 100])
% ylim([60 100])
xlabel('Re_{CO_{2},1D} [%]'); ylabel('Count')
box on; grid off; set(gca,"LineWidth",2,"FontName","times new roman","FontSize",15)

subplot(1,4,3)
hold on
histogram(ACRB3.Fvals1D(:,3),20,'FaceColor','r','FaceAlpha',0.31)
xlim([0 1e3])
% ylim([60 100])
xlabel('W_{eq,1D} [kWh/tonne]'); ylabel('Count')
box on; grid off; set(gca,"LineWidth",2,"FontName","times new roman","FontSize",15)

subplot(1,4,4)
hold on
histogram(ACRB3.Fvals1D(:,4),20,'FaceColor','r','FaceAlpha',0.31)
xlim([0 2])
% ylim([60 100])
xlabel('Pr_{CO_{2},1D} [mol/m^{3}s]'); ylabel('Count')
box on; grid off; set(gca,"LineWidth",2,"FontName","times new roman","FontSize",15)

set(gcf,"Position",[100,100,1100 300],"Units","pixels")
%%
figure(88)
b = [ACRB3.timeeBAAM;z13X.timeeBAAM];
a = [ACRB3.timekBAAM;z13X.timekBAAM];
c = [ACRB3.time1D;z13X.time1D];
% c = a.*12+b;
ydata = [ a b c];
hold on
xLabs = categorical(["0-D Kin" "0-D Eq"  "1-D Kin" ]);
x = repmat(xLabs,length(a),1);
boxS = swarmchart(x,ydata,5);
boxK1 = boxchart(x(:,1),ydata(:,1),"BoxFaceColor",'b',"MarkerStyle","none","BoxEdgeColor",'k',"BoxFaceAlpha",0.1);
boxK2 = boxchart(x(:,2),ydata(:,2),"BoxFaceColor",'k',"MarkerStyle","none","BoxEdgeColor",'k',"BoxFaceAlpha",0.1);
boxK3 = boxchart(x(:,3),ydata(:,3),"BoxFaceColor",'r',"MarkerStyle","none","BoxEdgeColor",'k',"BoxFaceAlpha",0.1);
% boxK3 = boxchart(x(:,3),ydata(:,3),"BoxFaceColor",'r',"MarkerStyle","none","BoxEdgeColor",'k',"BoxFaceAlpha",0.1);
boxS(1).Marker = '.';
boxS(1).SizeData = 30;
boxS(1).MarkerEdgeColor = 'b';
boxS(1).MarkerFaceColor = 'b';
boxS(2).SizeData = 30;
boxS(2).Marker = '.';
boxS(2).MarkerEdgeColor = 'k';
boxS(2).MarkerFaceColor = 'k';
boxS(3).SizeData = 30;
boxS(3).Marker = '.';
boxS(3).MarkerEdgeColor = 'r';
boxS(3).MarkerFaceColor = 'r';
plot(x(1,1),median(a),'-','MarkerSize',5,'Marker',"o",'MarkerFaceColor','g','MarkerEdgeColor','g')
plot(x(1,2),median(b),'-','MarkerSize',5,'Marker',"o",'MarkerFaceColor','g','MarkerEdgeColor','g')
plot(x(1,3),median(c),'-','MarkerSize',5,'Marker',"o",'MarkerFaceColor','g','MarkerEdgeColor','g')
mean(a)
mean(b)
mean(c)
ylabel("Simulation time [s]")
ylim([0.05 1000])
box on; grid off; set(gca,"LineWidth",2,"FontName","times new roman","FontSize",15,'YScale','log','GridLineWidth',1)

set(gcf,"Position",[100,100,350 300],"Units","pixels")
%%
% figure(24)
% 
% % meanVals = [mean(b).*ones(length(a),1) mean(a).*ones(length(a),1) mean(c).*ones(length(a),1)];
% % meanVals = [mean(b) mean(a) mean(c)];
% ODEnum = [6 1 6*10];
% meanVals = [ODEnum(2).*ones(length(a),1) ODEnum(1).*ones(length(a),1) ODEnum(3).*ones(length(a),1)];
% 
% % loglog(ODEnum,meanVals)
% 
% ydata = [b a c];
% hold on 
% nVals = logspace(0,3,1000);
% allvals = nVals.^2.*((c-b)./60^2)+b;
% plot(log10(nVals'),allvals','-g','LineWidth',0.1,'Color',[0 0 0 0.1]);
% % plot(log10(meanVals'),allvals','-g','LineWidth',0.1,'Color',[0 0 0 0.1]);
% % plot(log10(meanVals'),allvals','-g','LineWidth',0.1,'Color',[0 0 0 0.1]);
% scatter(log10(meanVals(:,1)'),ydata(:,1)',20,'filled','ok','MarkerFaceAlpha',0.2)
% scatter(log10(meanVals(:,2)'),ydata(:,2)',20,'filled','ob','MarkerFaceAlpha',0.2)
% scatter(log10(meanVals(:,3)'),ydata(:,3)',20,'filled','or','MarkerFaceAlpha',0.2)
% % boxS = swarmchart(log10(meanVals),ydata,5);
% 
% boxK1 = boxchart(log10(meanVals(:,1)),ydata(:,1),"BoxFaceColor",'b',"MarkerStyle","none","BoxEdgeColor",'k',"BoxFaceAlpha",0.1);
% boxK2 = boxchart(log10(meanVals(:,2)),ydata(:,2),"BoxFaceColor",'k',"MarkerStyle","none","BoxEdgeColor",'k',"BoxFaceAlpha",0.1);
% boxK3 = boxchart(log10(meanVals(:,3)),ydata(:,3),"BoxFaceColor",'r',"MarkerStyle","none","BoxEdgeColor",'k',"BoxFaceAlpha",0.1);
% 
% plot(log10(ODEnum(1)),mean(a),'-','MarkerSize',5,'Marker',"o",'MarkerFaceColor','g','MarkerEdgeColor','g')
% plot(log10(ODEnum(2)),mean(b),'-','MarkerSize',5,'Marker',"o",'MarkerFaceColor','g','MarkerEdgeColor','g')
% plot(log10(ODEnum(3)),mean(c),'-','MarkerSize',5,'Marker',"o",'MarkerFaceColor','g','MarkerEdgeColor','g')
% plot(log10(nVals),(nVals).^2.*((mean(c)-mean(b))./60^2)+mean(b),'-g','LineWidth',2)
% 
% % boxS(1).Marker = '.';
% % boxS(1).SizeData = 30;
% % boxS(1).MarkerEdgeColor = 'k';
% % boxS(1).MarkerFaceColor = 'k';
% % boxS(2).SizeData = 30;
% % boxS(2).Marker = '.';
% % boxS(2).MarkerEdgeColor = 'b';
% % boxS(2).MarkerFaceColor = 'b';
% % boxS(3).SizeData = 30;
% % boxS(3).Marker = '.';
% % boxS(3).MarkerEdgeColor = 'r';
% % boxS(3).MarkerFaceColor = 'r';
% 
% box on; grid off; set(gca,"LineWidth",2,"FontName","times new roman","FontSize",15,'YScale','log','XScale','linear','GridLineWidth',1,'XTick',[0 1 2],'YTick',[0.1 1 10 100 1000],'XTickLabel',["1" "10" "100"])
% 
% ylabel("Simulation time [s]")
% xlabel("Number of ODEs [-]")
% 
% set(gcf,"Position",[100,100,347 308],"Units","pixels")
% 
% ylim([0.05 1000])
% xlim([-0.5 2.33])
% 
% figure(243)
% subplot(1,2,1)
% hold on
% plot(log10(nVals),2.^(nVals)-1,'-b','LineWidth',2)
% plot(log10(nVals),(nVals).^2,'-r','LineWidth',2)
% plot(log10(nVals),(nVals.*log(nVals))+1,'-m','LineWidth',2)
% plot(log10(nVals),(log(nVals))+1,'-g','LineWidth',2)
% 
% box on; grid off; set(gca,"LineWidth",2,"FontName","times new roman","FontSize",15,'YScale','log','XScale','linear','GridLineWidth',1,'XTick',[0 1 2 3 4 ],'YTick',[0.1 1 10 100 1000 10000 100000 1000000 10000000],'XTickLabel',["1" "10" "100" "1000" "10000" "100000" "100000"])
% 
% ylabel("Simulation time [s]")
% xlabel("Number of inputs, N [-]")
% 
% set(gcf,"Position",[100,100,347 308],"Units","pixels")
% 
% ylim([1 12000000])
% xlim([0 3])
% 
% subplot(1,2,2)
% hold on
% 
% plot(log10(nVals),(nVals-0.5)./((1+nVals-0.5)),'-k','LineWidth',3)
% % boxS(1).Marker = '.';
% % boxS(1).SizeData = 30;
% % boxS(1).MarkerEdgeColor = 'k';
% % boxS(1).MarkerFaceColor = 'k';
% % boxS(2).SizeData = 30;
% % boxS(2).Marker = '.';
% % boxS(2).MarkerEdgeColor = 'b';
% % boxS(2).MarkerFaceColor = 'b';
% % boxS(3).SizeData = 30;
% % boxS(3).Marker = '.';
% % boxS(3).MarkerEdgeColor = 'r';
% % boxS(3).MarkerFaceColor = 'r';
% 
% box on; grid off; set(gca,"LineWidth",2,"FontName","times new roman","FontSize",15,'YScale','linear','XScale','linear','GridLineWidth',1,'XTick',[0 1 2 3 4 ],'XTickLabel',["1" "10" "100" "1000" "10000" "100000" "100000"])
% 
% ylabel("Model Accuracy [-]")
% xlabel("Number of inputs, N [-]")
% 
% set(gcf,"Position",[100,100,747 308],"Units","pixels")
% 
% ylim([0 1])
% xlim([0 3])

%%
% figure(26)
% 
% meanVals = [mean(b).*ones(length(a),1) mean(a).*ones(length(a),1) mean(c).*ones(length(a),1)];
% meanVals = [mean(b) mean(a) mean(c)];
% z13X.R2vals = [[z13X.R20DePu;z13X.R20DeRe;z13X.R20DeWC;z13X.R20DePr] [z13X.R20DkPu;z13X.R20DkRe;z13X.R20DkWC;z13X.R20DkPr] [1;1;1;1]];
% ACRB3.R2vals = [[ACRB3.R20DePu;ACRB3.R20DeRe;ACRB3.R20DeWC;ACRB3.R20DePr] [ACRB3.R20DkPu;ACRB3.R20DkRe;ACRB3.R20DkWC;ACRB3.R20DkPr] [1;1;1;1]];
% R2vals = [z13X.R2vals;ACRB3.R2vals];
% % z13X.Residuals = [z13X.]
% % R2valse = 1-()
% ODEnum = [6 1 6*10];
% meanVals = [ODEnum(2).*ones(8,1) ODEnum(1).*ones(8,1) ODEnum(3).*ones(8,1)];
% 
% % loglog(ODEnum,meanVals)
% 
% ydata = [b a c];
% 
% ydata = []
% hold on 
% nVals = logspace(0,3,1000);
% allvals = nVals.^2.*((c-b)./60^2)+b;
% % plot(log10(nVals'),allvals','-g','LineWidth',0.1,'Color',[0 0 0 0.1]);
% % plot(log10(meanVals'),allvals','-g','LineWidth',0.1,'Color',[0 0 0 0.1]);
% % plot(log10(meanVals'),allvals','-g','LineWidth',0.1,'Color',[0 0 0 0.1]);
% bar(mean(log10(meanVals(:,1)')),mean(R2vals(:,1))','FaceColor','k','FaceAlpha',0.2,'BarWidth',0.5)
% bar(mean(log10(meanVals(:,2)')),mean(R2vals(:,2))','FaceColor','b','FaceAlpha',0.2,'BarWidth',0.5)
% bar(mean(log10(meanVals(:,3)')),mean(R2vals(:,3))','FaceColor','r','FaceAlpha',0.2,'BarWidth',0.5)
% scatter(log10(meanVals(:,1)'),R2vals(:,1)',60,'filled','ok','MarkerFaceAlpha',1)
% scatter(log10(meanVals(:,2)'),R2vals(:,2)',60,'filled','ob','MarkerFaceAlpha',1)
% scatter(log10(meanVals(:,3)'),R2vals(:,3)',60,'filled','or','MarkerFaceAlpha',1)
% 
% % boxS = swarmchart(log10(meanVals),ydata,5);
% 
% boxK1 = boxchart(log10(meanVals(:,1)),R2vals(:,1),"BoxFaceColor",'b',"MarkerStyle","none","BoxEdgeColor",'k',"BoxFaceAlpha",0.1);
% boxK2 = boxchart(log10(meanVals(:,2)),R2vals(:,2),"BoxFaceColor",'k',"MarkerStyle","none","BoxEdgeColor",'k',"BoxFaceAlpha",0.1);
% boxK3 = boxchart(log10(meanVals(:,3)),R2vals(:,3),"BoxFaceColor",'r',"MarkerStyle","none","BoxEdgeColor",'k',"BoxFaceAlpha",0.1);
% 
% % plot(log10(ODEnum(1)),mean(a),'-','MarkerSize',5,'Marker',"o",'MarkerFaceColor','g','MarkerEdgeColor','g')
% % plot(log10(ODEnum(2)),mean(b),'-','MarkerSize',5,'Marker',"o",'MarkerFaceColor','g','MarkerEdgeColor','g')
% % plot(log10(ODEnum(3)),mean(c),'-','MarkerSize',5,'Marker',"o",'MarkerFaceColor','g','MarkerEdgeColor','g')
% % plot(log10(nVals),(nVals).^2.*((mean(c)-mean(b))./60^2)+mean(b),'-g','LineWidth',2)
% 
% % boxS(1).Marker = '.';
% % boxS(1).SizeData = 30;
% % boxS(1).MarkerEdgeColor = 'k';
% % boxS(1).MarkerFaceColor = 'k';
% % boxS(2).SizeData = 30;
% % boxS(2).Marker = '.';
% % boxS(2).MarkerEdgeColor = 'b';
% % boxS(2).MarkerFaceColor = 'b';
% % boxS(3).SizeData = 30;
% % boxS(3).Marker = '.';
% % boxS(3).MarkerEdgeColor = 'r';
% % boxS(3).MarkerFaceColor = 'r';
% 
% box on; grid off; set(gca,"LineWidth",2,"FontName","times new roman","FontSize",15,'YScale','linear','XScale','linear','GridLineWidth',1,'XTick',[0 1 2],'XTickLabel',["1" "10" "100"])
% 
% ylabel("nRMSE [-]")
% xlabel("Number of ODEs [-]")
% 
% set(gcf,"Position",[100,100,347 308],"Units","pixels")
% 
% % ylim([0.0 1.])
% xlim([-0.5 2.33])


%%
figure(2976)

meanVals = [mean(b).*ones(length(a),1) mean(a).*ones(length(a),1) mean(c).*ones(length(a),1)];
meanVals = [mean(b) mean(a) mean(c)];
z13X.R2vals = [[z13X.R20DePu;z13X.R20DeRe;z13X.R20DeWC;z13X.R20DePr] [z13X.R20DkPu;z13X.R20DkRe;z13X.R20DkWC;z13X.R20DkPr] [1;1;1;1]];
ACRB3.R2vals = [[ACRB3.R20DePu;ACRB3.R20DeRe;ACRB3.R20DeWC;ACRB3.R20DePr] [ACRB3.R20DkPu;ACRB3.R20DkRe;ACRB3.R20DkWC;ACRB3.R20DkPr] [1;1;1;1]];
R2vals = [z13X.R2vals+ACRB3.R2vals]./2;
% z13X.Residuals = [z13X.]
% R2valse = 1-()
ODEnum = [6 1 6*10];
meanVals = [ODEnum(2).*ones(4,1) ODEnum(1).*ones(4,1) ODEnum(3).*ones(4,1)];


xLabs = categorical(["0-D Kin" "0-D Eq"  "1-D Kin" ]);
x = repmat(xLabs,4,1);

% loglog(ODEnum,meanVals)

ydata = [b a c];

ydata = []
hold on 
nVals = logspace(0,3,1000);
allvals = nVals.^2.*((c-b)./60^2)+b;
% plot(log10(nVals'),allvals','-g','LineWidth',0.1,'Color',[0 0 0 0.1]);
% plot(log10(meanVals'),allvals','-g','LineWidth',0.1,'Color',[0 0 0 0.1]);
% plot(log10(meanVals'),allvals','-g','LineWidth',0.1,'Color',[0 0 0 0.1]);
bar("0-D Eq",mean(R2vals(:,1))','FaceColor','k','FaceAlpha',0.2,'BarWidth',0.6)
bar("0-D Kin",mean(R2vals(:,2))','FaceColor','b','FaceAlpha',0.2,'BarWidth',0.6)
bar("1-D Kin",mean(R2vals(:,3))','FaceColor','r','FaceAlpha',0.2,'BarWidth',0.6)
% scatter(x(:,2),R2vals(:,1)',60,'filled','ok','MarkerFaceAlpha',1)
% scatter(x(:,1),R2vals(:,2)',60,'filled','ob','MarkerFaceAlpha',1)
% scatter(x(:,3),R2vals(:,3)',60,'filled','or','MarkerFaceAlpha',1)

% boxS = swarmchart(log10(meanVals),ydata,5);

boxK1 = boxchart(x(:,2),R2vals(:,1),"BoxFaceColor",'k',"MarkerStyle","none","BoxEdgeColor",'k',"BoxFaceAlpha",0.1,"BoxWidth",0.4);
boxK2 = boxchart(x(:,1),R2vals(:,2),"BoxFaceColor",'b',"MarkerStyle","none","BoxEdgeColor",'b',"BoxFaceAlpha",0.1,"BoxWidth",0.4);
boxK3 = boxchart(x(:,3),R2vals(:,3),"BoxFaceColor",'r',"MarkerStyle","none","BoxEdgeColor",'r',"BoxFaceAlpha",0.1,"BoxWidth",0.4);

% plot(log10(ODEnum(1)),mean(a),'-','MarkerSize',5,'Marker',"o",'MarkerFaceColor','g','MarkerEdgeColor','g')
% plot(log10(ODEnum(2)),mean(b),'-','MarkerSize',5,'Marker',"o",'MarkerFaceColor','g','MarkerEdgeColor','g')
% plot(log10(ODEnum(3)),mean(c),'-','MarkerSize',5,'Marker',"o",'MarkerFaceColor','g','MarkerEdgeColor','g')
% plot(log10(nVals),(nVals).^2.*((mean(c)-mean(b))./60^2)+mean(b),'-g','LineWidth',2)

% boxS(1).Marker = '.';
% boxS(1).SizeData = 30;
% boxS(1).MarkerEdgeColor = 'k';
% boxS(1).MarkerFaceColor = 'k';
% boxS(2).SizeData = 30;
% boxS(2).Marker = '.';
% boxS(2).MarkerEdgeColor = 'b';
% boxS(2).MarkerFaceColor = 'b';
% boxS(3).SizeData = 30;
% boxS(3).Marker = '.';
% boxS(3).MarkerEdgeColor = 'r';
% boxS(3).MarkerFaceColor = 'r';

box on; grid off; set(gca,"LineWidth",2,"FontName","times new roman","FontSize",15,'YScale','linear','XScale','linear','GridLineWidth',1,'XTickLabel',xLabs)

ylabel("1 – nRMSE [-]")
% xlabel("Number of ODEs [-]")

set(gcf,"Position",[100,100,347 308],"Units","pixels")

ylim([0.0 1.1])
yline(1,'--r','LineWidth',0.5)
% xlim([-0.5 2.33])


