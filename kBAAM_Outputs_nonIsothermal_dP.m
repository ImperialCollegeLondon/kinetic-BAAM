%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Imperial College London, United Kingdom
% Multiphase Systems Laboratory
% Year:     2025
% MATLAB:   R2024a
% Authors:  Hassan Azzan (HA)
%
% Purpose:
% Solves non-isothermal kinetic batch adsorber analogue model (k-BAAM) and outputs Purity
% and recovery of the heavy product
%
% Last modified:
% - 2025-12-18, HA: Add total material balance and pressure drop
% - 2025-10-09, HA: Add wall energy balance
% - 2025-10-08, HA: Add reverse engineering optimization method
% - 2025-09-17, HA: Initial creation
%
% Input arguments:
%   - parameters: contains adsorbent properties and process parameters
%   - varargin: input 1 as a variable input argument to produce steady
%               state profile plots
%
% Output arguments:
%   - KPIS: heavy product recovery and purity
%
% Dependencies:
%   - kBAAM_ODEs_nonIsothermal_ND_dP.m
%   - DSL.m
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [KPIs] = kBAAM_Outputs_nonIsothermal_dP(parameters,varargin)
tic
if ~isfield(parameters,'pressureDrop')
    parameters.pressureDrop = 1;
end

if ~isfield(parameters,'equilibrium')
    parameters.equilibrium = 0;
end

if ~isfield(parameters,'cCSTR')
    parameters.cCSTR = 0;
end


if ~isfield(parameters,'testBT')
    parameters.testBT = 0;
end


if ~isfield(parameters,'testEvac')
    parameters.testEvac = 0;
end

if ~isfield(parameters,'normPlot')
    parameters.normPlot = 0;
end


if ~isfield(parameters,'forwardEvac')
    parameters.forwardEvac = 0;
end


if ~isfield(parameters,'SSLSTA')
    parameters.SSLSTA = 0;
end

if ~isfield(parameters,'plot0D')
    parameters.plot0D = 1;
end


if ~isfield(parameters,'rigid')
    parameters.rigid = 1;
end


%% Time durations of adsorption, blowdown, evacuation and pressurization steps
if nargin == 2
    outputType = size(varargin{1});
    theta = varargin{1};
end

try
    if parameters.outputType == "opt"
        theta = parameters.xRef.*theta;
    end
catch
end

if nargin == 2
    if parameters.processType == "VSA"
        parameters.v_in = theta(1);
        parameters.p_I = theta(2);
        parameters.t_ads = theta(3);
        parameters.t_blo = theta(4);
        parameters.t_evac = theta(5);
        parameters.p_L =  theta(6);
        % if parameters.outputType == "opt"
        %     parameters.p_L =  10.^parameters.p_L;
        %     parameters.p_I =  10.^parameters.p_I;
        % end
    elseif parameters.processType == "PVSA"
        parameters.v_in = theta(1);
        parameters.p_I = theta(2);
        parameters.t_ads = theta(3);
        parameters.t_blo = theta(4);
        parameters.t_evac = theta(5);
        parameters.p_L = theta(6);
        parameters.p_H = theta(7);
        % if parameters.outputType == "opt"
        %     parameters.p_L =  10.^parameters.p_L;
        %     parameters.p_I =  10.^parameters.p_I;
        %     parameters.p_H =  10.^parameters.p_H;
        %     % parameters.v_in = parameters.v_in.*parameters.p_H./1e5;
        % end
    elseif parameters.processType == "AdsorbentVSA"
        parameters.v_in = 0.4;
        parameters.p_I = theta(1);
        parameters.t_ads = theta(2);
        parameters.t_blo = theta(3);
        parameters.t_evac = 800;
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
        if parameters.outputType == "opt"
            parameters.p_I = 10.^parameters.p_I;
        end
    elseif parameters.processType == "AdsorbentPVSA"
        parameters.v_in = 0.4;
        parameters.p_I = theta(1);
        parameters.t_ads = theta(2);
        parameters.t_blo = theta(3);
        parameters.t_evac = 2000;
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
        parameters.p_H = theta(8);
        if parameters.outputType == "opt"
            parameters.p_I = 10.^parameters.p_I;
            parameters.p_H = 10.^parameters.p_H;
            % parameters.F_in = parameters.F_in.*parameters.p_H./1e5;
        end
    elseif parameters.processType == "AdsorbentVSAb0"
        parameters.v_in = 0.4;
        parameters.p_I = theta(1);
        parameters.t_ads = theta(2);
        parameters.t_blo = theta(3);
        parameters.t_evac = 800;
        parameters.qsb_1 = theta(4);
        parameters.qsb_2 = theta(4);
        parameters.qsd_1 = 0;
        parameters.qsd_2 = 0;
        parameters.bo_1 = theta(5);
        parameters.bo_2 = theta(6);
        parameters.do_1 = 0;
        parameters.do_2 = 0;
        parameters.delUb_1 = -theta(7);
        parameters.delUb_2 = -theta(7);
        parameters.delUd_1 = 0;
        parameters.delUd_2 = 0;
        if parameters.outputType == "opt"
            parameters.p_I = 10.^parameters.p_I;
        end
    elseif parameters.processType == "AdsorbentPVSAb0"
        parameters.v_in = 0.4;
        parameters.p_I = theta(1);
        parameters.t_ads = theta(2);
        parameters.t_blo = theta(3);
        parameters.t_evac = 2000;
        parameters.qsb_1 = theta(4);
        parameters.qsb_2 = theta(4);
        parameters.qsd_1 = 0;
        parameters.qsd_2 = 0;
        parameters.bo_1 = theta(5);
        parameters.bo_2 = theta(6);
        parameters.do_1 = 0;
        parameters.do_2 = 0;
        parameters.delUb_1 = -theta(7);
        parameters.delUb_2 = -theta(7);
        parameters.delUd_1 = 0;
        parameters.delUd_2 = 0;
        parameters.p_H = theta(8);
        if parameters.outputType == "opt"
            parameters.p_I = 10.^parameters.p_I;
            parameters.p_H = 10.^parameters.p_H;
            % parameters.F_in = parameters.F_in.*parameters.p_H./1e5;
        end
    elseif parameters.processType == "Resin"
        parameters.v_in = theta(6);
        parameters.y1_in = 0.0004;
        parameters.p_I = theta(1);
        parameters.p_L = theta(5);
        parameters.t_ads = theta(2);
        parameters.t_blo = theta(3);
        % parameters.qsb_1 = theta(4);
        % parameters.qsb_2 = parameters.qsb_1;
        parameters.qsb_2 = 5.3446e-02.*parameters.qsb_1./2.38;
        parameters.t_evac = theta(4);
        parameters.qsd_1 = 0;
        parameters.qsd_2 = 0;
        if ~parameters.fixResins
            % parameters.bo_1 = 1.21e-16.*parameters.qsb_1;
            parameters.bo_1 = 1.06e-16.*parameters.qsb_1; % correlation including Lewatit
        end
        parameters.bo_2 = 1.0137e-05;
        parameters.do_1 = 0;
        parameters.do_2 = 0;
        parameters.delUb_1 = -100e3;
        parameters.delUb_2 = -1.3912e+04;
        parameters.delUd_1 = 0;
        parameters.delUd_2 = 0;
        parameters.heating = 1;
        parameters.Theat = 373;
        parameters.rho_s = 1123;
        parameters.cp_s = 1300;
        if parameters.outputType == "opt"
            parameters.p_I = 10.^parameters.p_I;
            parameters.p_L = 10.^parameters.p_L;
            parameters.p_H = 1e5;
            % parameters.F_in = parameters.F_in.*parameters.p_H./1e5;
        end
        % parameters.LDF = -0.0298*parameters.qsb_1+0.0509;
        if ~parameters.fixResins
            % parameters.LDF = 0.0717.*exp(-1.316.*parameters.qsb_1);
            parameters.LDF = 0.132.*exp(-2.076.*parameters.qsb_1); % correlation including Lewatit
        end
    end
end

Rg = 8.3145;

if parameters.processType == "Resin"
    dt = 0.04;
else
    dt = 0.05; % [s]
end

if parameters.heating
    % parameters.Theat = 493;
    parameters.heatPowerDensity = 5e3; % [W/m2]
    parameters.h_in = 8.6*3;
    parameters.Theat = parameters.Theat+6;
else
    parameters.Theat = 0;
end

t_ads   = 0:dt:parameters.t_ads; % time vector for adsorption step [s]
t_blo   = 0:dt:parameters.t_blo; % time vector for blowdown step [s]
t_evac  = 0:dt:parameters.t_evac; % time vector for evacuation step [s]
t_press = 0:dt:parameters.t_press; % time vector for pressurization step [s]

parameters.r_in = (parameters.V_column./(parameters.Lbyr.*pi)).^(1./3); % inner radius of column [m]
% try
%     parameters.r_out = parameters.r_in + parameters.t_wall; % outer radius of column [m]
% catch
    parameters.r_out = parameters.r_in + 0.0175; % outer radius of column [m]
% end
parameters.v_in = parameters.v_in./parameters.e_bed; % superficial to interstitial
parameters.L = parameters.Lbyr.*parameters.r_in; % length of column [m]
parameters.A_in = parameters.r_in.^2.*pi; % inner cross sectional area of column [m2]
parameters.volFlowin = parameters.v_in.*parameters.A_in.*parameters.e_bed; % inlet volumetric flowrate [m3/s]
if ~parameters.pressureDrop
    parameters.deltaP = 150./4.*1./parameters.rp.^2.*((1-parameters.e_bed)./parameters.e_bed).^2.*1.72e-5.*parameters.v_in.*parameters.L./2;
else
    parameters.deltaP = 0;
end
% parameters.p_H = parameters.p_H;
parameters.F_in = parameters.volFlowin.*(parameters.p_H+parameters.deltaP)./(Rg.*parameters.T_feed); % inlet molar flowrate (ideal gas) [mol/s]

parameters.refVals = [1,(parameters.qsb_1+parameters.qsd_1), (parameters.qsb_1+parameters.qsd_1), parameters.T_feed, parameters.T_feed, parameters.p_H]; % vector of reference values of state variables for non-dimensionalization

% Reference values for non-dimensionalization
% parameters.volFluxRef = parameters.F_in./parameters.V_column; % reference molar flowrate per unit volume [mol/m3s]
% parameters.timeRef = 1e5./(Rg.*parameters.T_feed.*parameters.volFluxRef);  % reference time [s]
parameters.timeRef = parameters.V_column./parameters.volFlowin;  % reference time [s]
parameters.TRef = parameters.T_feed; % reference temperature [K]
parameters.TwRef = parameters.T_feed; % reference temperature [K]
parameters.PRef = parameters.p_H; % reference pressure [Pa]


% pressure profiles and derivatives for overall material balance
parameters.P_ads = @(t)parameters.p_H+parameters.deltaP; % pressure profile for adsorption (constant)
%
% parameters.lambdab = 0.5;  % rate constant for vacuum pump, lambda [1/s]
% parameters.lambdae = 0.1;  % rate constant for vacuum pump, lambda [1/s]
parameters.P_initH = parameters.P_ads(parameters.t_ads); % pressure at the end of ads
parameters.P_blo = @(t)parameters.p_I+(parameters.P_initH-parameters.p_I)*exp(-parameters.lambda*t); % pressure profile for blowdown
parameters.P_initL = parameters.P_blo(parameters.t_blo); % pressure at the end of blowdown
parameters.P_evac = @(t)parameters.p_L+(parameters.P_initL-parameters.p_L)*exp(-parameters.lambda*t); % pressure profile for evacuation
parameters.P_initR = parameters.P_evac(parameters.t_evac); % pressure at the end of evacuation
parameters.P_press = @(t)parameters.P_initH+(parameters.P_initR-parameters.P_initH)*exp(-parameters.lambda*t); % pressure profile for pressurization
parameters.dPdt_blo = @(t)-parameters.lambda*(parameters.P_initH-parameters.p_I)*exp(-parameters.lambda*t); % time derivative of pressure profile for blowdown
parameters.dPdt_evac =  @(t)-parameters.lambda*(parameters.P_initL-parameters.p_L)*exp(-parameters.lambda*t); % time derivative of pressure profile for evacuation
parameters.dPdt_press = @(t)-parameters.lambda*(parameters.P_initR-parameters.P_initH)*exp(-parameters.lambda*t); % time derivative of pressure profile for pressurization

%% Initial condition for matrix of solution states
% y1Init = [parameters.y1_in]; % initial mole fraction of component 1 in bed [-]
% y1Init = [1e-4]; % initial mole fraction of component 1 in bed [-]
%% Cyclic steady state simulation
y1Init = [0.99]; % initial mole fraction of component 1 in bed [-]
Tinit = parameters.T_feed;
if ~parameters.SSLSTA
    parameters.qRef = parameters.qsb_1+parameters.qsd_1; % reference adsorbed amount [mol/kg]
    [q1Init, q2Init] = DSL(parameters.p_L, y1Init, parameters.T_feed, parameters.qsb_1, parameters.qsd_1, parameters.qsb_2, parameters.qsd_2, parameters.bo_1, parameters.do_1, parameters.bo_2, parameters.do_2, parameters.delUb_1, parameters.delUd_1, parameters.delUb_2, parameters.delUd_2); % initial adsorbed amounts in bed [mol/kg]
else
    if parameters.rigid
        parameters.SSLSTA1(4:6) = parameters.SSLSTA1(1:3);
        parameters.SSLSTA2(4:6) = parameters.SSLSTA2(1:3);
    end
    parameters.qRef =  parameters.SSLSTA1(4); % reference adsorbed amount [mol/kg]
    parameters.refVals(2) = parameters.qRef ;
    [q1Init, q2Init] = SSLSTA(parameters.p_L, y1Init, parameters.T_feed, parameters); % initial adsorbed amounts in bed [mol/kg]
end
X0 = [y1Init; q1Init; q2Init; parameters.T_feed; parameters.T_feed;parameters.p_L]./parameters.refVals'; % dimensionless vector of initial states
temp_check = zeros(5, 1);

cycle = 0; % initialize cycle number
warning('off','all')
tic

if parameters.testBT
    max_no_Cycles = 1;
else
    max_no_Cycles = 100; % Maximum number of cycles to run to test CSS
end
% parameters.P_initH = parameters.p_H+2*150./4.*1./parameters.rp.^2.*((1-parameters.e_bed)./parameters.e_bed).^2.*1.72e-5.*parameters.v_in.*parameters.L./2;
 parameters.y1_LPP = 0.02; % initialize LP composition
try
while  cycle < max_no_Cycles && mean(temp_check) < 1  && parameters.p_H > parameters.p_I && parameters.p_I > parameters.p_L

    options = odeset('RelTol', 1e-6, 'AbsTol', 1e-6, 'MaxOrder', 2); % sets the levels of relative tolerance, absolute tolerance and maximum order of the ode15 for our system
    % optionsPres = odeset('RelTol', 0.8e-4, 'AbsTol', 0.8e-4, 'MaxOrder', 2); % sets the levels of relative tolerance, absolute tolerance and maximum order of the ode15 for our system
    cycle = cycle+1;


    parameters.loadingFraction = 1;
    parameters.P_press = @(t)parameters.P_initH+(parameters.P_initR-parameters.P_initH)*exp(-parameters.lambda*t); % pressure profile for pressurization
    parameters.dPdt_press = @(t)-parameters.lambda*(parameters.P_initR-parameters.P_initH)*exp(-parameters.lambda*t); % time derivative of pressure profile for pressurization

    [t4, X4] = ode15s(@(t,X) kBAAM_ODEs_nonIsothermal_ND_dP(t,X,parameters,'pres'), t_press./(parameters.timeRef), X0, options);%t1 is the time point at which the solution is evaluated, X1 is the solution states for adsorption step
    t4 = t4.*parameters.timeRef;
    X4(X4<0) = 0;
    X4(X4(:,1)>1,1) = 1;
    if parameters.testBT && ~parameters.testEvac
        parameters.P_initH = parameters.p_H+150./4.*1./parameters.rp.^2.*((1-parameters.e_bed)./parameters.e_bed).^2.*1.72e-5.*parameters.v_in.*parameters.L./1;
        y1Init = [0.000001]; % initial mole fraction of component 1 in bed [-]
        Tinit = parameters.T_feed;

        if ~parameters.SSLSTA
            [q1Init, q2Init] = DSL(parameters.P_initH, y1Init, parameters.T_feed, parameters.qsb_1, parameters.qsd_1, parameters.qsb_2, parameters.qsd_2, parameters.bo_1, parameters.do_1, parameters.bo_2, parameters.do_2, parameters.delUb_1, parameters.delUd_1, parameters.delUb_2, parameters.delUd_2); % initial adsorbed amounts in bed [mol/kg]
        else
            [q1Init, q2Init] = SSLSTA(parameters.P_initH, y1Init, parameters.T_feed, parameters); % initial adsorbed amounts in bed [mol/kg]
        end
        X0 = [y1Init; q1Init; q2Init; parameters.T_feed; parameters.T_feed;parameters.P_initH]./parameters.refVals'; % dimensionless vector of initial states
    else
        X0 = X4(end,:)'; %This sets up the initial condition for the next step, by taking the final state from the previous step.
    end
    X0Evac = X4(end,:)';
    X4 = X4.*parameters.refVals;
    Tinit = X4(end,4);
    if ~parameters.SSLSTA
        [q1max, q2max] = DSL(X0(6).*parameters.PRef, parameters.y1_in, X0(4).*parameters.TRef, parameters.qsb_1, parameters.qsd_1, parameters.qsb_2, parameters.qsd_2, parameters.bo_1, parameters.do_1, parameters.bo_2, parameters.do_2, parameters.delUb_1, parameters.delUd_1, parameters.delUb_2, parameters.delUd_2); % initial adsorbed amounts in bed [mol/kg]
    else
        [q1max, q2max] = SSLSTA(X0(6).*parameters.PRef, parameters.y1_in, X0(4).*parameters.TRef, parameters); % initial adsorbed amounts in bed [mol/kg]
    end
    parameters.q1init = X0(2).*parameters.qRef ;
    parameters.y1init = X0(1) ;
    parameters.q2init = X0(3).*parameters.qRef./q2max;



    [t1, X1] = ode15s(@(t,X) kBAAM_ODEs_nonIsothermal_ND_dP(t,X,parameters,'ads'), t_ads./(parameters.timeRef), X0, options); %t1 is the time point at which the solution is evaluated, X1 is the solution states for adsorption step
    t1 = t1.*parameters.timeRef;
    X1(X1<0) = 0;
    X1(X1(:,1)>1,1) = 1;
    X0 = X1(end,:)';
    X1 = X1.*parameters.refVals;
    if parameters.processType == "Resin"
        if t1(end) < 0.95*parameters.t_ads
            parameters.t_ads = t1(end); % Update evacuation time to length of simulation instead of discarding. Integration seems to fail when flowrate in close to 0 and mole fraction close to 1.
            t_ads  = 0:dt:parameters.t_ads;
            % cycle = 1;
        end
    end

    % if parameters.pressureDrop
    %     parameters.P_initH = X1(end,6);
    % end

    [q1max, q2max] = DSL(X1(end,6), parameters.y1_in, X1(end,4), parameters.qsb_1, parameters.qsd_1, parameters.qsb_2, parameters.qsd_2, parameters.bo_1, parameters.do_1, parameters.bo_2, parameters.do_2, parameters.delUb_1, parameters.delUd_1, parameters.delUb_2, parameters.delUd_2); % initial adsorbed amounts in bed [mol/kg]
    parameters.loadingFraction = (X1(end,2)-X1(1,2))./((q1max-X1(1,2)));
    % parameters.loadingFraction = 1;
    parameters.q1init = X1(end,2);
    parameters.q2init = X1(end,3);
    parameters.y1init = X1(end,1) ;
    Fin_ads = parameters.volFlowin.*(2.*X1(:,6)-parameters.p_H)./(Rg.*parameters.T_feed);
    Fout_ads = Fin_ads - (1 - parameters.e_bed) * parameters.V_column * parameters.rho_s * (gradient(X1(:,2),dt) + gradient(X1(:,3),dt)) + (parameters.e_bed .* X1(:,6)./ (Rg.*X1(:,4).^2)) .* gradient(X1(:,4),dt) .* parameters.V_column - (parameters.e_bed ./ (Rg.*X1(:,4))) .* gradient(X1(:,6),dt) .* parameters.V_column;
    % Fout_ads(Fout_ads<0) = 0;
    F_1_out_ads = Fout_ads.*X1(:,1);
    mol_1_out_ads = trapz(t1,F_1_out_ads); moltot_out_ads = trapz(t1,Fout_ads);
    parameters.y1_LPP = mol_1_out_ads./moltot_out_ads;

    [t2, X2] = ode15s(@(t,X) kBAAM_ODEs_nonIsothermal_ND_dP(t,X,parameters,'blo'), t_blo./(parameters.timeRef), X0, options);%t1 is the time point at which the solution is evaluated, X1 is the solution states for adsorption step
    t2 = t2.*parameters.timeRef;
    X2(X2<0) = 0;
    X2(X2(:,1)>1,1) = 1;
    X0 = X2(end,:)';
    X2 = X2.*parameters.refVals;


    % if ~parameters.forwardEvac
    % else
    %     [q1max, q2max] = DSL(X1(end,6), parameters.y1_in, X1(end,4), parameters.qsb_1, parameters.qsd_1, parameters.qsb_2, parameters.qsd_2, parameters.bo_1, parameters.do_1, parameters.bo_2, parameters.do_2, parameters.delUb_1, parameters.delUd_1, parameters.delUb_2, parameters.delUd_2); % initial adsorbed amounts in bed [mol/kg]
    %     parameters.loadingFraction = (X1(end,2)-X1(1,2))./((q1max-X1(1,2)));
    % end


    % if parameters.testEvac && parameters.testBT
    % X0=X0Evac;
    % end
    [t3, X3] = ode15s(@(t,X) kBAAM_ODEs_nonIsothermal_ND_dP(t,X,parameters,'evac'), t_evac./(parameters.timeRef), X0, options) ;%t1 is the time point at which the solution is evaluated, X1 is the solution states for adsorption step
    t3 = t3.*parameters.timeRef;
    X3(X3<0) = 0;
    X3(X3(:,1)>1,1) = 1;
    X0 = X3(end,:)'; %This sets up the initial condition for the next step, by taking the final state from the previous step.%end means the last row of X1, containing all the state variables at %the t final , : takes all the columns (for all types of the state%variables
    X3 = X3.*parameters.refVals;
    % parameters.loadingFraction = 1;

    if parameters.processType == "Resin"
        if t3(end) < 0.95*parameters.t_evac
            parameters.t_evac = t3(end); % Update evacuation time to length of simulation instead of discarding. Integration seems to fail when flowrate in close to 0 and mole fraction close to 1.
            t_evac  = 0:dt:parameters.t_evac;
            parameters.P_initR = parameters.P_evac(parameters.t_evac);
            % cycle = 1;
        end
    end

    n_1_ads = (X1(end,1).*X1(end,6) * parameters.V_column * parameters.e_bed / (Rg * X1(end,4))) + X1(end,2)* parameters.V_column * (1-parameters.e_bed).*parameters.rho_s ;
    n_1_bd =   (X2(end,1) .*X2(end,6) * parameters.V_column * parameters.e_bed / (Rg * X2(end,4))) + X2(end,2)* parameters.V_column * (1-parameters.e_bed).*parameters.rho_s ;
    n_1_evac = (X3(end,1) .*X3(end,6) * parameters.V_column * parameters.e_bed / (Rg * X3(end,4))) + X3(end,2)* parameters.V_column * (1-parameters.e_bed).*parameters.rho_s ;
    n_2_bd = ((1-X2(end,1)).*X2(end,6) * parameters.V_column * parameters.e_bed / (Rg * X2(end,4))) + X2(end,3)* parameters.V_column * (1-parameters.e_bed).*parameters.rho_s ;
    n_2_evac = ((1-X3(end,1)).*X3(end,6) * parameters.V_column * parameters.e_bed / (Rg * X3(end,4))) + X3(end,3)* parameters.V_column * (1-parameters.e_bed).*parameters.rho_s ;
    n_1_pres = (X4(end,1).*X4(end,6) * parameters.V_column * parameters.e_bed / (Rg * X4(end,4))) + X4(end,2)* parameters.V_column * (1-parameters.e_bed).*parameters.rho_s ;
    n_1_presInit = (X4(1,1).*X4(1,6) * parameters.V_column * parameters.e_bed / (Rg * X4(1,4))) + X4(1,2)* parameters.V_column * (1-parameters.e_bed).*parameters.rho_s ;
    n_1_adsin = trapz(t1,parameters.volFlowin.*(2.*X1(:,6)-parameters.p_H)./(8.314.*parameters.T_feed) * parameters.y1_in);


    cycle_time = (parameters.t_ads + parameters.t_blo + parameters.t_evac + parameters.t_press);

    productivity = (n_1_bd - n_1_evac) /(parameters.V_column.*cycle_time);

    %%  Energy Calculation
    Fout_bd = 0 - (1 - parameters.e_bed) .* parameters.V_column * parameters.rho_s .* (gradient(X2(:,2),dt) + gradient(X2(:,3),dt))-(parameters.e_bed ./ (Rg.*X2(:,4))) .* gradient(X2(:,6),dt) .* parameters.V_column - (parameters.e_bed .* X2(:,6)./ (Rg.*X2(:,4).^2)) .* gradient(X2(:,4),dt) .* parameters.V_column ;
    Fout_bd(Fout_bd<0) = 0;
    Fout_evac = 0 - (1 - parameters.e_bed) .* parameters.V_column * parameters.rho_s .* (gradient(X3(:,2),dt) + gradient(X3(:,3),dt))-(parameters.e_bed ./ (Rg.*X3(:,4))) .* gradient(X3(:,6),dt) .* parameters.V_column - (parameters.e_bed .* X3(:,6)./ (Rg.*X3(:,4).^2)) .* gradient(X3(:,4),dt) .* parameters.V_column ;
    Fin_pres = (1 - parameters.e_bed) .* parameters.V_column * parameters.rho_s .*  (gradient(X4(:,2),dt) + gradient(X4(:,3),dt))  + (parameters.e_bed ./ (Rg.*X4(:,4))) .* gradient(X4(:,6),dt) .* parameters.V_column + (parameters.e_bed .* X4(:,6)./ (Rg.*X4(:,4).^2)) .* gradient(X4(:,4),dt) .* parameters.V_column ;  %calculates the inlet flow rate of the pressurization step

    if parameters.pressType == "LPP"
    F_1_in_pres = Fin_pres.*parameters.y1_LPP;
    else
    F_1_in_pres = Fin_pres.*parameters.y1_in;
    end
    mol_1_in_pres = trapz(t4,F_1_in_pres); 

    F_1_out_bd = Fout_bd.*X2(:,1);
    mol_1_out_bd = trapz(t2,F_1_out_bd); 

    F_1_out_evac = Fout_evac.*X3(:,1);
    mol_1_out_evac = trapz(t3,F_1_out_evac); 


    vin_pres = Fin_pres.*Rg.*parameters.T_feed./(parameters.P_press(t4))./parameters.A_in./parameters.e_bed;
    vout_evac = Fout_evac.*Rg.*X3(:,4)./(parameters.P_evac(t3))./parameters.A_in./parameters.e_bed;
    vout_bd = Fout_bd.*Rg.*X2(:,4)./(parameters.P_blo(t2))./parameters.A_in./parameters.e_bed;
    % vout_evac = Fout_evac.*Rg.*X3(:,4)./(X3(:,6))./parameters.A_in./parameters.e_bed;
    % vout_bd = Fout_bd.*Rg.*X2(:,4)./(X2(:,6))./parameters.A_in./parameters.e_bed;

    voutHalf_bd = (-2./(parameters.L)).*(4./150./1.72e-5).*(parameters.e_bed./(1 - parameters.e_bed)).^2.*parameters.rp.^2.*(parameters.P_blo(t2)-X2(:,6));
    voutHalf_evac = (-2./(parameters.L)).*(4./150./1.72e-5).*(parameters.e_bed./(1 - parameters.e_bed)).^2.*parameters.rp.^2.*(parameters.P_evac(t3)-X3(:,6));


    t_cycle = [t1; t2 + t1(end); t3 + t1(end) + t2(end); t4 + t1(end) + t2(end) + t3(end)]; %Shift the t2 time vector forward in time so it starts immediately after t1 ends
    X_cycle = [X1; X2; X3; X4];
    F_cycleOut = [Fout_ads;Fout_bd;Fout_evac;zeros(length(t4),1)];
    F_cycleIn = [Fin_ads;zeros(length(t2),1);zeros(length(t3),1);Fin_pres];

    % recovery_percentage = 100 * (n_1_bd - n_1_evac) / ((n_1_ads - n_1_evac + mol_1_out_ads));
    % recovery_percentage = 100 * (n_1_bd - n_1_evac) / (n_1_pres - n_1_evac + n_1_adsin);
    if parameters.pressType == "LPP"
        % mole_LP_recycle = n_1_pres-n_1_presInit;
        mole_LP_recycle = mol_1_in_pres;
    else
        mole_LP_recycle = 0;
    end
    recovery_percentage = 100 * (n_1_bd - n_1_evac) / (mol_1_out_bd + mol_1_out_evac + mol_1_out_ads-mole_LP_recycle);
    % recovery_percentage = 100 * (n_1_bd - n_1_evac) / ((n_1_ads - n_1_evac + mol_1_out_ads-mole_LP_recycle));
    % recovery_percentage = 100 * (n_1_bd - n_1_evac) / (mol_1_in_pres + n_1_adsin);

    purity_percentage = 100 * (n_1_bd - n_1_evac) / (n_1_bd - n_1_evac + max(0,n_2_bd - n_2_evac));

    % if round(recovery_percentage,1)>100 || round(purity_percentage,1)>100
    % if round(recovery_percentage,1)>100  
        % recovery_percentage = 0;
        % purity_percentage = 0;
        skipFlag = 0;
    % else
    %     skipFlag = 0;
    % end


    eta_bd = 0.8.*(19.55.*parameters.P_blo(t2).*1e-5./(1+19.55.*parameters.P_blo(t2).*1e-5));
    eta_evac = 0.8.*(19.55.*parameters.P_evac(t3).*1e-5./(1+19.55.*parameters.P_evac(t3).*1e-5));
    % eta_bd = 0.8.*(19.55.*parameters.p_I.*1e-5./(1+19.55.*parameters.p_I.*1e-5));
    % eta_evac = 0.8.*(19.55.*parameters.p_L.*1e-5./(1+19.55.*parameters.p_L.*1e-5));
    eta_press = 0.72;
    eta_ads = 0.72;



    if parameters.heating
        heatFlag = X3(:,4) < parameters.Theat;
        Qheat = trapz(t3,parameters.heatPowerDensity.*(parameters.Theat-X3(:,4))./(parameters.Theat-parameters.T_feed).*parameters.r_out.*2.*parameters.L.*heatFlag); % external heat flux if heating is used
        EC_HEAT = Qheat;
    else
        EC_HEAT = 0;
    end
    EC_BD   = trapz(t2,1./eta_bd    .*vout_bd          .*parameters.A_in.*parameters.e_bed.*parameters.P_blo(t2).*(1.4./0.4).*((1e5./min(1e5,parameters.P_blo(t2))).^(0.4./1.4)-1));
    EC_EVAC = trapz(t3,1./eta_evac  .*vout_evac        .*parameters.A_in.*parameters.e_bed.*parameters.P_evac(t3).*(1.4./0.4).*((1e5./min(1e5,parameters.P_evac(t3))).^(0.4./1.4)-1));
    % EC_FAN = trapz(t1, 1./eta_ads   .*parameters.v_in  .*parameters.A_in.*parameters.e_bed.*X1(:,6).*(1.4./0.4).*((max(1e5,X1(:,6))./1e5).^(0.4./1.4)-1));
    EC_FAN = trapz(t1, 1./eta_ads   .*parameters.v_in  .*parameters.A_in.*parameters.e_bed.*(2.*X1(:,6)-parameters.p_H).*(1.4./0.4).*((max(1e5,(2.*X1(:,6)-parameters.p_H))./1e5).^(0.4./1.4)-1));
    EC_PRES = trapz(t4,1./eta_press .*vin_pres         .*parameters.A_in.*parameters.e_bed.*parameters.P_press(t4).*(1.4./0.4).*((max(1e5,parameters.P_press(t4))./1e5).^(0.4./1.4)-1));

    SEC = (EC_PRES + EC_BD + EC_EVAC + EC_HEAT + EC_FAN)./((n_1_bd - n_1_evac).*0.044)./3600; % kWh/tonne

    recovery_percentageValues(cycle) = recovery_percentage;
    purity_percentageValues(cycle) = purity_percentage;
    productivity_Values(cycle) = productivity;
    SEC_Values(cycle) = SEC;

    process_indicators = [purity_percentageValues; recovery_percentageValues; productivity_Values; SEC_Values];

    if cycle > 11
        for i = 1:4
            for k = 0:10
                if abs(100*(process_indicators(i, cycle-k) - process_indicators(i, cycle-5))/process_indicators(i, cycle-5)) <= 0.04 && ~skipFlag
                    temp_check(k+1) = 1;
                else
                    temp_check(k+1) = 0;
                end
            end
        end
    end

    simTime = toc;
    if parameters.processType == "Resin"
    else
        if simTime > 200
            temp_check(:) = 1;
        end
    end

end




% if round(recovery_percentage,1)>100 || round(purity_percentage,1)>100
% if round(recovery_percentage,1)>100  
%     SEC = 100e6;
%     purity_percentage = 0;
%     recovery_percentage = 0;
%     productivity = 0;
% end




if parameters.processType == "Resin"
    phi_pen = [0, 0];
    phi_pen(1) = 0.80.*((1.*max(0,(95-purity_percentage)))).^2 + (max(0,(0-recovery_percentage))).^2;
    phi_pen(2) = 0.30.*((1.*max(0,(95-purity_percentage))).^2 + (max(0,(0-recovery_percentage))).^2);
    if parameters.OptType == "Const"
        KPIs = [(-productivity+phi_pen(1)) 3600.*50*(SEC.*2.77778e-7+phi_pen(2))];
    else
        KPIs = [ -purity_percentage -recovery_percentage];
    end
else
    phi_pen = [0, 0];
    phi_pen(1) = 0.80.*((1.*max(0,(95-purity_percentage)))).^2 + (max(0,(90-recovery_percentage))).^2;
    phi_pen(2) = 0.30.*((1.*max(0,(95-purity_percentage))).^2 + (max(0,(90-recovery_percentage))).^2);
    if parameters.OptType == "Const"
        KPIs = [(-productivity+phi_pen(1)) 10*(SEC.*2.77778e-7.*3600+phi_pen(2))];
    else
        KPIs = [ -purity_percentage -recovery_percentage];
    end
end

if parameters.OptType == "sampling"
    KPIs = [purity_percentage, recovery_percentage, SEC, productivity,simTime,cycle];
end



catch
    SEC = 100e6;
    purity_percentage = 0;
    recovery_percentage = 0;
    productivity = 0;
    simTime = 1e3;

    if parameters.processType == "Resin"
        if parameters.OptType == "Const"
            phi_pen = [0, 0];
            phi_pen(1) = 0.80.*((1.*max(0,(95-purity_percentage)))).^2 + (max(0,(0-recovery_percentage))).^2;
            phi_pen(2) = 0.30.*((1.*max(0,(95-purity_percentage))).^2 + (max(0,(0-recovery_percentage))).^2);
            KPIs = [(-productivity+phi_pen(1)) 3600.*50*(SEC.*2.77778e-7+phi_pen(2))];
        else
        KPIs = [ -purity_percentage -recovery_percentage];
        end
    else
        if parameters.OptType == "Const"
            phi_pen = [0, 0];
            phi_pen(1) = 0.80.*((1.*max(0,(95-purity_percentage)))).^2 + (max(0,(90-recovery_percentage))).^2;
            phi_pen(2) = 0.30.*((1.*max(0,(95-purity_percentage))).^2 + (max(0,(90-recovery_percentage))).^2);
            KPIs = [(-productivity+phi_pen(1)) 1*(SEC.*2.77778e-7+phi_pen(2))];
        else
        KPIs = [ -purity_percentage -recovery_percentage];
        end
    end


    if parameters.OptType == "sampling"
        KPIs = [purity_percentage, recovery_percentage, SEC, productivity,simTime,cycle];
    end


end

%%
if parameters.outputType == "plot" 
    t_cycle = [t1; t2 + t1(end); t3 + t1(end) + t2(end); t4 + t1(end) + t2(end) + t3(end)]; %Shift the t2 time vector forward in time so it starts immediately after t1 ends
    X_cycle = [X1; X2; X3; X4];
    F_cycleOut = [Fout_ads;Fout_bd;Fout_evac;zeros(length(t4),1)];
    F_cycleIn = [Fin_ads;zeros(length(t2),1);zeros(length(t3),1);Fin_pres];

    if parameters.normPlot
        t0 = parameters.timeRef;
        p0 = parameters.p_H;
        q0 = parameters.refVals(3);
        T0 = parameters.T_feed;
    else
        t0 = 1;
        p0 = 1;
        q0 = 1;
        T0 = 1;
    end
    % P_ads = X1(:,6);
    P1 = X1(:,6);
    P2 = X2(:,6);
    P3 = X3(:,6);
    P4 = X4(:,6);
    t_ads_end  = t1(end);
    t_blo_end  = t_ads_end + t2(end);
    t_evac_end = t_blo_end + t3(end);
    P_cycle = [P1; P2; P3; P4];
    if parameters.cCSTR || ~parameters.pressureDrop
        linestyleVal = '-.';
    else
        linestyleVal = '-';
    end
    P_cycle2 = [ones(length(t1),1).*parameters.P_ads(t1); parameters.P_blo(t2); parameters.P_evac(t3); parameters.P_press(t4)];
    figure(1);
    subplot(4,2,2)
    hold on; xlabel('time [s]'); hold on;
    if parameters.plot0D; 
    plot(t_cycle./t0, P_cycle./1e5,'-', 'Color','#0B0','LineWidth', 2, 'DisplayName','0-D','LineStyle',linestyleVal);% ylabel('P [bar]'); xlabel('time [s]'); hold on;%check unit
    % plot(t_cycle, P_cycle2./1e5,':', 'Color','k','LineWidth', 2, 'DisplayName','0-D');% ylabel('P [bar]'); xlabel('time [s]'); hold on;%check unit
    end
    title('Mean Column Pressure [bar]')
    xline(t_ads_end./t0,  'k--', 'LineWidth', 0.8,'HandleVisibility','off');
    % legend;
    xline(t_blo_end./t0,  'k--', 'LineWidth', 0.8 , 'HandleVisibility','off');
    xline(t_evac_end./t0, 'k--', 'LineWidth', 0.8, 'HandleVisibility','off');
    box on; grid off; set(gca,"LineWidth",2,"FontName","times new roman","FontSize",20)
    xlim([0 t_cycle(end)./t0])


    subplot(4,2,3)
    hold on; xlabel('time [s]'); hold on;
    T = parameters.T_feed;
    if ~parameters.SSLSTA
        [q1_starvals, q2_starvals] = DSL(P_cycle, X_cycle(:,1), X_cycle(:,4), parameters.qsb_1, parameters.qsd_1, parameters.qsb_2, parameters.qsd_2, parameters.bo_1, parameters.do_1, parameters.bo_2, parameters.do_2, parameters.delUb_1, parameters.delUd_1, parameters.delUb_2, parameters.delUd_2);
        [q1_starvalsAds1, q2_starvalsAds] = DSL(X1(:,6), X1(:,1), X1(:,4), parameters.qsb_1, parameters.qsd_1, parameters.qsb_2, parameters.qsd_2, parameters.bo_1, parameters.do_1, parameters.bo_2, parameters.do_2, parameters.delUb_1, parameters.delUd_1, parameters.delUb_2, parameters.delUd_2);
        [q1_starvalsAds, q2_starvalsAds1] = DSL(X1(:,6), ones(length(t_ads),1).*parameters.y1_in, X1(:,4), parameters.qsb_1, parameters.qsd_1, parameters.qsb_2, parameters.qsd_2, parameters.bo_1, parameters.do_1, parameters.bo_2, parameters.do_2, parameters.delUb_1, parameters.delUd_1, parameters.delUb_2, parameters.delUd_2);
    else
        [q1_starvals, q2_starvals] = SSLSTA(P_cycle, X_cycle(:,1), X_cycle(:,4), parameters);
        [q1_starvalsAds1, q2_starvalsAds] = SSLSTA(X1(:,6), X1(:,1), X1(:,4), parameters);
        [q1_starvalsAds, q2_starvalsAds1] = SSLSTA(X1(:,6), ones(length(t_ads),1).*parameters.y1_in, X1(:,4), parameters);
    end

    % [q1_starvalsAds1, q2_starvalsAds] = DSL(parameters.P_ads(t_ads), X1(:,1), T, parameters.qsb_1, parameters.qsd_1, parameters.qsb_2, parameters.qsd_2, parameters.bo_1, parameters.do_1, parameters.bo_2, parameters.do_2, parameters.delUb_1, parameters.delUd_1, parameters.delUb_2, parameters.delUd_2);
    q1_starvals(1:length(t1)) = q1_starvalsAds;
    q2_starvals(1:length(t_ads)) = q2_starvalsAds;
    % plot(t_cycle./t0, q1_starvals,'k--','LineWidth', 2,'DisplayName','q1*');
    if parameters.plot0D; 
        plot(t_cycle./t0, X_cycle(:,2)./1, 'b-','LineWidth', 2 , 'DisplayName','0-D','LineStyle',linestyleVal); hold on;
    end
    % ylabel('q_{CO_{2}} [mol/kg]'); xlabel('time [s]')
    title('Mean Adsorbed amount of CO_{2} [mol/kg]')
    % legend;
    xline(t_ads_end./t0,  'k--', 'LineWidth', 0.8,'HandleVisibility','off');
    xline(t_blo_end./t0,  'k--', 'LineWidth', 0.8 ,'HandleVisibility','off');
    xline(t_evac_end./t0, 'k--', 'LineWidth', 0.8 ,'HandleVisibility','off');
    box on; grid off; set(gca,"LineWidth",2,"FontName","times new roman","FontSize",20)
    xlim([0 t_cycle(end)./t0])

    subplot(4,2,4)
    hold on; xlabel('time [s]'); hold on;
    % plot(t_cycle./t0, q2_starvals, 'k--','LineWidth', 2,'DisplayName','q2*');
    if parameters.plot0D; 
    plot(t_cycle./t0, X_cycle(:,3)./1, 'r-','LineWidth', 2 , 'DisplayName','0-D','LineStyle',linestyleVal); hold on;
    end
    % ylabel('q_{N_{2}} [mol/kg]'); xlabel('time [s]')
    title('Mean Adsorbed amount of N_{2} [mol/kg]')
    xline(t_ads_end./t0,  'k--', 'LineWidth', 0.8, 'HandleVisibility','off');
    xline(t_blo_end./t0,  'k--', 'LineWidth', 0.8 , 'HandleVisibility','off');
    xline(t_evac_end./t0, 'k--', 'LineWidth', 0.8 , 'HandleVisibility','off');
    xlim([0 t_cycle(end)./t0])
    box on; grid off; set(gca,"LineWidth",2,"FontName","times new roman","FontSize",20)
    % legend;

    subplot(4,2,5)
    hold on; xlabel('time [s]'); hold on;
    if parameters.plot0D; 
    plot(t_cycle./t0, X_cycle(:,4)./1, 'm-','LineWidth', 2, 'DisplayName','0-D','LineStyle',linestyleVal);% ylabel('T [K]'); xlabel('time [s]'); hold on;\
    end
    title('Mean Column Temperature [K]')
    xline(t_ads_end./t0,  'k--', 'LineWidth', 0.8,'HandleVisibility','off');
    xline(t_blo_end./t0,  'k--', 'LineWidth', 0.8, 'HandleVisibility','off');
    xline(t_evac_end./t0, 'k--', 'LineWidth', 0.8, 'HandleVisibility','off');
    box on; grid off; set(gca,"LineWidth",2,"FontName","times new roman","FontSize",20)
    set(gca, 'YScale', 'linear')
    xlim([0 t_cycle(end)])
    % legend

    subplot(4,2,6)
    hold on
    if parameters.plot0D; 
    plot(t_cycle./t0, X_cycle(:,5)./1, 'g-','LineWidth', 2, 'DisplayName','0-D','LineStyle',linestyleVal);% ylabel('T_{w} [K]'); xlabel('time [s]'); hold on;
    end
    title('Mean Wall Temperature [K]'); xlabel('time [s]'); hold on;
    xline(t_ads_end./t0,  'k--', 'LineWidth', 0.8,'HandleVisibility','off');
    xline(t_blo_end./t0,  'k--', 'LineWidth', 0.8, 'HandleVisibility','off');
    xline(t_evac_end./t0, 'k--', 'LineWidth', 0.8, 'HandleVisibility','off');
    box on; grid off; set(gca,"LineWidth",2,"FontName","times new roman","FontSize",20)
    set(gca, 'YScale', 'linear')
    xlim([0 t_cycle(end)./t0])
    % legend

    subplot(4,2,1)
    hold on
    if parameters.plot0D; 
    plot(t_cycle./t0, X_cycle(:,1), 'k-','LineWidth', 2, 'DisplayName','0-D','LineStyle',linestyleVal);% ylabel('y_{CO_{2}}')
    end
    title('Mole fraction of CO_{2} in the Outlet [-]'); xlabel('time [s]'); hold on;
    xline(t_ads_end./t0,  'k--', 'LineWidth', 0.8,'HandleVisibility','off');
    xline(t_blo_end./t0,  'k--', 'LineWidth', 0.8, 'HandleVisibility','off');
    xline(t_evac_end./t0, 'k--', 'LineWidth', 0.8, 'HandleVisibility','off');
    box on; grid off; set(gca,"LineWidth",2,"FontName","times new roman","FontSize",20)
    % set(gca, 'YScale', 'l7 og')
    xlim([0 t_cycle(end)./t0])
    ylim([-0.05 1.05])
    % legend

    subplot(4,2,7)
    hold on; xlabel('time [s]'); hold on;
    if parameters.plot0D; 
    plot(t_cycle./t0, F_cycleOut, 'r-','LineWidth', 2, 'DisplayName','0-D','LineStyle',linestyleVal);% ylabel('F_{total}'); xlabel('time [s]'); hold on;
    plot(t_cycle./t0, F_cycleIn, 'b-','LineWidth', 2, 'DisplayName','0-D','LineStyle',linestyleVal);% ylabel('F_{total}'); xlabel('time [s]'); hold on;
    end
    title('Molar Flowrate [mol/s]')
    xline(t_ads_end./t0,  'k--', 'LineWidth', 0.8,'HandleVisibility','off');
    xline(t_blo_end./t0,  'k--', 'LineWidth', 0.8, 'HandleVisibility','off');
    xline(t_evac_end./t0, 'k--', 'LineWidth', 0.8, 'HandleVisibility','off');
    box on; grid off; set(gca,"LineWidth",2,"FontName","times new roman","FontSize",20)
    % set(gca, 'YScale', 'l7 og')
    xlim([0 t_cycle(end)./t0])
    % ylim([-0.05 1.05])
    % legend

    set(gcf,"Position",[100,100,1200 896],"Units","pixels")
    process_indicators = process_indicators';
    CCScycle = cycle

    c1in  = trapz(t_cycle,F_cycleIn);
    c1out  = trapz(t_cycle,F_cycleOut);
    c1in.*parameters.y1_in - trapz(t_cycle,F_cycleOut.*X_cycle(:,1))
    KPIs = process_indicators;
    [purity_percentage, recovery_percentage, SEC, productivity,simTime,cycle]

    figure(834)
    q1vals = X_cycle(:,2);
    q2vals = X_cycle(:,3);
    dq1dt = gradient(q1vals,dt);
    dq2dt = gradient(q2vals,dt);

    Dp = parameters.Dm/parameters.tau; % Effective pore diffusivity [m2/s]

    k01 = (15*parameters.epsilon_p*Dp/(parameters.rp^2))
    cvals = X_cycle((1:length(t1)),1).*X_cycle((1:length(t1)),6)./(8.314.*X_cycle((1:length(t1)),4));
    subplot(1,2,1)
    hold on
    % plot(t1,dq1dt.*q1_starvalsAds1./cvals,'-b','LineWidth',2)
    plot(t_cycle./t0,dq1dt.*1./1.*1./1,'b','LineWidth',2,'LineStyle',linestyleVal)
    box on; grid off; set(gca,"LineWidth",2,"FontName","times new roman","FontSize",20)
    ylabel('\partial{\it{Q}}_1/\partial{\it{t}} [molkg^{-1}s^{-1}]')
    xlabel('time [s]')

    subplot(1,2,2)
    hold on
    box on; grid off; set(gca,"LineWidth",2,"FontName","times new roman","FontSize",20)
    plot(t_cycle./t0,dq2dt.*1./1,'r','LineWidth',2,'LineStyle',linestyleVal)
    ylabel('\partial{\it{Q}}_2/\partial{\it{t}} [molkg^{-1}s^{-1}]')
    xlabel('time [s]')
    ylim([-1e-3 0.1e-3])

    % subplot(2,2,3)
    % hold on
    % plot(q1vals./(q1_starvalsAds(end)),dq1dt./(q1_starvalsAds-q1vals),'-b','LineWidth',2)
    % box on; grid off; set(gca,"LineWidth",2,"FontName","times new roman","FontSize",20)
    %
    % subplot(2,2,4)
    % hold on
    % box on; grid off; set(gca,"LineWidth",2,"FontName","times new roman","FontSize",20)
    % plot(q2vals./(q1_starvalsAds1(end)),dq2dt./(q2_starvalsAds-q2vals),'-r','LineWidth',2)



    if parameters.testBT
        KPIs = [purity_percentage, recovery_percentage, SEC, productivity,simTime,cycle ,X1(1,6)];
    end
elseif parameters.OptType ~= "sampling"
    % Save decision variable to a .txt file
    if exist('rawData','dir')
        fileID = fopen(['rawData',filesep,parameters.fileName,'.txt'],'a+');
        if parameters.processType == "AdsorbentVSA" || parameters.processType == "AdsorbentPVSA"
            fprintf(fileID,'%12.9f %12.9f %12.9f %12.9f %12.9f %12.9f %12.9f %12.9f %12.9f %12.9f %12.9f %12.9f %12.9f %12.9f %12.9f %12.9f %12.9f %12.9f \n', ...
                parameters.p_H, parameters.p_I, parameters.p_L, parameters.F_in, parameters.t_ads, parameters.t_blo, parameters.t_evac,purity_percentage,recovery_percentage,productivity, SEC, theta(1:7));
        elseif parameters.processType == "Resin"
            fprintf(fileID,'%12.9f %12.9f %12.9f %12.9f %12.9f %12.9f %12.9f %12.9f %12.9f %12.9f %12.9f %12.9f %12.9f %12.9f %12.9f %12.9f %12.9f \n', ...
                parameters.p_H, parameters.p_I, parameters.p_L, parameters.F_in, parameters.t_ads, parameters.t_blo, parameters.t_evac,purity_percentage,recovery_percentage,productivity, SEC, parameters.V_column, parameters.v_in, theta(1:3), parameters.qsb_1);
        else
            fprintf(fileID,'%12.9f %12.9f %12.9f %12.9f %12.9f %12.9f %12.9f %12.9f %12.9f %12.9f %12.9f %12.9f %12.9f \n', ...
                parameters.p_H, parameters.p_I, parameters.p_L, parameters.F_in, parameters.t_ads, parameters.t_blo, parameters.t_evac,purity_percentage,recovery_percentage,productivity, SEC, parameters.V_column, parameters.v_in);
        end
        fclose(fileID);
    else
        mkdir rawData
        fileID = fopen(['rawData',filesep,parameters.fileName,'.txt'],'a+');
        if parameters.processType == "AdsorbentVSA" || parameters.processType == "AdsorbentPVSA"
            fprintf(fileID,'%12.9f %12.9f %12.9f %12.9f %12.9f %12.9f %12.9f %12.9f %12.9f %12.9f %12.9f %12.9f %12.9f %12.9f %12.9f %12.9f %12.9f %12.9f \n', ...
                parameters.p_H, parameters.p_I, parameters.p_L, parameters.F_in, parameters.t_ads, parameters.t_blo, parameters.t_evac,purity_percentage,recovery_percentage,productivity, SEC, theta(1:7));
        elseif parameters.processType == "Resin"
            fprintf(fileID,'%12.9f %12.9f %12.9f %12.9f %12.9f %12.9f %12.9f %12.9f %12.9f %12.9f %12.9f %12.9f %12.9f %12.9f %12.9f %12.9f %12.9f \n', ...
                parameters.p_H, parameters.p_I, parameters.p_L, parameters.F_in, parameters.t_ads, parameters.t_blo, parameters.t_evac,purity_percentage,recovery_percentage,productivity, SEC, parameters.V_column, parameters.v_in, theta(1:3), parameters.qsb_1);
        else
            fprintf(fileID,'%12.9f %12.9f %12.9f %12.9f %12.9f %12.9f %12.9f %12.9f %12.9f %12.9f %12.9f %12.9f %12.9f \n', ...
                parameters.p_H, parameters.p_I, parameters.p_L, parameters.F_in, parameters.t_ads, parameters.t_blo, parameters.t_evac,purity_percentage,recovery_percentage,productivity, SEC, parameters.V_column, parameters.v_in);
        end
        fclose(fileID);
    end
end

warning('on','all')

end