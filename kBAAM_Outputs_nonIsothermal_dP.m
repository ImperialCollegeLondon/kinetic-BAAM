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
% and recovery of the heavy product. The blowdown step uses a 2-node (series CSTR)
% spatial discretisation where P, T, and Tw are shared between nodes, and y1/q1/q2
% are resolved independently per node.
%
% Blowdown state vector (9 states):
%   Node 1 (feed end): y1_1, q1_1/qRef, q2_1/qRef
%   Node 2 (product end): y1_2, q1_2/qRef, q2_2/qRef
%   Shared: T/TRef, Tw/TwRef, P/PRef
%
% Initial conditions for blowdown nodes:
%   Node 1: at equilibrium with y1_in at (P_ads_end, T_ads_end)
%   Node 2: same y1, q1, q2 as at the start of adsorption (from previous cycle)
%
% Last modified:
% - 2026-05-01, HA: Add 2-node blowdown using kBAAM_ODEs_nonIsothermal_ND_dP_blo2node
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
%   - kBAAM_ODEs_nonIsothermal_ND_dP_blo2node.m
%   - DSL.m
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [KPIs] = kBAAM_Outputs_nonIsothermal_dP(parameters,varargin)

%% Default parameter settings (assign missing fields to defaults)
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

if ~isfield(parameters,'heating')
    parameters.heating = 0;
end

if ~isfield(parameters,'amine')
    parameters.amine = 0;
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

if ~isfield(parameters,'ResinSens')
    parameters.ResinSens = 0;
end

if ~isfield(parameters,'LDFtest')
    parameters.LDFtest = 0;
end

if ~isfield(parameters,'LDFFactor')
    parameters.LDFFactor = 1;
end

if ~isfield(parameters,'useMassMatrix')
    parameters.useMassMatrix = 0;
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
        parameters.p_H = 101325;
    elseif parameters.processType == "PVSA"
        parameters.v_in = theta(1);
        parameters.p_I = theta(2);
        parameters.t_ads = theta(3);
        parameters.t_blo = theta(4);
        parameters.t_evac = theta(5);

        if length(theta) >6
            parameters.p_L = theta(6);
            parameters.p_H = theta(7);
        else
            parameters.p_L = parameters.p_L;
            parameters.p_H = theta(6);
        end
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
        end
    elseif parameters.processType == "Resin"
        parameters.v_in = theta(6);
        parameters.y1_in = 0.0004;
        parameters.p_I = theta(1);
        parameters.p_L = theta(5);
        parameters.t_ads = theta(2);
        parameters.t_blo = theta(3);
        parameters.qsb_2 = 5.3446e-02.*parameters.qsb_1./2.38;
        parameters.t_evac = theta(4);
        parameters.qsd_1 = 0;
        parameters.qsd_2 = 0;
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
        end
        if ~parameters.fixResins
            parameters.LDF = 0.132.*exp(-2.076.*parameters.qsb_1); % correlation including Lewatit
            parameters.bo_1 = 1.06e-16.*parameters.qsb_1; % correlation including Lewatit
        end
    end
end

%% Process parameters and step time vectors
Rg = 8.3145; % universal gas constant [J/mol/K]

if parameters.processType == "Resin" || parameters.processType == "ResinSens"
    dt = 0.04;
else
    dt = 0.05; % [s]
end

if parameters.heating
    parameters.heatPowerDensity = 5e3; % [W/m2]
    parameters.h_in = 8.6*3;
    parameters.Theat = parameters.Theat+6;
else
    parameters.Theat = 0;
end

if parameters.ResinSens
    parameters.qsb_1 = parameters.qsb_1.*theta(7);
    parameters.qsb_2 = parameters.qsb_2.*theta(8);
    parameters.LDF = parameters.LDF.*theta(9);
    parameters.delUb_1 = parameters.delUb_1.*theta(10);
    parameters.rho_s = parameters.rho_s.*theta(11);
    if ~parameters.fixResins
        parameters.LDF = 0.132.*exp(-2.076.*parameters.qsb_1); % correlation including Lewatit
        parameters.bo_1 = 1.06e-16.*parameters.qsb_1; % correlation including Lewatit
    end
end

t_ads   = 0:dt:parameters.t_ads; % time vector for adsorption step [s]
t_blo   = 0:dt:parameters.t_blo; % time vector for blowdown step [s]
t_evac  = 0:dt:parameters.t_evac; % time vector for evacuation step [s]
t_press = 0:dt:parameters.t_press; % time vector for pressurization step [s]

%% Column geometry, inlet conditions, and reference values
parameters.r_in = (parameters.V_column./(parameters.Lbyr.*pi)).^(1./3); % inner radius of column [m]
parameters.r_out = parameters.r_in + 0.0175; % outer radius of column [m]
parameters.v_in = parameters.v_in./parameters.e_bed; % superficial to interstitial
parameters.L = parameters.Lbyr.*parameters.r_in; % length of column [m]
parameters.A_in = parameters.r_in.^2.*pi; % inner cross sectional area of column [m2]
parameters.volFlowin = parameters.v_in.*parameters.A_in.*parameters.e_bed; % inlet volumetric flowrate [m3/s]
parameters.deltaP = 0;
parameters.F_in = parameters.volFlowin.*(parameters.p_H+parameters.deltaP)./(Rg.*parameters.T_feed); % inlet molar flowrate (ideal gas) [mol/s]

parameters.refVals = [1,(parameters.qsb_1+parameters.qsd_1), (parameters.qsb_1+parameters.qsd_1), parameters.T_feed, parameters.T_feed, parameters.p_H]; % vector of reference values of state variables for non-dimensionalization

% Reference values for non-dimensionalization
parameters.timeRef = parameters.V_column./parameters.volFlowin;  % reference time [s]
parameters.TRef = parameters.T_feed; % reference temperature [K]
parameters.TwRef = parameters.T_feed; % reference temperature [K]
parameters.PRef = parameters.p_H; % reference pressure [Pa]


% pressure profiles and derivatives for overall material balance
parameters.P_ads = @(t)parameters.p_H+parameters.deltaP; % pressure profile for adsorption (constant)
parameters.P_initH = parameters.P_ads(parameters.t_ads); % pressure at the end of ads
parameters.P_blo = @(t)parameters.p_I+(parameters.P_initH-parameters.p_I)*exp(-parameters.lambda*t); % pressure profile for blowdown
parameters.P_initL = parameters.P_blo(parameters.t_blo); % pressure at the end of blowdown
parameters.P_evac = @(t)parameters.p_L+(parameters.p_I-parameters.p_L)*exp(-parameters.lambda*t); % pressure profile for evacuation
parameters.P_initR = parameters.P_evac(parameters.t_evac); % pressure at the end of evacuation
parameters.P_press = @(t)parameters.P_initH+(parameters.P_initR-parameters.P_initH)*exp(-parameters.lambda*t); % pressure profile for pressurization
parameters.dPdt_blo = @(t)-parameters.lambda*(parameters.P_initH-parameters.p_I)*exp(-parameters.lambda*t); % time derivative of pressure profile for blowdown
parameters.dPdt_evac =  @(t)-parameters.lambda*(parameters.P_initL-parameters.p_L)*exp(-parameters.lambda*t); % time derivative of pressure profile for evacuation
parameters.dPdt_press = @(t)-parameters.lambda*(parameters.P_initR-parameters.P_initH)*exp(-parameters.lambda*t); % time derivative of pressure profile for pressurization


%% Initial condition for matrix of solution states
y1Init = 0.99; % initial mole fraction of component 1 in bed [-]
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

%% Pre-compute step-invariant constants (avoid recomputation in ODE)
parameters.isIsothermal = (parameters.modelType == "isothermal");
parameters.isResin = (parameters.processType == "Resin" || parameters.processType == "ResinSens");
parameters.Ab = (1 - parameters.e_bed) / parameters.e_bed;
parameters.darcyK = (4/150/1.72e-5) * (parameters.e_bed/(1-parameters.e_bed))^2 * parameters.rp^2;
parameters.tRef_qRef = parameters.timeRef / parameters.qRef;
parameters.Ab_rhos_qRef_tRef = parameters.Ab * parameters.rho_s * parameters.qRef / parameters.timeRef;
parameters.inv_tRef = 1 / parameters.timeRef;
parameters.PRef_tRef = parameters.PRef / parameters.timeRef;
parameters.R_TRef = 8.3145 * parameters.TRef;
parameters.cpg_eV = parameters.cp_g / (parameters.V_column * parameters.e_bed);
parameters.two_hin_rin_e = 2 * parameters.h_in / (parameters.r_in * parameters.e_bed);
parameters.wall_coeff1 = 2 * parameters.h_in * parameters.r_in / (parameters.r_out^2 - parameters.r_in^2);
parameters.wall_coeff2 = 2 * parameters.h_out * parameters.r_out / (parameters.r_out^2 - parameters.r_in^2);
parameters.wall_prefactor = parameters.timeRef / parameters.TwRef / (parameters.rho_w * parameters.cp_w);
% Flow / Darcy composite constants (used in every ODE call)
parameters.darcy_PRef_2overL = (2 / parameters.L) * parameters.darcyK * parameters.PRef; % (2/L)*darcyK*PRef [m/s/Pa * Pa = m/s]
parameters.Fin_ads_prefactor  = parameters.volFlowin * parameters.PRef / (8.3145 * parameters.TRef); % volFlowin*PRef/(R*TRef) [mol/s]
parameters.Fout_prefactor     = parameters.A_in * parameters.e_bed / 8.3145;                         % A*e/R [mol*K/J = mol/(Pa·m)]

if parameters.pressureDrop
    odeFunc = @kBAAM_ODEs_nonIsothermal_ND_dP;
else
    odeFunc = @kBAAM_ODEs_nonIsothermal_ND_nodP;
end

cycle = 0; % initialize cycle number
warning('off','all')
tic

if parameters.testBT
    max_no_Cycles = 1;
else
    max_no_Cycles = 200; % Maximum number of cycles to run to test CSS
end
parameters.y1_LPP = 0.15; % initialize LP composition
options = odeset('RelTol', 1e-5, 'AbsTol', 1e-5, 'MaxOrder', 2);

% [delH1,~] = computeDSLHeatUnary(linspace(0,10e5,1000), 1, parameters.T_feed, parameters.PRef, parameters.TRef, parameters);
% [~,delH2] = computeDSLHeatUnary(linspace(0,10e5,1000), 0, parameters.T_feed, parameters.PRef, parameters.TRef, parameters);
% 
% [q1vals, ~] = DSL(linspace(0,20e5,1000), 1, parameters.T_feed, parameters.qsb_1, parameters.qsd_1, parameters.qsb_2, parameters.qsd_2, parameters.bo_1, parameters.do_1, parameters.bo_2, parameters.do_2, parameters.delUb_1, parameters.delUd_1, parameters.delUb_2, parameters.delUd_2); % initial adsorbed amounts in bed [mol/kg]
% [~, q2vals] = DSL(linspace(0,20e5,1000), 0, parameters.T_feed, parameters.qsb_1, parameters.qsd_1, parameters.qsb_2, parameters.qsd_2, parameters.bo_1, parameters.do_1, parameters.bo_2, parameters.do_2, parameters.delUb_1, parameters.delUd_1, parameters.delUb_2, parameters.delUd_2); % initial adsorbed amounts in bed [mol/kg]
% 
% parameters.heat1 = [q1vals' delH1'];
% parameters.heat2 = [q2vals' delH2'];

try
    while  cycle < max_no_Cycles && mean(temp_check) < 1  && parameters.p_H > parameters.p_I && parameters.p_I > parameters.p_L

        cycle = cycle+1;

        parameters.loadingFraction = 1;
        parameters.P_press = @(t)parameters.P_initH+(parameters.P_initR-parameters.P_initH)*exp(-parameters.lambda*t); % pressure profile for pressurization
        parameters.dPdt_press = @(t)-parameters.lambda*(parameters.P_initR-parameters.P_initH)*exp(-parameters.lambda*t); % time derivative of pressure profile for pressurization

        [t4, X4] = ode15s(@(t,X) odeFunc(t,X,parameters,'pres'), t_press./(parameters.timeRef), X0, options);%t1 is the time point at which the solution is evaluated, X1 is the solution states for adsorption step
        t4 = t4.*parameters.timeRef;
        X4(X4<0) = 0;
        X4(X4(:,1)>1,1) = 1;
        if parameters.testBT && ~parameters.testEvac
            parameters.P_initH = parameters.p_H+150./4.*1./parameters.rp.^2.*((1-parameters.e_bed)./parameters.e_bed).^2.*1.72e-5.*parameters.v_in.*parameters.L./1;
            y1Init = 0.000001; % initial mole fraction of component 1 in bed [-]
            if ~parameters.SSLSTA
                [q1Init, q2Init] = DSL(parameters.P_initH, y1Init, parameters.T_feed, parameters.qsb_1, parameters.qsd_1, parameters.qsb_2, parameters.qsd_2, parameters.bo_1, parameters.do_1, parameters.bo_2, parameters.do_2, parameters.delUb_1, parameters.delUd_1, parameters.delUb_2, parameters.delUd_2); % initial adsorbed amounts in bed [mol/kg]
            else
                [q1Init, q2Init] = SSLSTA(parameters.P_initH, y1Init, parameters.T_feed, parameters); % initial adsorbed amounts in bed [mol/kg]
            end
            X0 = [y1Init; q1Init; q2Init; parameters.T_feed; parameters.T_feed;parameters.P_initH]./parameters.refVals'; % dimensionless vector of initial states
        else
            X0 = X4(end,:)'; %This sets up the initial condition for the next step, by taking the final state from the previous step.
        end
        X4 = X4.*parameters.refVals;
        parameters.q1init = X0(2).*parameters.qRef ;
        parameters.y1init = X0(1);
        parameters.q2init = X0(3).*parameters.qRef;

        [q1max, ~] = DSL(X4(end,6), parameters.y1_in, X4(end,4), parameters.qsb_1, parameters.qsd_1, parameters.qsb_2, parameters.qsd_2, parameters.bo_1, parameters.do_1, parameters.bo_2, parameters.do_2, parameters.delUb_1, parameters.delUd_1, parameters.delUb_2, parameters.delUd_2); % initial adsorbed amounts in bed [mol/kg]
        [q10, ~] = DSL(X4(end,6), X4(end,1), X4(end,4), parameters.qsb_1, parameters.qsd_1, parameters.qsb_2, parameters.qsd_2, parameters.bo_1, parameters.do_1, parameters.bo_2, parameters.do_2, parameters.delUb_1, parameters.delUd_1, parameters.delUb_2, parameters.delUd_2); % initial adsorbed amounts in bed [mol/kg]
        parameters.loadingFraction = (X4(end,2)-X4(1,2))./((q1max-q10));


        [t1, X1] = ode15s(@(t,X) odeFunc(t,X,parameters,'ads'), t_ads./(parameters.timeRef), X0, options); %t1 is the time point at which the solution is evaluated, X1 is the solution states for adsorption step
        t1 = t1.*parameters.timeRef;
        X1(X1<0) = 0;
        X1(X1(:,1)>1,1) = 1;
        X0 = X1(end,:)';
        X1 = X1.*parameters.refVals;
        if parameters.processType == "Resin" || parameters.processType == "ResinSens"
            if t1(end) < 0.95*parameters.t_ads
                parameters.t_ads = t1(end); % Update evacuation time to length of simulation instead of discarding. Integration seems to fail when flowrate in close to 0 and mole fraction close to 1.
                t_ads  = 0:dt:parameters.t_ads;
            end
        end

        % if parameters.pressureDrop
        %     parameters.P_initH = X1(end,6);
        %     parameters.P_blo = @(t)parameters.p_I+(parameters.P_initH-parameters.p_I)*exp(-parameters.lambda*t); % pressure profile for blowdown
        % end

        [q1max, ~] = DSL(X1(end,6), parameters.y1_in, X1(end,4), parameters.qsb_1, parameters.qsd_1, parameters.qsb_2, parameters.qsd_2, parameters.bo_1, parameters.do_1, parameters.bo_2, parameters.do_2, parameters.delUb_1, parameters.delUd_1, parameters.delUb_2, parameters.delUd_2); % initial adsorbed amounts in bed [mol/kg]
        [q10, ~] = DSL(X4(end,6), X4(end,1), X4(end,4), parameters.qsb_1, parameters.qsd_1, parameters.qsb_2, parameters.qsd_2, parameters.bo_1, parameters.do_1, parameters.bo_2, parameters.do_2, parameters.delUb_1, parameters.delUd_1, parameters.delUb_2, parameters.delUd_2); % initial adsorbed amounts in bed [mol/kg]
        parameters.loadingFraction = (X1(end,2)-q10)./((q1max-X1(1,2)));
        parameters.q1init = X1(end,2);
        parameters.q2init = X1(end,3);
        parameters.y1init = X1(end,1) ;
        % Adsorption inlet flow (Ergun-corrected pressure at bed exit) [mol/s]
        Fin_ads = parameters.volFlowin.*(2.*X1(:,6)-parameters.p_H)./(Rg.*parameters.T_feed);
        if parameters.pressureDrop
            % Adsorption outlet flowrate from overall material balance [mol/s]:
            v_outA = (2./parameters.L) .* parameters.darcyK .* (X1(:,6) - parameters.p_H);
            Fout_ads  = parameters.p_H .* parameters.A_in .* parameters.e_bed ./ (Rg .* X1(:,4)) .* v_outA;
            Fout_ads(Fout_ads<0) = 0;
        else
            % Adsorption outlet flowrate from overall material balance [mol/s]:
            %   Fout = Fin - adsorption sink - gas compression + thermal expansion
            Fout_ads = Fin_ads ...
                - (1 - parameters.e_bed) .* parameters.V_column .* parameters.rho_s .* (gradient(X1(:,2),dt) + gradient(X1(:,3),dt)) ...  % adsorption sink
                - (parameters.e_bed ./ (Rg*X1(:,4))) .* gradient(X1(:,6),dt) .* parameters.V_column ...                                  % gas compression term
                + (parameters.e_bed .* X1(:,6) ./ (Rg*(X1(:,4)).^2)) .* gradient(X1(:,4),dt) .* parameters.V_column ;                   % thermal expansion term
            Fout_ads(Fout_ads<0) = 0;
        end
        % CO2 moles and average mole fraction in ads-step effluent (used for LPP pressurisation)
        F_1_out_ads = Fout_ads.*X1(:,1);
        mol_1_out_ads = trapz(t1,F_1_out_ads); moltot_out_ads = trapz(t1,Fout_ads);
        parameters.y1_LPP = mol_1_out_ads./moltot_out_ads; % avg CO2 mole fraction in LPP gas [-]

        f_w = max(0.05, min(0.95, parameters.loadingFraction));  % same clamp as ODE

        %% BLOWDOWN
        if parameters.pressureDrop*(f_w<0.95)
            %% ---- Blowdown: 2-node integration ----
            options2 = odeset('RelTol', 1e-5, 'AbsTol', 1e-5, 'MaxOrder', 2);
            % refVals for the 10-state blowdown vector:
            %   [1, qRef, qRef, 1, qRef, qRef, TRef, TwRef, PRef, PRef]
            refVals_blo = [1, parameters.qRef, parameters.qRef, ...    % node 1: y1_1, q1_1, q2_1
                1, parameters.qRef, parameters.qRef, ...    % node 2: y1_2, q1_2, q2_2
                parameters.TRef, parameters.TwRef, ...      % shared T, Tw
                parameters.PRef, parameters.PRef];          % P1, P2

            % Node 1 initial condition: equilibrium with y1_in at end-of-ads (P, T)
            P_ads_end = parameters.p_H;  % [Pa]
            T_ads_end = X1(end,4);  % [K]
            if ~parameters.SSLSTA
                [q1_n1, q2_n1] = DSL(P_ads_end, parameters.y1_in, T_ads_end, parameters.qsb_1, parameters.qsd_1, parameters.qsb_2, parameters.qsd_2, parameters.bo_1, parameters.do_1, parameters.bo_2, parameters.do_2, parameters.delUb_1, parameters.delUd_1, parameters.delUb_2, parameters.delUd_2);
            else
                [q1_n1, q2_n1] = SSLSTA(P_ads_end, parameters.y1_in, T_ads_end, parameters);
            end

            % Node 2 initial condition: y1 at end of adsorption; q1, q2 at start of adsorption
            % (product end sees end-of-ads composition but start-of-ads loading)
            y1_n2 = X1(end,1); % [mol frac] — end-of-ads mole fraction
            q1_n2 = X1(1,2);   % [mol/kg]  — start-of-ads loading
            q2_n2 = X1(1,3);   % [mol/kg]  — start-of-ads loading

            % Both nodes start at the same pressure (end-of-ads column pressure)
            X0_blo = [parameters.y1_in; q1_n1; q2_n1; ...   % node 1
                y1_n2;            q1_n2; q2_n2; ...   % node 2
                T_ads_end; X1(end,5); ...              % shared T, Tw
                P_ads_end; P_ads_end] ...              % P1, P2 (same at t=0)
                ./ refVals_blo';

            [t2, X2_10] = ode15s(@(t,X) kBAAM_ODEs_nonIsothermal_ND_dP_blo2node(t,X,parameters,'blo'), t_blo./(parameters.timeRef), X0_blo, options2);
            t2 = t2 .* parameters.timeRef;
            X2_10(X2_10 < 0) = 0;
            X2_10(X2_10(:,1) > 1, 1) = 1;
            X2_10(X2_10(:,4) > 1, 4) = 1;

            % Scale back to dimensional
            X2_10dim = X2_10 .* refVals_blo;

            % Collapse to 6-state for evac IC: loadingFraction-weighted average,
            % P = volume-weighted average of P1 and P2
            f_w = max(0.05, min(0.9, parameters.loadingFraction));  % same clamp as ODE
            y1_blo_end = f_w.*X2_10dim(end,1) + (1-f_w).*X2_10dim(end,4);
            q1_blo_end = f_w.*X2_10dim(end,2) + (1-f_w).*X2_10dim(end,5);
            q2_blo_end = f_w.*X2_10dim(end,3) + (1-f_w).*X2_10dim(end,6);
            T_blo_end  = X2_10dim(end,7);
            Tw_blo_end = X2_10dim(end,8);
            P_blo_end  = f_w.*X2_10dim(end,9) + (1-f_w).*X2_10dim(end,10);  % vol-weighted avg P
            X0 = [y1_blo_end; q1_blo_end; q2_blo_end; T_blo_end; Tw_blo_end; P_blo_end] ./ parameters.refVals';

            % Build a pseudo 6-column X2 (time × 6) for KPI/energy calculations below.
            % Columns: [y1_avg, q1_avg, q2_avg, T, Tw, P_avg]
            P_avg_blo = f_w.*X2_10dim(:,9) + (1-f_w).*X2_10dim(:,10);
            X2 = [ f_w.*X2_10dim(:,1) + (1-f_w).*X2_10dim(:,4), ...
                f_w.*X2_10dim(:,2) + (1-f_w).*X2_10dim(:,5), ...
                f_w.*X2_10dim(:,3) + (1-f_w).*X2_10dim(:,6), ...
                X2_10dim(:,7), X2_10dim(:,8), P_avg_blo];

            % Outlet composition during blowdown is node-2 (product-end) composition
            y1_bd_out = X2_10dim(:,4);
        else
            [t2, X2] = ode15s(@(t,X) odeFunc(t,X,parameters,'blo'), t_blo./(parameters.timeRef), X0, options); %t1 is the time point at which the solution is evaluated, X1 is the solution states for adsorption step
            t2 = t2.*parameters.timeRef;
            X2(X2<0) = 0;
            X2(X2(:,1)>1,1) = 1;
            X0 = X2(end,:)';
            X2 = X2.*parameters.refVals;
            f_w = 1;
            % Outlet composition during blowdown is node-2 (product-end) composition
            y1_bd_out = X2(:,1);
        end

        [t3, X3] = ode15s(@(t,X) odeFunc(t,X,parameters,'evac'), t_evac./(parameters.timeRef), X0, options);%t1 is the time point at which the solution is evaluated, X1 is the solution states for adsorption step
        t3 = t3.*parameters.timeRef;
        X3(X3<0) = 0;
        X3(X3(:,1)>1,1) = 1;
        X0 = X3(end,:)'; %This sets up the initial condition for the next step, by taking the final state from the previous step.%end means the last row of X1, containing all the state variables at %the t final , : takes all the columns (for all types of the state%variables
        X3 = X3.*parameters.refVals;

        if parameters.processType == "Resin" || parameters.processType == "ResinSens"
            if t3(end) < 0.95*parameters.t_evac
                parameters.t_evac = t3(end); % Update evacuation time to length of simulation instead of discarding. Integration seems to fail when flowrate in close to 0 and mole fraction close to 1.
                t_evac  = 0:dt:parameters.t_evac;
            end
        end

        %% Mole inventories at end of each step (gas-phase + adsorbed-phase) [mol]
        n_1_ads = (X1(end,1).*X1(end,6) * parameters.V_column * parameters.e_bed / (Rg * X1(end,4))) + X1(end,2)* parameters.V_column * (1-parameters.e_bed).*parameters.rho_s ;
        f_blo = f_w;  % same clamp as ODE
        n_1_bd = (X3(1,1) .*X3(1,6) * parameters.V_column * parameters.e_bed / (Rg * X3(1,4))) + X3(1,2)* parameters.V_column * (1-parameters.e_bed).*parameters.rho_s ;
        n_1_evac = (X3(end,1) .*X3(end,6) * parameters.V_column * parameters.e_bed / (Rg * X3(end,4))) + X3(end,2)* parameters.V_column * (1-parameters.e_bed).*parameters.rho_s ;
        n_2_bd = ((1-X3(1,1)).*X3(1,6) * parameters.V_column * parameters.e_bed / (Rg * X3(1,4))) + X3(1,3)* parameters.V_column * (1-parameters.e_bed).*parameters.rho_s ;
        n_2_evac = ((1-X3(end,1)).*X3(end,6) * parameters.V_column * parameters.e_bed / (Rg * X3(end,4))) + X3(end,3)* parameters.V_column * (1-parameters.e_bed).*parameters.rho_s ;
        n_1_pres = (X4(end,1).*X4(end,6) * parameters.V_column * parameters.e_bed / (Rg * X4(end,4))) + X4(end,2)* parameters.V_column * (1-parameters.e_bed).*parameters.rho_s ;
        n_1_presInit = (X4(1,1).*X4(1,6) * parameters.V_column * parameters.e_bed / (Rg * X4(1,4))) + X4(1,2)* parameters.V_column * (1-parameters.e_bed).*parameters.rho_s ;

        cycle_time = (parameters.t_ads + parameters.t_blo + parameters.t_evac + parameters.t_press);

        productivity = (n_1_bd - n_1_evac) /(parameters.V_column.*cycle_time);

        %% Energy Calculation
        % Step flowrates from overall material balance (Fin=0 for blo/evac; Fout=0 for pres)
        if parameters.pressureDrop.*(f_w<0.95)
            % Blowdown outlet flowrate [mol/s]: Fin=0 (inlet valve closed)
            % Blowdown outlet flow exits from node-2 (product end); use P2 directly
            % v_out distance = node-2 centre to outlet = L2/2, so factor = 2/L2
            L2_blo   = (1 - f_blo) * parameters.L;
            v_outB   = (2./L2_blo) .* parameters.darcyK .* (X2_10dim(:,10) - parameters.P_blo(t2));
            Fout_bd  = X2_10dim(:,10) .* parameters.A_in .* parameters.e_bed ./ (Rg .* X2_10dim(:,7)) .* v_outB;
            Fout_bd(Fout_bd<0) = 0; % enforce non-negative
            % Evacuation outlet flowrate [mol/s]: Fin=0; dP/dt compression term suppressed (pump-dominated)
            v_outE = (2./parameters.L) .* parameters.darcyK .* (X3(:,6) - parameters.P_evac(t3));
            Fout_evac  = X3(:,6) .* parameters.A_in .* parameters.e_bed ./ (Rg .* X3(:,4)) .* v_outE;
            Fout_evac(Fout_evac<0) = 0; % enforce non-negative
            % Pressurisation inlet flowrate [mol/s]: Fout=0 (outlet valve closed)
            v_outP = -(2./parameters.L) .* parameters.darcyK .* (X4(:,6) - parameters.P_press(t4));
            Fin_pres  = parameters.P_press(t4) .* parameters.A_in .* parameters.e_bed ./ (Rg .* parameters.T_feed) .* v_outP;
        else
            %% Energy Calculation
            % Step flowrates from overall material balance (Fin=0 for blo/evac; Fout=0 for pres)

            % Blowdown outlet flowrate [mol/s]: Fin=0 (inlet valve closed)
            Fout_bd = 0 - (1 - parameters.e_bed) .* parameters.V_column * parameters.rho_s .* (gradient(X2(:,2),dt) + gradient(X2(:,3),dt)) ...
                - (parameters.e_bed ./ (Rg.*X2(:,4))) .* gradient(X2(:,6),dt) .* parameters.V_column ...
                + (parameters.e_bed .* X2(:,6)./ (Rg.*X2(:,4).^2)) .* gradient(X2(:,4),dt) .* parameters.V_column ;
            % v_outB = (2./parameters.L) .* parameters.darcyK .* (X2(:,6) - parameters.P_blo(t2));
            % Fout_bd  = X2(:,6) .* parameters.A_in .* parameters.e_bed ./ (Rg .* X2(:,4)) .* v_outB;
            Fout_bd(Fout_bd<0) = 0; % enforce non-negative
            % Evacuation outlet flowrate [mol/s]: Fin=0; dP/dt compression term suppressed (pump-dominated)
            Fout_evac = 0 - (1 - parameters.e_bed) .* parameters.V_column * parameters.rho_s .* (gradient(X3(:,2),dt) + gradient(X3(:,3),dt)) ...
                - (parameters.e_bed ./ (Rg.*X3(:,4))) .* gradient(X3(:,6),dt) .* parameters.V_column ...
                + (parameters.e_bed .* X3(:,6)./ (Rg.*X3(:,4).^2)) .* gradient(X3(:,4),dt) .* parameters.V_column ;
            % Pressurisation inlet flowrate [mol/s]: Fout=0 (outlet valve closed)
            Fin_pres =  (1 - parameters.e_bed) .* parameters.V_column * parameters.rho_s .* (gradient(X4(:,2),dt) + gradient(X4(:,3),dt)) ...
                + (parameters.e_bed ./ (Rg.*X4(:,4))) .* gradient(X4(:,6),dt) .* parameters.V_column ...
                - (parameters.e_bed .* X4(:,6)./ (Rg.*X4(:,4).^2)) .* gradient(X4(:,4),dt) .* parameters.V_column ;

        end



        if parameters.pressType == "LPP"
            mole_LP_recycle = n_1_pres-n_1_presInit;
            recovery_percentage = 100 * (n_1_bd - n_1_evac) /  ((n_1_ads - n_1_evac + mol_1_out_ads-max(0,mole_LP_recycle)));
        else
            recovery_percentage = 100 * (n_1_bd - n_1_evac) / ((n_1_ads - n_1_evac + mol_1_out_ads));
            % recovery_percentage = 100 * (trapz(t3,Fout_evac.*X3(:,1))) / ((trapz(t3,Fout_evac.*X3(:,1)) +trapz(t2,Fout_bd.*y1_bd_out) + mol_1_out_ads));
            % recovery_percentage = 100 * (n_1_bd - n_1_evac) / ((trapz(t1,Fin_ads.*parameters.y1_in) + trapz(t4,Fin_pres.*parameters.y1_in)));
        end
        purity_percentage = 100 * (n_1_bd - n_1_evac) / (n_1_bd - n_1_evac + max(0,n_2_bd - n_2_evac));

        eta_bd = 0.8.*(19.55.*parameters.P_blo(t2).*1e-5./(1+19.55.*parameters.P_blo(t2).*1e-5));
        eta_evac = 0.8.*(19.55.*parameters.P_evac(t3).*1e-5./(1+19.55.*parameters.P_evac(t3).*1e-5));
        eta_press = 0.72;
        eta_ads = 0.72;

        P_atm = 101325; % Atmospheric pressure [Pa]

        if parameters.heating
            heatFlag = X3(:,4) < parameters.Theat;
            Qheat = trapz(t3,parameters.heatPowerDensity.*(parameters.Theat-X3(:,4))./(parameters.Theat-parameters.T_feed).*parameters.r_out.*2.*parameters.L.*heatFlag); % external heat flux if heating is used
            EC_HEAT = Qheat;
        else
            EC_HEAT = 0;
        end

        % Adiabatic compression work [J]: W = (F*R*T/eta) * (gamma/(gamma-1)) * ((P_out/P_in)^((gamma-1)/gamma) - 1)
        % BD/EVAC: column gas temperature used (gas is already inside column)
        % PRES/FAN: T_feed used (compressor/fan draws in ambient feed gas)
        EC_BD   = trapz(t2,1./eta_bd    .*Fout_bd  .*Rg.*X2(:,4).*(1.4./0.4).*((P_atm./min(P_atm,parameters.P_blo(t2))).^(0.4./1.4)-1));  % vacuum pump work, BD step [J]
        EC_EVAC = trapz(t3,1./eta_evac  .*Fout_evac.*Rg.*X3(:,4).*(1.4./0.4).*((P_atm./min(P_atm,parameters.P_evac(t3))).^(0.4./1.4)-1)); % vacuum pump work, evac step [J]
        EC_PRES = trapz(t4,1./eta_press .*Fin_pres .*Rg.*parameters.T_feed.*(1.4./0.4).*((max(P_atm,parameters.P_press(t4))./P_atm).^(0.4./1.4)-1)); % compressor work, pres step [J]
        EC_FAN  = trapz(t1, 1./eta_ads  .*Fin_ads  .*Rg.*parameters.T_feed.*(1.4./0.4).*((max(P_atm,(2.*X1(:,6)-parameters.p_H))./P_atm).^(0.4./1.4)-1));  % fan work, ads step [J]

        SEC = (EC_PRES + EC_BD + EC_EVAC + EC_HEAT + EC_FAN)./((n_1_bd - n_1_evac).*0.04401)./3600; % specific energy consumption [kWh/tonne CO2]

        recovery_percentageValues(cycle) = recovery_percentage;
        purity_percentageValues(cycle) = purity_percentage;
        productivity_Values(cycle) = productivity;
        SEC_Values(cycle) = SEC;


        t_cycle = [t1; t2 + t1(end); t3 + t1(end) + t2(end); t4 + t1(end) + t2(end) + t3(end)]; %Shift the t2 time vector forward in time so it starts immediately after t1 ends
        X_cycle = [X1; X2; X3; X4];
        F_cycleOut = [Fout_ads;Fout_bd;Fout_evac;zeros(length(t4),1)];
        F_cycleIn = [Fin_ads;zeros(length(t2),1);zeros(length(t3),1);Fin_pres];

        mol1in  = trapz(t_cycle,F_cycleIn.*1);
        % mol1out  = trapz(t_cycle,F_cycleOut.*X_cycle(:,1));
        mol1out  = trapz(t_cycle,F_cycleOut.*1);
        MBerror = mol1in - mol1out;
        MBerrorVals(cycle) = MBerror;

        process_indicators = [purity_percentageValues; recovery_percentageValues; productivity_Values; SEC_Values];

        if parameters.outputType == "plot"
            process_indicators'
        end

        if cycle > 6
            for i = 1:4
                for k = 0:5
                    if abs(100*(process_indicators(i, cycle-k) - process_indicators(i, cycle-5))/process_indicators(i, cycle-5)) <= 0.02
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

    if round(recovery_percentage,1) < 0  ||  recovery_percentage > 100 ||  purity_percentage > 100 || round(purity_percentage,1) < 0
        SEC = 100e6;
        purity_percentage = 0;
        recovery_percentage = 0;
        productivity = 0;
    end

    if parameters.processType == "Resin" || parameters.amine
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

    if parameters.processType == "Resin" || parameters.amine
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
    X2(:,1) = y1_bd_out;
    X_cycle = [X1; X2; X3; X4];
    F_cycleOut = [Fout_ads;Fout_bd;Fout_evac;zeros(length(t4),1)];
    F_cycleIn = [Fin_ads;zeros(length(t2),1);zeros(length(t3),1);Fin_pres];

    if parameters.normPlot
        t0 = parameters.timeRef;
    else
        t0 = 1;
    end
    P1 = X1(:,6);
    P2 = X2(:,6);  % shared P (column 6 of collapsed X2)
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
    if parameters.plot0D
        plot(t_cycle./t0, P_cycle./1e5,'-', 'Color','#0B0','LineWidth', 1.5, 'DisplayName','0-D','LineStyle',linestyleVal);% ylabel('P [bar]'); xlabel('time [s]'); hold on;%check unit
        % plot(t_cycle./t0, P_cycle2./1e5,'--', 'Color','#0B0','LineWidth', 1.5, 'DisplayName','0-D','LineStyle','--');% ylabel('P [bar]'); xlabel('time [s]'); hold on;%check unit
    end
    title('Mean Column Pressure [bar]')
    xline(t_ads_end./t0,  'k--', 'LineWidth', 0.8,'HandleVisibility','off');
    xline(t_blo_end./t0,  'k--', 'LineWidth', 0.8 , 'HandleVisibility','off');
    xline(t_evac_end./t0, 'k--', 'LineWidth', 0.8, 'HandleVisibility','off');
    box on; grid off; set(gca,"LineWidth",2,"FontName","times new roman","FontSize",20)
    xlim([0 t_cycle(end)./t0])


    subplot(4,2,3)
    hold on; xlabel('time [s]'); hold on;
    if ~parameters.SSLSTA
        [q1_starvals, q2_starvals] = DSL(P_cycle, X_cycle(:,1), X_cycle(:,4), parameters.qsb_1, parameters.qsd_1, parameters.qsb_2, parameters.qsd_2, parameters.bo_1, parameters.do_1, parameters.bo_2, parameters.do_2, parameters.delUb_1, parameters.delUd_1, parameters.delUb_2, parameters.delUd_2);
        [q1_starvalsAds1, q2_starvalsAds] = DSL(X1(:,6), X1(:,1), X1(:,4), parameters.qsb_1, parameters.qsd_1, parameters.qsb_2, parameters.qsd_2, parameters.bo_1, parameters.do_1, parameters.bo_2, parameters.do_2, parameters.delUb_1, parameters.delUd_1, parameters.delUb_2, parameters.delUd_2);
        [q1_starvalsAds, q2_starvalsAds1] = DSL(X1(:,6), ones(length(t_ads),1).*parameters.y1_in, X1(:,4), parameters.qsb_1, parameters.qsd_1, parameters.qsb_2, parameters.qsd_2, parameters.bo_1, parameters.do_1, parameters.bo_2, parameters.do_2, parameters.delUb_1, parameters.delUd_1, parameters.delUb_2, parameters.delUd_2);
    else
        [q1_starvals, q2_starvals] = SSLSTA(P_cycle, X_cycle(:,1), X_cycle(:,4), parameters);
        [q1_starvalsAds1, q2_starvalsAds] = SSLSTA(X1(:,6), X1(:,1), X1(:,4), parameters);
        [q1_starvalsAds, q2_starvalsAds1] = SSLSTA(X1(:,6), ones(length(t_ads),1).*parameters.y1_in, X1(:,4), parameters);
    end
    q1_starvals(1:length(t1)) = q1_starvalsAds;
    q2_starvals(1:length(t_ads)) = q2_starvalsAds;
    % plot(t_cycle./t0, q1_starvals,'k--','LineWidth', 3,'DisplayName','q1*');
    if parameters.plot0D;
        plot(t_cycle./t0, X_cycle(:,2)./1, 'b-','LineWidth', 1.5 , 'DisplayName','0-D','LineStyle',linestyleVal); hold on;
    end
    title('Mean Adsorbed amount of CO_{2} [mol/kg]')
    xline(t_ads_end./t0,  'k--', 'LineWidth', 0.8,'HandleVisibility','off');
    xline(t_blo_end./t0,  'k--', 'LineWidth', 0.8 ,'HandleVisibility','off');
    xline(t_evac_end./t0, 'k--', 'LineWidth', 0.8 ,'HandleVisibility','off');
    box on; grid off; set(gca,"LineWidth",2,"FontName","times new roman","FontSize",20)
    xlim([0 t_cycle(end)./t0])

    subplot(4,2,4)
    hold on; xlabel('time [s]'); hold on;
    % plot(t_cycle./t0, q2_starvals, 'k--','LineWidth', 3,'DisplayName','q2*');
    if parameters.plot0D;
        plot(t_cycle./t0, X_cycle(:,3)./1, 'r-','LineWidth', 1.5 , 'DisplayName','0-D','LineStyle',linestyleVal); hold on;
    end
    title('Mean Adsorbed amount of N_{2} [mol/kg]')
    xline(t_ads_end./t0,  'k--', 'LineWidth', 0.8, 'HandleVisibility','off');
    xline(t_blo_end./t0,  'k--', 'LineWidth', 0.8 , 'HandleVisibility','off');
    xline(t_evac_end./t0, 'k--', 'LineWidth', 0.8 , 'HandleVisibility','off');
    xlim([0 t_cycle(end)./t0])
    box on; grid off; set(gca,"LineWidth",2,"FontName","times new roman","FontSize",20)

    subplot(4,2,5)
    hold on; xlabel('time [s]'); hold on;
    if parameters.plot0D;
        plot(t_cycle./t0, X_cycle(:,4)./1, 'm-','LineWidth', 1.5, 'DisplayName','0-D','LineStyle',linestyleVal);% ylabel('T [K]'); xlabel('time [s]'); hold on;\
    end
    title('Mean Column Temperature [K]')
    xline(t_ads_end./t0,  'k--', 'LineWidth', 0.8,'HandleVisibility','off');
    xline(t_blo_end./t0,  'k--', 'LineWidth', 0.8, 'HandleVisibility','off');
    xline(t_evac_end./t0, 'k--', 'LineWidth', 0.8, 'HandleVisibility','off');
    box on; grid off; set(gca,"LineWidth",2,"FontName","times new roman","FontSize",20)
    set(gca, 'YScale', 'linear')
    xlim([0 t_cycle(end)])

    subplot(4,2,6)
    hold on
    if parameters.plot0D
        plot(t_cycle./t0, X_cycle(:,5)./1, 'g-','LineWidth', 1.5, 'DisplayName','0-D','LineStyle',linestyleVal);% ylabel('T_{w} [K]'); xlabel('time [s]'); hold on;
    end
    title('Mean Wall Temperature [K]'); xlabel('time [s]'); hold on;
    xline(t_ads_end./t0,  'k--', 'LineWidth', 0.8,'HandleVisibility','off');
    xline(t_blo_end./t0,  'k--', 'LineWidth', 0.8, 'HandleVisibility','off');
    xline(t_evac_end./t0, 'k--', 'LineWidth', 0.8, 'HandleVisibility','off');
    box on; grid off; set(gca,"LineWidth",2,"FontName","times new roman","FontSize",20)
    set(gca, 'YScale', 'linear')
    xlim([0 t_cycle(end)./t0])

    subplot(4,2,1)
    hold on
    if parameters.plot0D
        plot(t_cycle./t0, X_cycle(:,1), 'k-','LineWidth', 1.5, 'DisplayName','0-D','LineStyle',linestyleVal);% ylabel('y_{CO_{2}}')
    end
    title('Mole fraction of CO_{2} in the Outlet [-]'); xlabel('time [s]'); hold on;
    xline(t_ads_end./t0,  'k--', 'LineWidth', 0.8,'HandleVisibility','off');
    xline(t_blo_end./t0,  'k--', 'LineWidth', 0.8, 'HandleVisibility','off');
    xline(t_evac_end./t0, 'k--', 'LineWidth', 0.8, 'HandleVisibility','off');
    box on; grid off; set(gca,"LineWidth",2,"FontName","times new roman","FontSize",20)
    xlim([0 t_cycle(end)./t0])
    ylim([-0.05 1.05])

    subplot(4,2,7)
    hold on; xlabel('time [s]'); hold on;
    if parameters.plot0D
        plot(t_cycle./t0, F_cycleOut, 'r-','LineWidth', 1.5, 'DisplayName','0-D','LineStyle',linestyleVal);% ylabel('F_{total}'); xlabel('time [s]'); hold on;
        plot(t_cycle./t0, F_cycleIn, 'b-','LineWidth', 1.5, 'DisplayName','0-D','LineStyle',linestyleVal);% ylabel('F_{total}'); xlabel('time [s]'); hold on;
    end
    title('Molar Flowrate [mol/s]')
    xline(t_ads_end./t0,  'k--', 'LineWidth', 0.8,'HandleVisibility','off');
    xline(t_blo_end./t0,  'k--', 'LineWidth', 0.8, 'HandleVisibility','off');
    xline(t_evac_end./t0, 'k--', 'LineWidth', 0.8, 'HandleVisibility','off');
    box on; grid off; set(gca,"LineWidth",2,"FontName","times new roman","FontSize",20)
    xlim([0 t_cycle(end)./t0])

    subplot(4,2,8)
    hold on; xlabel('time [s]'); hold on;
    if parameters.plot0D;
        if parameters.pressType == "LPP"
            F_cycleIn_1 = [Fin_ads.*parameters.y1_in;zeros(length(t2),1);zeros(length(t3),1);Fin_pres.*parameters.y1_LPP];
        else
            F_cycleIn_1 = [Fin_ads.*parameters.y1_in;zeros(length(t2),1);zeros(length(t3),1);Fin_pres.*parameters.y1_in];
        end
        plot(t_cycle./t0, F_cycleOut.*X_cycle(:,1), 'r-','LineWidth', 1.5, 'DisplayName','0-D','LineStyle',linestyleVal);% ylabel('F_{total}'); xlabel('time [s]'); hold on;
        plot(t_cycle./t0, F_cycleIn_1, 'b-','LineWidth', 1.5, 'DisplayName','0-D','LineStyle',linestyleVal);% ylabel('F_{total}'); xlabel('time [s]'); hold on;
    end
    title('Molar Flowrate of CO_{2} [mol/s]')
    xline(t_ads_end./t0,  'k--', 'LineWidth', 0.8,'HandleVisibility','off');
    xline(t_blo_end./t0,  'k--', 'LineWidth', 0.8, 'HandleVisibility','off');
    xline(t_evac_end./t0, 'k--', 'LineWidth', 0.8, 'HandleVisibility','off');
    box on; grid off; set(gca,"LineWidth",2,"FontName","times new roman","FontSize",20)
    xlim([0 t_cycle(end)./t0])

    set(gcf,"Position",[100,100,1200 896],"Units","pixels")
    process_indicators = process_indicators';
    CCScycle = cycle;

    figure(834)
    q1vals = X_cycle(:,2);
    q2vals = X_cycle(:,3);
    dq1dt = gradient(q1vals,dt);
    dq2dt = gradient(q2vals,dt);




    Dp = parameters.Dm/parameters.tau; % Effective pore diffusivity [m2/s]

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

    if parameters.testBT
        KPIs = [purity_percentage, recovery_percentage, SEC, productivity,simTime,cycle ,X1(1,6)];
    else
        KPIs = process_indicators;
    end
elseif parameters.OptType ~= "sampling"
    % Append KPI results to rawData/<fileName>.txt
    if ~exist('rawData','dir')
        mkdir rawData
    end
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

warning('on','all')

end