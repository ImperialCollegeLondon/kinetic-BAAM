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
%   - kBAAM_ODEs_nonIsothermal_ND.m
%   - DSL.m
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [KPIs] = kBAAM_Outputs_nonIsothermal(parameters,varargin)

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
        if parameters.outputType == "opt"
            parameters.p_I =  10.^parameters.p_I;
        end
    elseif parameters.processType == "PVSA"
        parameters.v_in = theta(1);
        parameters.p_I = theta(2);
        parameters.t_ads = theta(3);
        parameters.t_blo = theta(4);
        parameters.t_evac = theta(5);
        % parameters.V_column = theta(6);
        parameters.p_H = theta(6);
        if parameters.outputType == "opt"
            parameters.p_I =  10.^parameters.p_I;
            parameters.p_H =  10.^parameters.p_H;
            % parameters.v_in = parameters.v_in.*parameters.p_H./1e5;
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
        % parameters.bo_1 = 1.21e-16.*parameters.qsb_1;
        parameters.bo_1 = 1.06e-16.*parameters.qsb_1;
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
        % parameters.LDF = 0.0717.*exp(-1.316.*parameters.qsb_1);
        parameters.LDF = 0.132.*exp(-2.076.*parameters.qsb_1);
    end
end

Rg = 8.3145;

if parameters.processType == "Resin"
    dt = 0.5;
else
    dt = 0.1; % [s]
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

parameters.refVals = [1,(parameters.qsb_1+parameters.qsd_1), (parameters.qsb_1+parameters.qsd_1), parameters.T_feed, parameters.T_feed]; % vector of reference values of state variables for non-dimensionalization

parameters.Lbyr = 6; % aspect ratio length by inner radius
parameters.r_in = (parameters.V_column./(parameters.Lbyr.*pi)).^(1./3); % inner radius of column [m]
try
    parameters.r_out = parameters.r_in + parameters.t_wall; % outer radius of column [m]
catch
    parameters.r_out = parameters.r_in + 0.003; % outer radius of column [m]
end
parameters.L = parameters.Lbyr.*parameters.r_in; % length of column [m]
A_in = parameters.r_in.^2.*pi; % inner cross sectional area of column [m2]
volFlowin = parameters.v_in.*A_in; % inlet volumetric flowrate [m3/s]
parameters.F_in = volFlowin.*parameters.p_H./(Rg.*parameters.T_feed); % inlet molar flowrate (ideal gas) [mol/s]
volFluxRef = parameters.F_in./parameters.V_column; % reference molar flowrate per unit volume [mol/m3s]
timeRef = parameters.p_H./(Rg.*parameters.T_feed.*volFluxRef);  % reference time [s]


deltaP = 150./4.*1./parameters.rp.^2.*((1-parameters.e_bed)./parameters.e_bed).^2.*1.72e-5.*parameters.v_in.*parameters.L;

% pressure profiles and derivatives for overall material balance
parameters.lambda = 0.5;  % rate constant for vacuum pump, lambda [1/s]
parameters.P_ads = @(t)parameters.p_H; % pressure profile for adsorption (constant)
parameters.P_blo = @(t)parameters.p_I+(parameters.p_H-parameters.p_I)*exp(-parameters.lambda*t); % pressure profile for blowdown
parameters.P_initL = parameters.P_blo(parameters.t_blo); % pressure at the end of blowdown
parameters.P_evac = @(t)parameters.p_L+(parameters.P_initL-parameters.p_L)*exp(-parameters.lambda*t); % pressure profile for evacuation
parameters.P_initR = parameters.P_evac(parameters.t_evac); % pressure at the end of evacuation
parameters.P_press = @(t)parameters.p_H+(parameters.P_initR-parameters.p_H)*exp(-3*t); % pressure profile for pressurization
parameters.dPdt_blo = @(t)-parameters.lambda*(parameters.p_H-parameters.p_I)*exp(-parameters.lambda*t); % time derivative of pressure profile for blowdown
parameters.dPdt_evac =  @(t)-parameters.lambda*(parameters.p_I-parameters.p_L)*exp(-parameters.lambda*t); % time derivative of pressure profile for evacuation
parameters.dPdt_press = @(t)-parameters.lambda*(parameters.p_L-parameters.p_H)*exp(-parameters.lambda*t); % time derivative of pressure profile for pressurization

%% Initial condition for matrix of solution states
y1Init = [parameters.y1_in]; % initial mole fraction of component 1 in bed [-]
[q1Init, q2Init] = DSL(parameters.p_H, y1Init, parameters.T_feed, parameters.qsb_1, parameters.qsd_1, parameters.qsb_2, parameters.qsd_2, parameters.bo_1, parameters.do_1, parameters.bo_2, parameters.do_2, parameters.delUb_1, parameters.delUd_1, parameters.delUb_2, parameters.delUd_2); % initial adsorbed amounts in bed [mol/kg]
X0 = [y1Init; q1Init; q2Init; parameters.T_feed; parameters.T_feed]./parameters.refVals'; % dimensionless vector of initial states

%% Cyclic steady state simulation
max_no_Cycles = 100; % Maximum number of cycles to run to test CSS

temp_check = zeros(5, 1);

cycle = 1; % initialize cycle number
warning('off','all')
tic
% try
    while  cycle < max_no_Cycles && mean(temp_check) < 1
        options = odeset('RelTol', 1e-4, 'AbsTol', 1e-4, 'MaxOrder', 2); % sets the levels of relative tolerance, absolute tolerance and maximum order of the ode15 for our system


        [t1, X1] = ode15s(@(t,X) kBAAM_ODEs_nonIsothermal_ND(t,X,parameters,'ads'), t_ads./(timeRef), X0, options); %t1 is the time point at which the solution is evaluated, X1 is the solution states for adsorption step
        t1 = t1.*timeRef;
        X1(X1<0) = 0;
        X1(X1(:,1)>1,1) = 1;
        X0 = X1(end,:)';
        X1 = X1.*parameters.refVals;

        Fout_ads = parameters.F_in - (1 - parameters.e_bed) * parameters.V_column * parameters.rho_s * (gradient(X1(:,2),dt.*timeRef) + gradient(X1(:,3),dt.*timeRef));
        Fout_ads(Fout_ads<0) = 0;
        F_1_out_ads = Fout_ads.*X1(:,1);

        mol_1_out_ads = trapz(t1.*timeRef,F_1_out_ads); moltot_out_ads = trapz(t1.*timeRef,Fout_ads);
        parameters.y1_LPP = mol_1_out_ads./moltot_out_ads;

        [t2, X2] = ode15s(@(t,X) kBAAM_ODEs_nonIsothermal_ND(t,X,parameters,'blo'), t_blo./(timeRef), X0, options);%t1 is the time point at which the solution is evaluated, X1 is the solution states for adsorption step
        t2 = t2.*timeRef;
        X2(X2<0) = 0;
        X2(X2(:,1)>1,1) = 1;
        X0 = X2(end,:)';
        X2 = X2.*parameters.refVals;

        [t3, X3] = ode15s(@(t,X) kBAAM_ODEs_nonIsothermal_ND(t,X,parameters,'evac'), t_evac./(timeRef), X0, options) ;%t1 is the time point at which the solution is evaluated, X1 is the solution states for adsorption step
        t3 = t3.*timeRef;
        X3(X3<0) = 0;
        X3(X3(:,1)>1,1) = 1;
        X0 = X3(end,:)'; %This sets up the initial condition for the next step, by taking the final state from the previous step.%end means the last row of X1, containing all the state variables at %the t final , : takes all the columns (for all types of the state%variables
        X3 = X3.*parameters.refVals;

        if t3(end) < 0.95*parameters.t_evac
            parameters.t_evac = t3(end); % Update evacuation time to length of simulation instead of discarding. Integration seems to fail when flowrate in close to 0 and mole fraction close to 1.
            t_evac  = 0:dt:parameters.t_evac;
            parameters.P_initR = parameters.P_evac(parameters.t_evac);
            cycle = 1;
        end

        [t4, X4] = ode15s(@(t,X) kBAAM_ODEs_nonIsothermal_ND(t,X,parameters,'pres'), t_press./(timeRef), X0, options);%t1 is the time point at which the solution is evaluated, X1 is the solution states for adsorption step
        t4 = t4.*timeRef;
        X4(X4<0) = 0;
        X4(X4(:,1)>1,1) = 1;
        X0 = X4(end,:)'; %This sets up the initial condition for the next step, by taking the final state from the previous step.
        X4 = X4.*parameters.refVals;

        n_1_bd = (X2(end,1)* parameters.P_blo(t2(end)) * parameters.V_column * parameters.e_bed / (Rg * X2(end,4))) + X2(end,2)* parameters.V_column * (1-parameters.e_bed).*parameters.rho_s ;
        n_1_ads = (X1(end,1)* parameters.p_H * parameters.V_column * parameters.e_bed / (Rg * X1(end,4))) + X1(end,2)* parameters.V_column * (1-parameters.e_bed).*parameters.rho_s ;
        n_1_evac = (X3(end,1)* parameters.P_evac(t3(end)) * parameters.V_column * parameters.e_bed / (Rg * X3(end,4))) + X3(end,2)* parameters.V_column * (1-parameters.e_bed).*parameters.rho_s ;
        n_2_bd = ((1-X2(end,1))* parameters.P_blo(t2(end)) * parameters.V_column * parameters.e_bed / (Rg * X2(end,4))) + X2(end,3)* parameters.V_column * (1-parameters.e_bed).*parameters.rho_s ;
        n_2_evac = ((1-X3(end,1))* parameters.P_evac(t3(end)) * parameters.V_column * parameters.e_bed / (Rg * X3(end,4))) + X3(end,3)* parameters.V_column * (1-parameters.e_bed).*parameters.rho_s ;
        n_1_pres = (X4(end,1)* parameters.P_press(t4(end)) * parameters.V_column * parameters.e_bed / (Rg * X4(end,4))) + X4(end,2)* parameters.V_column * (1-parameters.e_bed).*parameters.rho_s ;
        % n_1_adsin = parameters.F_in * parameters.y1_in * parameters.t_ads;
        % n_1_adsin = trapz(t1,parameters.v_in.*A_in.*parameters.p_H./(8.314.*X1(:,4)) * parameters.y1_in);
        n_1_adsin=n_1_ads-n_1_pres+mol_1_out_ads;
        if parameters.pressType == "LPP"
            recovery_percentage = 100 * (n_1_bd - n_1_evac) / (n_1_adsin);
        else
            recovery_percentage = 100 * (n_1_bd - n_1_evac) / (n_1_pres - n_1_evac + n_1_adsin);
        end
        purity_percentage = 100 * (n_1_bd - n_1_evac) / (n_1_bd - n_1_evac + n_2_bd - n_2_evac);

        if recovery_percentage>100
            recovery_percentage = 0;
            purity_percentage = 0;
        end

        if purity_percentage>100
            recovery_percentage = 0;
            purity_percentage = 0;
        end

        cycle_time = (parameters.t_ads + parameters.t_blo + parameters.t_evac + parameters.t_press);

        productivity = (n_1_bd - n_1_evac) /(parameters.V_column.*cycle_time);

        %%  Energy Calculation
        Fout_bd = 0 - (1 - parameters.e_bed) .* parameters.V_column * parameters.rho_s .* (gradient(X2(:,2),dt) + gradient(X2(:,3),dt))-(parameters.e_bed ./ (Rg.*X2(:,4))) .* parameters.dPdt_blo(t2) .* parameters.V_column ;
        Fout_bd(Fout_bd<0) = 0;
        vout_bd = Fout_bd.*Rg.*X2(:,4)./(parameters.P_blo(t2))./A_in;
        Fout_evac = 0 - (1 - parameters.e_bed) .* parameters.V_column * parameters.rho_s .* (gradient(X3(:,2),dt) + gradient(X3(:,3),dt))-(parameters.e_bed ./ (Rg.*X3(:,4))) .* parameters.dPdt_evac(t3) .* parameters.V_column ;
        vout_evac = Fout_evac.*Rg.*X3(:,4)./(parameters.P_evac(t3))./A_in;
        Fin_pres = (1 - parameters.e_bed) .* parameters.V_column * parameters.rho_s .*  (gradient(X4(:,2),dt) + gradient(X4(:,3),dt))  + (parameters.e_bed ./ (Rg.*X4(:,4))) .* parameters.dPdt_press(t4) .* parameters.V_column;  %calculates the inlet flow rate of the pressurization step
        vin_pres = Fin_pres.*Rg.*X4(:,4)./(parameters.P_press(t4))./A_in;
        % Fin_ads = parameters.F_in;  %calculates the inlet flow rate of the pressurization step

        eta_bd = 0.8.*(19.55.*parameters.P_blo(t2).*1e-5./(1+19.55.*parameters.P_blo(t2).*1e-5));
        eta_evac = 0.8.*(19.55.*parameters.P_evac(t3).*1e-5./(1+19.55.*parameters.P_evac(t3).*1e-5));
        % eta_press = 0.8.*(19.55.*parameters.P_press(t4).*1e-5./(1+19.55.*parameters.P_press(t4).*1e-5));
        % eta_ads = 0.8.*(19.55.*parameters.p_H.*1e-5./(1+19.55.*parameters.p_H.*1e-5));
        eta_press = 0.8;
        eta_ads = 0.8;


        % EC_BD   = trapz(t2,1./eta_bd.*Fout_bd  .*Rg.*X2(:,4).*(1.4./0.4).*((1e5./min(1e5,parameters.P_blo(t2))).^(0.4./1.4)-1));
        % EC_EVAC = trapz(t3,1./eta_evac.*Fout_evac.*Rg.*X3(:,4).*(1.4./0.4).*((1e5./min(1e5,parameters.P_evac(t3))).^(0.4./1.4)-1));
        % EC_PRES = trapz(t4,1./eta_press.*Fin_pres .*Rg.*X4(:,4).*(1.4./0.4).*((max(1e5,parameters.P_press(t4))./1e5).^(0.4./1.4)-1));
        % EC_ADS = (t1(end).*1./eta_ads.*Fin_ads .*Rg.*parameters.T_feed.*(1.4./0.4).*((max(1e5,parameters.p_H)./1e5).^(0.4./1.4)-1));
        EC_BD   = trapz(t2,1./eta_bd.*vout_bd.*A_in.*parameters.e_bed.*parameters.P_blo(t2).*(1.4./0.4).*((1e5./min(1e5,parameters.P_blo(t2))).^(0.4./1.4)-1));
        EC_EVAC = trapz(t3,1./eta_evac.*vout_evac.*A_in.*parameters.e_bed.*parameters.P_evac(t3).*(1.4./0.4).*((1e5./min(1e5,parameters.P_evac(t3))).^(0.4./1.4)-1));
        EC_FAN = t1(end).*1./eta_press.*parameters.v_in.*A_in.*parameters.e_bed.*(parameters.p_H+deltaP).*(1.4./0.4).*((max(1e5,(parameters.p_H+deltaP))./1e5).^(0.4./1.4)-1);
        
        

        if parameters.heating
            % EC_HEAT = 120.*(parameters.r_out).*pi.*2.*parameters.L.*t3(end);
            % EC_HEAT = trapz(t3,1200.*(parameters.r_out).*pi.*2.*parameters.L.*(X3(:,4)-X3(:,5)));
            heatFlag = X3(:,4) < parameters.Theat;
            Qheat = trapz(t3,parameters.heatPowerDensity.*(parameters.Theat-X3(:,4))./(parameters.Theat-parameters.T_feed).*parameters.r_out.*2.*parameters.L.*heatFlag); % external heat flux if heating is used
            EC_HEAT = Qheat;
        else
            EC_HEAT = 0;
        end
        EC_PRES = trapz(t4,1./eta_press.*vin_pres.*A_in.*parameters.e_bed.*parameters.P_press(t4).*(1.4./0.4).*((max(1e5,parameters.P_press(t4))./1e5).^(0.4./1.4)-1));
        EC_ADS = (t1(end).*1./eta_ads.*parameters.v_in.*A_in.*parameters.e_bed.*parameters.p_H.*(1.4./0.4).*((max(1e5,parameters.p_H)./1e5).^(0.4./1.4)-1));

        SEC = (EC_PRES + EC_ADS + EC_BD + EC_EVAC + EC_HEAT + EC_FAN)./((n_1_bd - n_1_evac).*0.044);

        recovery_percentageValues(cycle) = recovery_percentage;
        purity_percentageValues(cycle) = purity_percentage;
        productivity_Values(cycle) = productivity;
        SEC_Values(cycle) = SEC;

        process_indicators = [purity_percentageValues; recovery_percentageValues; productivity_Values; SEC_Values];

        if cycle > 5
            for i = 1:4
                for k = 0:4
                    if abs(100*(process_indicators(i, cycle-k) - process_indicators(i, cycle-5))/process_indicators(i, cycle-5)) <= 0.5
                        temp_check(k+1) = 1;
                    else
                        temp_check(k+1) = 0;
                    end
                end
            end
        end

        cycle = cycle+1;

    end



 
    if purity_percentage < 0 || recovery_percentage < 0 
        SEC = 100e6;
        purity_percentage = 0;
        recovery_percentage = 0;
        productivity = 0;
    end




    if parameters.processType == "Resin"
        phi_pen = [0, 0];
        phi_pen(1) = 0.80.*((1.*max(0,(95-purity_percentage)))).^2 + (max(0,(0-recovery_percentage))).^2;
        phi_pen(2) = 0.30.*((1.*max(0,(95-purity_percentage))).^2 + (max(0,(0-recovery_percentage))).^2);
        if parameters.OptType == "Const"
            KPIs = [(-productivity+phi_pen(1)) 50*(SEC.*2.77778e-7+phi_pen(2))];
        else
            KPIs = [-recovery_percentage -purity_percentage];
        end
    else
        phi_pen = [0, 0];
        phi_pen(1) = 0.80.*((1.*max(0,(95-purity_percentage)))).^2 + (max(0,(90-recovery_percentage))).^2;
        phi_pen(2) = 0.30.*((1.*max(0,(95-purity_percentage))).^2 + (max(0,(90-recovery_percentage))).^2);
        if parameters.OptType == "Const"
            KPIs = [(-productivity+phi_pen(1)) 10*(SEC.*2.77778e-7+phi_pen(2))];
        else
            KPIs = [-recovery_percentage -purity_percentage];
        end
    end



% catch
%     SEC = 100e6;
%     purity_percentage = 0;
%     recovery_percentage = 0;
%     productivity = 0;
% 
%     if parameters.processType == "Resin"
%         if parameters.OptType == "Const"
%             KPIs = [(-productivity+phi_pen(1)) 50*(SEC.*2.77778e-7+phi_pen(2))];
%         else
%             KPIs = [-recovery_percentage -purity_percentage];
%         end
%     else
%         if parameters.OptType == "Const"
%             KPIs = [(-productivity+phi_pen(1)) 10*(SEC.*2.77778e-7+phi_pen(2))];
%         else
%             KPIs = [-recovery_percentage -purity_percentage];
%         end
%     end
% end
%%

if parameters.outputType == "plot"

    t_cycle = [t1; t2 + t1(end); t3 + t1(end) + t2(end); t4 + t1(end) + t2(end) + t3(end)]; %Shift the t2 time vector forward in time so it starts immediately after t1 ends
    X_cycle = [X1; X2; X3; X4];
    P_ads = parameters.P_ads(t1);
    P1 = repmat(P_ads, size(t1));
    P2 = parameters.P_blo(t2);
    P3 = parameters.P_evac(t3);
    P4 = parameters.P_press(t4);
    t_ads_end  = t1(end);
    t_blo_end  = t_ads_end + t2(end);
    t_evac_end = t_blo_end + t3(end);
    P_cycle = [P1; P2; P3; P4];
    figure(1);
    subplot(6,1,1)
    hold on
    plot(t_cycle, P_cycle./1e5,'-', 'Color','#0B0','LineWidth', 2); ylabel('P [bar]'); xlabel('time [s]'); hold on;%check unit
    title('Column Pressure (Cyclic Steady State)')
    xline(t_ads_end,  'k--', 'LineWidth', 1, 'HandleVisibility','off');
    xline(t_blo_end,  'k--', 'LineWidth', 1 , 'HandleVisibility','off');
    xline(t_evac_end, 'k--', 'LineWidth', 1, 'HandleVisibility','off');
    box on; grid off; set(gca,"LineWidth",2,"FontName","CMU Serif","FontSize",15)
    xlim([0 t_cycle(end)])


    subplot(6,1,2)
    hold on
    T = parameters.T_feed;
    [q1_starvals, q2_starvals] = DSL(P_cycle, X_cycle(:,1), X_cycle(:,4), parameters.qsb_1, parameters.qsd_1, parameters.qsb_2, parameters.qsd_2, parameters.bo_1, parameters.do_1, parameters.bo_2, parameters.do_2, parameters.delUb_1, parameters.delUd_1, parameters.delUb_2, parameters.delUd_2);
    [q1_starvalsAds, q2_starvalsAds] = DSL(parameters.P_ads(t_ads), parameters.y1_in, X1(:,4), parameters.qsb_1, parameters.qsd_1, parameters.qsb_2, parameters.qsd_2, parameters.bo_1, parameters.do_1, parameters.bo_2, parameters.do_2, parameters.delUb_1, parameters.delUd_1, parameters.delUb_2, parameters.delUd_2);
    [q1_starvalsAds, q2_starvalsAds] = DSL(parameters.P_ads(t_ads), parameters.y1_in, X1(:,4), parameters.qsb_1, parameters.qsd_1, parameters.qsb_2, parameters.qsd_2, parameters.bo_1, parameters.do_1, parameters.bo_2, parameters.do_2, parameters.delUb_1, parameters.delUd_1, parameters.delUb_2, parameters.delUd_2);
    % [q1_starvalsAds, q2_starvalsAds] = DSL(parameters.P_ads(t_ads), X1(:,1), T, parameters.qsb_1, parameters.qsd_1, parameters.qsb_2, parameters.qsd_2, parameters.bo_1, parameters.do_1, parameters.bo_2, parameters.do_2, parameters.delUb_1, parameters.delUd_1, parameters.delUb_2, parameters.delUd_2);
    q1_starvals(1:length(t_ads)) = q1_starvalsAds;
    q2_starvals(1:length(t_ads)) = q2_starvalsAds;
    plot(t_cycle, q1_starvals,'k--','LineWidth', 2,'DisplayName','q1*');
    plot(t_cycle, X_cycle(:,2), 'b-','LineWidth', 2 , 'DisplayName','q1'); hold on;
    ylabel('q_{CO_{2}} [mol/kg]'); xlabel('time [s]')
    title('Adsorbed amount of CO_{2} (Cyclic Steady State)')
    legend;
    xline(t_ads_end,  'k--', 'LineWidth', 1,'HandleVisibility','off');
    xline(t_blo_end,  'k--', 'LineWidth', 1 ,'HandleVisibility','off');
    xline(t_evac_end, 'k--', 'LineWidth', 1 ,'HandleVisibility','off');
    box on; grid off; set(gca,"LineWidth",2,"FontName","CMU Serif","FontSize",15)
    xlim([0 t_cycle(end)])

    subplot(6,1,3)
    hold on
    plot(t_cycle, q2_starvals, 'k--','LineWidth', 2,'DisplayName','q2*');
    plot(t_cycle, X_cycle(:,3), 'r--','LineWidth', 2 , 'DisplayName','q2'); hold on;
    ylabel('q_{N_{2}} [mol/kg]'); xlabel('time [s]')
    title('Adsorbed amount of N_{2} (Cyclic Steady State)')
    xline(t_ads_end,  'k--', 'LineWidth', 1, 'HandleVisibility','off');
    xline(t_blo_end,  'k--', 'LineWidth', 1 , 'HandleVisibility','off');
    xline(t_evac_end, 'k--', 'LineWidth', 1 , 'HandleVisibility','off');
    xlim([0 t_cycle(end)])
    box on; grid off; set(gca,"LineWidth",2,"FontName","CMU Serif","FontSize",15)
    legend;

    subplot(6,1,4)
    hold on
    plot(t_cycle, X_cycle(:,4), 'm-','LineWidth', 2); ylabel('T [K]'); xlabel('time [s]'); hold on;
    title('Column Temperature (Cyclic Steady State)')
    xline(t_ads_end,  'k--', 'LineWidth', 1,'HandleVisibility','off');
    xline(t_blo_end,  'k--', 'LineWidth', 1, 'HandleVisibility','off');
    xline(t_evac_end, 'k--', 'LineWidth', 1, 'HandleVisibility','off');
    box on; grid off; set(gca,"LineWidth",2,"FontName","CMU Serif","FontSize",15)
    set(gca, 'YScale', 'linear')
    xlim([0 t_cycle(end)])

    subplot(6,1,5)
    hold on
    plot(t_cycle, X_cycle(:,5), 'g-','LineWidth', 2); ylabel('T [K]'); xlabel('time [s]'); hold on;
    title('Wall Temperature (Cyclic Steady State)')
    xline(t_ads_end,  'k--', 'LineWidth', 1,'HandleVisibility','off');
    xline(t_blo_end,  'k--', 'LineWidth', 1, 'HandleVisibility','off');
    xline(t_evac_end, 'k--', 'LineWidth', 1, 'HandleVisibility','off');
    box on; grid off; set(gca,"LineWidth",2,"FontName","CMU Serif","FontSize",15)
    set(gca, 'YScale', 'linear')
    xlim([0 t_cycle(end)])

    subplot(6,1,6)
    hold on
    plot(t_cycle, X_cycle(:,1), 'k-','LineWidth', 2); ylabel('y_{CO_{2}}'); xlabel('time [s]'); hold on;
    title('Mole fraction of CO_{2} (Cyclic Steady State)')
    xline(t_ads_end,  'k--', 'LineWidth', 1,'HandleVisibility','off');
    xline(t_blo_end,  'k--', 'LineWidth', 1, 'HandleVisibility','off');
    xline(t_evac_end, 'k--', 'LineWidth', 1, 'HandleVisibility','off');
    box on; grid off; set(gca,"LineWidth",2,"FontName","CMU Serif","FontSize",15)
    % set(gca, 'YScale', 'l7 og')
    xlim([0 t_cycle(end)])
    ylim([-0.05 1.05])


    set(gcf,"Position",[100,100,657 896],"Units","pixels")

    CCScycle = cycle

else
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