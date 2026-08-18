%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Imperial College London, United Kingdom
% Multiphase Systems Laboratory
% Year:     2025
% MATLAB:   R2024a
% Authors:  Hassan Azzan (HA), Ayca Yilmaz (AY)
%
% Purpose:
% Function that takes dimensionless t (time), X (column vector of states), parameters (structure
% of the parameter for the process), and stepName (the indicated step for
% the cycle) as inputs and returns vector of ODEs dXdt as an output. Same
% script is used to simulate both isothermal and non-isothermal models.
%
% Last modified:
% - 2025-10-09, HA: Add wall energy balance
% - 2025-09-21, HA: Initial creation
%
% Input arguments:
%   - t: dimensionless time
%   - X: column vector of dimensionless state variables
%   - parameters: contains adsorbent properties and process parameters
%   - stepName: step being simulated
%
% Output arguments:
%   - dXdt: vector of time derivatives (ODEs) of y1, q1, q2, T, Tw
%
% Dependencies:
%   - DSL.m
%   - LDFCoefficient.m
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function dXdt = kBAAM_ODEs_nonIsothermal_ND_nodP(t,X,parameters,stepName)
% y1 = max(X(1),eps) ; % Molar fraction of component 1 from the column vector X
y1 = X(1)  ; % Molar fraction of component 1 from the column vector X
q1 = X(2) ; % Adsorbed phase concentration of component 1 from the column vector X
q2 = X(3) ; % Adsorbed phase concentration of component 2 from the column vector X
T  = X(4) ; % Gas/Solid temperature from the column vector X
Tw = X(5) ; % Column wall temperature from the column vector X
P = X(6) ; % Column pressure from the column vector X
R = 8.3145; % Universal gas constant [J/molK]

dXdt = zeros(size(X)); % initialize vector of ODEs

% Reference values for non-dimensionalization
% volFluxRef = parameters.volFluxRef; % reference molar flowrate per unit volume [mol/m3s]
timeRef = parameters.timeRef;  % reference time [s]
qRef = parameters.qRef; % reference adsorbed amount [mol/kg]
TRef = parameters.TRef; % reference temperature [K]
TwRef = parameters.TwRef; % reference temperature [K]
PRef = parameters.PRef; % reference pressure [Pa]

Qheat = 0; % heat input due to external heating [W/m3]

dy1dt = 0;
dq1dt = 0;
dq2dt = 0;
dTdt = 0;
dTwdt = 0;
dPdt = 0;

% parameters.rp = 1e-6;


switch stepName
    case 'ads' %for the case when the stepName is indicated as the adsorption step ('ads')
        P_out = parameters.P_ads(t.*timeRef)./PRef; % dimensionless pressure inside the column as the adsorption pressure as a function of time

        % P = max(P,P_out);
        if ~parameters.pressureDrop
            P = P_out;
        end
        % parameters.h_out = 5;
        % parameters.h_in = 5;
        % Competitive equilibrium adsorbed amount and LDF coefficient for
        % both species at P, T, y1.
        % The driving force for adsorption step is given by the difference
        % between instantaneous adsorbed amount and the equilibrium amount
        % at FEED COMPOSITION (y1_in) at P,T.

        if ~parameters.SSLSTA
            [q1_starIn, q2_starIn] = DSL(P.*PRef, parameters.y1_in, T.*TRef, parameters.qsb_1, parameters.qsd_1, parameters.qsb_2, parameters.qsd_2, parameters.bo_1, parameters.do_1, parameters.bo_2, parameters.do_2, parameters.delUb_1, parameters.delUd_1, parameters.delUb_2, parameters.delUd_2);
            [q1_start, q2_start] = DSL(P.*PRef, y1, T.*TRef, parameters.qsb_1, parameters.qsd_1, parameters.qsb_2, parameters.qsd_2, parameters.bo_1, parameters.do_1, parameters.bo_2, parameters.do_2, parameters.delUb_1, parameters.delUd_1, parameters.delUb_2, parameters.delUd_2);
        else
            [q1_starIn, q2_starIn] = SSLSTA(P.*PRef, parameters.y1_in, T.*TRef, parameters);
            [q1_start, q2_start] = SSLSTA(P.*PRef, y1, T.*TRef, parameters);
        end

        [k1In, ~] = LDFCoefficient(P.*PRef,y1,T.*TRef, q1_starIn   ,q2_starIn   ,parameters);
        [k1t, k2t] = LDFCoefficient(P.*PRef,y1,T.*TRef, q1_start   ,q2_start    ,parameters);

        dq1dt =  timeRef./qRef.* k1In * (q1_starIn   - q1.*qRef); % calculates the time derivative of q1
        dq2dt =  timeRef./qRef.* k2t  * (q2_start    - q2.*qRef);

        % if parameters.cCSTR || parameters.processType == "Resin" || parameters.processType == "ResinSens"
        if parameters.cCSTR  
            dq1dt =  timeRef./qRef.* (k1t * (q1_start    - q1.*qRef)); % calculates the time derivative of q1
            dq2dt =  timeRef./qRef.* (k2t * (q2_start    - q2.*qRef)); % calculates the time derivative of q1
        end


        outletFlag = 1;

        % parameters.F_in = parameters.volFlowin.*(2.*P-P_out).*PRef./(R.*1.*TRef); % inlet molar flowrate (ideal gas) [mol/s]
        parameters.F_in = parameters.volFlowin.*(2.*P-P_out).*PRef./(R.*1.*TRef); % inlet molar flowrate (ideal gas) [mol/s]

        if parameters.pressureDrop
            voutHalf = (-1./(parameters.L)).*(4./150./1.72e-5).*(parameters.e_bed./(1 - parameters.e_bed)).^2.*parameters.rp.^2.*2.*(P_out.*PRef-P.*PRef);
            FoutHalf = ((P_out).*PRef.*parameters.A_in.*parameters.e_bed./(R.*(T).*TRef)).*voutHalf;
            dPdt = (timeRef./PRef).* ...
                + ( (P.*PRef)./(T.*TRef).*(TRef./timeRef) * dTdt ...
                -   R.*T.*TRef./1.*(1 - parameters.e_bed)./parameters.e_bed * parameters.rho_s * qRef./timeRef * (dq1dt + dq2dt) ...
                +   R.*T.*TRef./1.*((parameters.F_in - 2.*FoutHalf)./(parameters.e_bed.*parameters.V_column)));
        else
            dPdt = 0 ; % pressure derivative with respect to time [dimensionless]
        end
    case 'blo' % for the case when the stepName is indicated as the blowdown step ('blo')
        P_out = parameters.P_blo(t.*timeRef)./PRef; % dimensionless pressure inside the column as a function of time

        if ~parameters.pressureDrop
            P = P_out;
        end

        % Competitive equilibrium adsorbed amount and LDF coefficient for
        % both species at P, T, y1
        % The driving force for adsorption step is given by the difference
        % between instantaneous adsorbed amount and the equilibrium amount
        % at OUTLET COMPOSITION (y1) at P,T.
        if ~parameters.SSLSTA
            [q1_star , q2_star ] = DSL(P.*PRef, y1, T.*TRef, parameters.qsb_1, parameters.qsd_1, parameters.qsb_2, parameters.qsd_2, parameters.bo_1, parameters.do_1, parameters.bo_2, parameters.do_2, parameters.delUb_1, parameters.delUd_1, parameters.delUb_2, parameters.delUd_2);
            [q1_starMax , q2_starMax ] = DSL(P.*PRef,  parameters.y1init, T.*TRef, parameters.qsb_1, parameters.qsd_1, parameters.qsb_2, parameters.qsd_2, parameters.bo_1, parameters.do_1, parameters.bo_2, parameters.do_2, parameters.delUb_1, parameters.delUd_1, parameters.delUb_2, parameters.delUd_2);
        else
            [q1_star , q2_star ] = SSLSTA(P.*PRef, y1, T.*TRef, parameters);
        end

        [k1, k2] = LDFCoefficient(P.*PRef,y1,T.*TRef,q1_star,q2_star,parameters);
        [k1b, k2b] = LDFCoefficient(P.*PRef,parameters.y1init,T.*TRef,parameters.q1init,parameters.q2init,parameters);
        dq1dt = timeRef./qRef.* k1b * (q1_star - q1.*qRef); % calculates the time derivative of q1
        dq2dt = timeRef./qRef.* k2b * (q2_star - q2.*qRef); % calculates the time derivative of q2

        outletFlag = 1;

        parameters.F_in = 0;
        parameters.y1_in = 0;

        if parameters.pressureDrop
            voutHalf = (-1./(parameters.L) .*(4./150./1.72e-5).*(parameters.e_bed./(1 - parameters.e_bed)).^2.*parameters.rp.^2.*(P_out.*PRef-P.*PRef));
            FoutHalf = (P.*PRef.*parameters.A_in.*parameters.e_bed./(R.*(T).*TRef)).*voutHalf;
            dPdt = (timeRef./PRef).* ...
                + ( (P.*PRef)./(T.*TRef) * dTdt.*(TRef./timeRef) ...
                -   R.*T.*TRef./1.*(1 - parameters.e_bed)./parameters.e_bed * parameters.rho_s * qRef./timeRef * (dq1dt + dq2dt) ...
                +   R.*T.*TRef./1.*((parameters.F_in - 1.*FoutHalf)./(parameters.e_bed.*parameters.V_column)));
        else
            dPdt = parameters.dPdt_blo(t.*timeRef)./(PRef./timeRef); % calculates pressure derivative in blowdown tep with respect to time
        end
    case 'evac'
        P_out = parameters.P_evac(t.*timeRef)./PRef; % dimensionless pressure inside the column as a function of time

        if ~parameters.pressureDrop
            P = P_out;
        end
        % Competitive equilibrium adsorbed amount and LDF coefficient for
        % both species at P, T, y1
        % at OUTLET COMPOSITION (y1) at P,T.
        if ~parameters.SSLSTA
            [q1_star , q2_star ] = DSL(P.*PRef, y1, T.*TRef, parameters.qsb_1, parameters.qsd_1, parameters.qsb_2, parameters.qsd_2, parameters.bo_1, parameters.do_1, parameters.bo_2, parameters.do_2, parameters.delUb_1, parameters.delUd_1, parameters.delUb_2, parameters.delUd_2);
        else
            [q1_star , q2_star ] = SSLSTA(P.*PRef, y1, T.*TRef, parameters);
        end

        [k1, k2] = LDFCoefficient(P.*PRef,y1,T.*TRef,q1_star,q2_star,parameters);

        dq1dt = timeRef./qRef.* k1 * (q1_star - q1.*qRef); % calculates the time derivative of q1
        dq2dt = timeRef./qRef.* k2 * (q2_star - q2.*qRef); % calculates the time derivative of q2

        outletFlag = 1;

        if parameters.heating
            if T.*TRef < parameters.Theat
                Qheat = parameters.heatPowerDensity.*(parameters.Theat-T.*TRef)./(parameters.Theat-TRef)./(parameters.r_out-parameters.r_in); % external heat flux if heating is used
            else
                Qheat = 0;
            end
        end

        parameters.F_in = 0;
        parameters.y1_in = 0;
        if parameters.pressureDrop
            voutHalf = (-2./(parameters.L) .*(4./150./1.72e-5).*(parameters.e_bed./(1 - parameters.e_bed)).^2.*parameters.rp.^2.*(P_out.*PRef-P.*PRef));
            FoutHalf = (P.*PRef.*parameters.A_in.*parameters.e_bed./(R.*(T).*TRef)).*voutHalf;
            dPdt = (timeRef./PRef).* ...
                + ( P.*PRef./(T.*TRef) * dTdt.*(TRef./timeRef) ...
                -   R.*T.*TRef./1.*(1 - parameters.e_bed)./parameters.e_bed * parameters.rho_s * qRef./timeRef * (dq1dt + dq2dt) ...
                +   R.*T.*TRef./1.*((parameters.F_in - 1.*FoutHalf)./(parameters.e_bed.*parameters.V_column)));
        else
            dPdt = parameters.dPdt_evac(t.*timeRef)./(PRef./timeRef); % calculates pressure derivative in evacuation tep with respect to time
        end
    case 'pres'
        P_out = parameters.P_press(t.*timeRef)./PRef; % dimensionless pressure inside the column as a function of time

        if ~parameters.pressureDrop
            P = P_out;
        end

        % Competitive equilibrium adsorbed amount and LDF coefficient for
        % both species at P, T, y1
        % at INSTANTANEOUS COMPOSITION (y1) at P,T.
        if ~parameters.SSLSTA
            [q1_star , q2_star ] = DSL(P.*PRef, y1, T.*TRef, parameters.qsb_1, parameters.qsd_1, parameters.qsb_2, parameters.qsd_2, parameters.bo_1, parameters.do_1, parameters.bo_2, parameters.do_2, parameters.delUb_1, parameters.delUd_1, parameters.delUb_2, parameters.delUd_2);
        else
            [q1_star , q2_star ] = SSLSTA(P.*PRef, y1, T.*TRef, parameters);
        end

        [k1, k2] = LDFCoefficient(P.*PRef,y1, T.*TRef,q1_star,q2_star,parameters);

        dq1dt = timeRef./qRef.* k1 * (q1_star - q1.*qRef); % calculates the time derivative of q1
        dq2dt = timeRef./qRef.* k2 * (q2_star - q2.*qRef); % calculates the time derivative of q2

        outletFlag = 0;
        FoutHalf = 0;
        if parameters.pressureDrop

            vinHalf = -(-2./parameters.L .*(4./150./1.72e-5).*(parameters.e_bed./(1 - parameters.e_bed)).^2.*parameters.rp.^2.*(P_out.*PRef-P.*PRef));
            parameters.F_in = (P_out.*PRef.*parameters.A_in.*parameters.e_bed./(R.*TRef)).*vinHalf;
            dPdt = (timeRef./PRef).* ...
                + ( (P.*PRef)./(T.*TRef).*(TRef./timeRef) * dTdt ...
                -   R.*T.*TRef./1.*(1 - parameters.e_bed)./parameters.e_bed * parameters.rho_s * qRef./timeRef * (dq1dt + dq2dt) ...
                +   R.*T.*TRef./1.*((parameters.F_in - FoutHalf)./(parameters.e_bed.*parameters.V_column)));
        else
            dPdt = parameters.dPdt_press(t.*timeRef)./(PRef./timeRef);
        end
        if parameters.pressType == "LPP" % If the pressurization is LPP set the composition of the inlet to be that of the light product from Adsorption (an input)
            parameters.y1_in = parameters.y1_LPP;
        end
                Fout = 0;
end

    
if string(stepName) == "ads" && ~parameters.cCSTR
    y1val = parameters.y1_in;
else
    y1val = y1;
end

if parameters.SSLSTA
    [delH1,delH2] = computeSSLSTAHeatBinary(P.*PRef./1e5, y1, T.*TRef, parameters);
else
    [delH1,delH2] = computeDSLHeatUnary(P, y1val, T, PRef, TRef, parameters);
end

if parameters.processType == "Resin" || parameters.processType == "ResinSens"
    delH1 = -parameters.delUb_1;
    delH2 = -parameters.delUb_2;
end

switch stepName
    case 'ads'
        Fout = parameters.F_in - (1 - parameters.e_bed) * parameters.V_column * parameters.rho_s * (dq1dt + dq2dt)./(timeRef./qRef) ...
            - (parameters.e_bed / (R*T.*TRef)) * dPdt.*(PRef./timeRef) * parameters.V_column ...
            + (parameters.e_bed .* P / (R*(T.*TRef).^2)) * dTdt.*(TRef./timeRef) * parameters.V_column ; %calculates the outlet flow rate with respect to the time gradients of q1 and q2 (mol/s)
    case 'blo'
        Fout = parameters.F_in - (1 - parameters.e_bed) * parameters.V_column * parameters.rho_s * (dq1dt + dq2dt)./(timeRef./qRef) ...
            - (parameters.e_bed / (R*T.*TRef)) * dPdt.*(PRef./timeRef) * parameters.V_column ...
            + (parameters.e_bed .* P / (R*(T.*TRef).^2)) * dTdt.*(TRef./timeRef) * parameters.V_column ; %calculates the outlet flow rate with respect to the time gradients of q1 and q2 (mol/s)
    case 'evac'
        Fout = parameters.F_in - (1 - parameters.e_bed) * parameters.V_column * parameters.rho_s * (dq1dt + dq2dt)./(timeRef./qRef) ...
            - (parameters.e_bed / (R*T.*TRef)) * dPdt.*(PRef./timeRef) * parameters.V_column ...
            + (parameters.e_bed .* P / (R*(T.*TRef).^2)) * dTdt.*(TRef./timeRef) * parameters.V_column ; %calculates the outlet flow rate with respect to the time gradients of q1 and q2 (mol/s)   
    case 'pres'
        parameters.F_in = (1 - parameters.e_bed) * parameters.V_column * parameters.rho_s * (dq1dt + dq2dt)./(timeRef./qRef)  + (parameters.e_bed / (R*T.*TRef)) * dPdt.*(PRef./timeRef) * parameters.V_column - (parameters.e_bed .* P / (R*(T.*TRef).^2)) * dTdt.*(TRef./timeRef) * parameters.V_column ;  %calculates the inlet flow rate of the pressurization step
end

% Calculates the derivative of temperature with respect to time (K/s) (dimensionless)
if parameters.modelType == "isothermal"
    dTdt = 0; % set temperature derivative to 0 if isothermal
else
    coefft1 = timeRef./TRef./(((1 - parameters.e_bed)./ parameters.e_bed).*(parameters.rho_s.*parameters.cp_s + parameters.cp_a.* parameters.rho_s.*qRef.*(q1+q2))); % Reciprocal of the coefficient of the temperature derivative wrt to time
    dTdt = coefft1.*(...
        +  (parameters.F_in.*(TRef) - Fout.*T*TRef) ./(parameters.V_column.*parameters.e_bed).*parameters.cp_g ...
        -  parameters.cp_g./R.*PRef./timeRef.*dPdt ...
        -  (((1 - parameters.e_bed)./ parameters.e_bed).* parameters.cp_a.*parameters.rho_s.*T.*TRef.*qRef./timeRef.*(dq1dt + dq2dt)) ...
        +  (((1 - parameters.e_bed)./ parameters.e_bed).* parameters.rho_s.*qRef./timeRef.*(delH1.*dq1dt + delH2.*dq2dt)) ...
        -  (2.*parameters.h_in./parameters.r_in./parameters.e_bed.*(T.*TRef-Tw.*TRef)));
end

% Calculates the derivative of wall temperature with respect to time (K/s) (dimensionless)
if parameters.modelType == "isothermal"
    dTwdt = 0; % set wall temperature derivative to 0 if isothermal
else
    dTwdt = (timeRef./TwRef)./(parameters.rho_w.*parameters.cp_w).*(+2.*parameters.h_in.*parameters.r_in./(parameters.r_out.^2-parameters.r_in.^2).*(T.*TRef-Tw.*TRef) - 2.*parameters.h_out.*parameters.r_out./(parameters.r_out.^2-parameters.r_in.^2).*(Tw.*TRef-TRef)+Qheat);
end

dy1dt = (timeRef./1).* ...
    + ( y1./(T.*TRef).*(TRef./timeRef) * dTdt ...
    -   y1./(P.*PRef).*(PRef./timeRef) * dPdt ...
    -   R.*T.*TRef./(P.*PRef).*(1 - parameters.e_bed)./parameters.e_bed * parameters.rho_s * qRef./timeRef * (dq1dt) ...
    +   R.*T.*TRef./(P.*PRef).*((parameters.y1_in.*parameters.F_in-y1.*Fout  )./(parameters.e_bed.*parameters.V_column)));

% Pack ODEs to output vector
dXdt(1) = dy1dt;
dXdt(2) = dq1dt;
dXdt(3) = dq2dt;
dXdt(4) = dTdt;
dXdt(5) = dTwdt;
dXdt(6) = dPdt;
end