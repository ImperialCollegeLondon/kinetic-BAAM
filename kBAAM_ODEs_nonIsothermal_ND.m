%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Imperial College London, United Kingdom
% Multiphase Systems Laboratory
% Year:     2025
% MATLAB:   R2024a
% Authors:  Hassan Azzan (HA), Ayca Yilmaz (AY)
%
% Purpose:
% Function that takes t (time), X (column vector), parameters (structure
% of the parameter for the process), and stepName (the indicated step for
% the cycle) as inputs and returns vector of ODEs dXdt as an output
%
% Last modified:
% - 2025-09-21, HA: Initial creation
%
% Input arguments:
%   - t: time
%   - X: column vector of state variables
%   - parameters: contains adsorbent properties and process parameters
%   - stepName: step being simulated
%
% Output arguments:
%   - dXdt: vector of time derivatives (ODEs) of y1, q1, q2, T
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function dXdt = kBAAM_ODEs_nonIsothermal_ND(t,X,parameters,stepName)
y1 = X(1) ; % Molar fraction of component 1 (Co2) from the column vector X
q1 = X(2) ; % Adsorbed phase concentration of component 1 (CO2) from the column vector X
q2 = X(3) ; % Adsorbed phase concentration of component 2 (N2) from the column vector X
T  = X(4) ; % Gas/Solid temperature from the column vector X

R = 8.314; % Universal gas constant [J/molK]

dXdt = zeros(size(X)); % initialize vector of ODEs

% Reference values for non-dimensionalization
volFluxRef = parameters.F_in./parameters.V_column;
timeRef = parameters.p_H./(R.*parameters.T_feed.*volFluxRef); 
qRef = parameters.qsb_1+parameters.qsd_1;
TRef = parameters.T_feed;
PRef = parameters.p_H;


switch stepName
    case 'ads' %for the case when the stepName is indicated as the adsorption step ('ads')
        P = parameters.P_ads(t.*timeRef)./PRef; %defines the pressure inside the column as the adsorption pressure as a function of time

        % Competitive equilibrium adsorbed amount and LDF coefficient for
        % both species at P, T, y1
        [q1_star, q2_star] = DSL(P.*PRef, parameters.y1_in, TRef, parameters.qsb_1, parameters.qsd_1, parameters.qsb_2, parameters.qsd_2, parameters.bo_1, parameters.do_1, parameters.bo_2, parameters.do_2, parameters.delUb_1, parameters.delUd_1, parameters.delUb_2, parameters.delUd_2);
        [k1, k2] = LDFCoefficient(P.*PRef,y1,T.*TRef,q1_star,q2_star,parameters);

        dq1dt = timeRef./qRef.* k1 * (q1_star - q1.*qRef); % calculates the time derivative of q1
        dq2dt = timeRef./qRef.* k2 * (q2_star - q2.*qRef); % calculates the time derivative of q2
        dPdt = 0 ; % pressure derivative with respect to time

        Fout = parameters.F_in - (1 - parameters.e_bed) * parameters.V_column * parameters.rho_s * (dq1dt + dq2dt)./(timeRef./qRef); %calculates the outlet flow rate with respect to the time gradients of q1 and q2 (mol/s)
    case 'blo' % for the case when the stepName is indicated as the blowdown step ('blo')
        P = parameters.P_blo(t.*timeRef)./PRef; % defines the pressure inside the column as the blowdown pressure as a function of time

        % Competitive equilibrium adsorbed amount and LDF coefficient for
        % both species at P, T, y1
        [q1_star, q2_star] = DSL(P.*PRef, y1, T.*TRef, parameters.qsb_1, parameters.qsd_1, parameters.qsb_2, parameters.qsd_2, parameters.bo_1, parameters.do_1, parameters.bo_2, parameters.do_2, parameters.delUb_1, parameters.delUd_1, parameters.delUb_2, parameters.delUd_2);
        [k1, k2] = LDFCoefficient(P.*PRef,y1,T.*TRef,q1_star,q2_star,parameters);

        dq1dt = timeRef./qRef.* k1 * (q1_star - q1.*qRef); % calculates the time derivative of q1
        dq2dt = timeRef./qRef.* k2 * (q2_star - q2.*qRef); % calculates the time derivative of q2
        dPdt = parameters.dPdt_blo(t.*timeRef)./(PRef./timeRef); % calculates pressure derivative in blowdown tep with respect to time

        parameters.F_in = 0;
        Fout = -(1 - parameters.e_bed) * parameters.V_column * parameters.rho_s * (dq1dt + dq2dt)./(timeRef./qRef) - (parameters.e_bed / (R*T.*TRef)) * dPdt.*(parameters.p_H./timeRef) * parameters.V_column ; %calculates the outlet flow rate of the blowdown step
    case 'evac'
        P = parameters.P_evac(t.*timeRef)./PRef; % defines the pressure inside the column as the blowdown pressure as a function of time

        % Competitive equilibrium adsorbed amount and LDF coefficient for
        % both species at P, T, y1
        [q1_star, q2_star] = DSL(P.*PRef, y1, T.*TRef, parameters.qsb_1, parameters.qsd_1, parameters.qsb_2, parameters.qsd_2, parameters.bo_1, parameters.do_1, parameters.bo_2, parameters.do_2, parameters.delUb_1, parameters.delUd_1, parameters.delUb_2, parameters.delUd_2);
        [k1, k2] = LDFCoefficient(P.*PRef,y1,T.*TRef,q1_star,q2_star,parameters);

        dq1dt = timeRef./qRef.* k1 * (q1_star - q1.*qRef); % calculates the time derivative of q1
        dq2dt = timeRef./qRef.* k2 * (q2_star - q2.*qRef); % calculates the time derivative of q2
        dPdt = parameters.dPdt_evac(t.*timeRef)./(parameters.p_H./timeRef); % calculates pressure derivative in evacuation tep with respect to time

        if parameters.heating 
            Qexternal = 60.*r_in.*pi.*2.*L.*(T.*TRef-parameters.Theat); % external heat flux if heating is used 
        end

        parameters.F_in = 0;
        Fout = -(1 - parameters.e_bed) * parameters.V_column * parameters.rho_s * (dq1dt + dq2dt)./(timeRef./qRef) - (parameters.e_bed / (R*T.*TRef)) * dPdt.*(parameters.p_H./timeRef) * parameters.V_column ;
    case 'pres'
        P = parameters.P_press(t.*timeRef)./PRef;

        % Competitive equilibrium adsorbed amount and LDF coefficient for
        % both species at P, T, y1
        [q1_star, q2_star] = DSL(P.*PRef, y1, T.*TRef, parameters.qsb_1, parameters.qsd_1, parameters.qsb_2, parameters.qsd_2, parameters.bo_1, parameters.do_1, parameters.bo_2, parameters.do_2, parameters.delUb_1, parameters.delUd_1, parameters.delUb_2, parameters.delUd_2);
        [k1, k2] = LDFCoefficient(P.*PRef,y1,T.*TRef,q1_star,q2_star,parameters);

        dq1dt = timeRef./qRef.* k1 * (q1_star - q1.*qRef); % calculates the time derivative of q1
        dq2dt = timeRef./qRef.* k2 * (q2_star - q2.*qRef); % calculates the time derivative of q2
        dPdt = parameters.dPdt_press(t.*timeRef)./(PRef./timeRef);

        parameters.F_in = (1 - parameters.e_bed) * parameters.V_column * parameters.rho_s * (dq1dt + dq2dt)./(timeRef./qRef)  + (parameters.e_bed / (R*T.*TRef)) * dPdt.*(parameters.p_H./timeRef) * parameters.V_column;  %calculates the inlet flow rate of the pressurization step
        Fout = 0;

        if parameters.pressType == "LPP" % If the pressurization is LPP set the composition of the inlet to be that of the light product from Adsorption (an input)
            parameters.y1_in = parameters.y1_LPP;
        end
end

% Analytical computation of heat of adsorption as a function of P and T
b1 = parameters.bo_1 .* exp(-parameters.delUb_1 ./ (R .* T.*TRef));
d1 = parameters.do_1 .* exp(-parameters.delUd_1 ./ (R .* T.*TRef));

b2 = parameters.bo_2 .* exp(-parameters.delUb_2 ./ (R .* T.*TRef));
d2 = parameters.do_2 .* exp(-parameters.delUd_2 ./ (R .* T.*TRef));

c1 = y1.*P.*PRef./(R.*T.*TRef);
c2 = (1-y1).*PRef./(R.*T.*TRef);

delH1 = -(parameters.qsb_1.*b1.*parameters.delUb_1./(1+b1.*c1).^2 + parameters.qsd_1.*d1.*parameters.delUd_1./(1+d1.*c1).^2)./...
    (parameters.qsb_1.*b1./(1+b1.*c1).^2+parameters.qsd_1.*d1./(1+d1.*c1).^2);

delH2 = -(parameters.qsb_2.*b2.*parameters.delUb_2./((1+b2.*c2).^2) + parameters.qsd_2.*d2.*parameters.delUd_2./((1+d2.*c2).^2))./...
    (parameters.qsb_2.*b2./((1+b2.*c2).^2)+parameters.qsd_2.*d2./((1+d2.*c2).^2)); 

% delH1 = 31.19e3;
% delH2 = 16.38e3;


% Reciprocal of the coefficient of the temperature derivative wrt to time
coefft1 = timeRef./TRef./(((1 - parameters.e_bed)./ parameters.e_bed).*(parameters.rho_s.*parameters.cp_s + parameters.cp_a.*parameters.rho_s.*qRef.*(q1+q2)));

% Calculates the time derivative of temperature with respect to time (K/s)
% (dimensionless)

h_in = 8.6; % inside heat transfer coefficient for the bed (Haghpanah 2014)
Qexternal = 2.*h_in./parameters.r_in./parameters.e_bed.*(T.*TRef-TRef); % heat transfer with the ambient

dTdt = coefft1.*(...
    + 1*((parameters.F_in)./parameters.V_column.*parameters.cp_g.*(TRef-TRef.*T)) ...
    - 1*parameters.cp_g./R.*dPdt.*PRef./timeRef ...
    - 1*(((1 - parameters.e_bed)./ parameters.e_bed).* parameters.cp_a.*parameters.rho_s.*T.*TRef.*qRef./timeRef.*(dq1dt + dq2dt)) + ...
    + 1*(((1 - parameters.e_bed)./ parameters.e_bed).* parameters.rho_s.*qRef./timeRef.*(delH1.*dq1dt + delH2.*dq2dt))...
    - 1*Qexternal);
if parameters.modelType == "isothermal"
    dTdt = 0;
end

% Calculates the time derivative of mole fraction of CO2 with respect to time (1/s)
% (dimensionless)
dy1dt = T./(parameters.e_bed.*P).*(-parameters.e_bed.*y1./T.*dPdt + P.*y1./(T.^2).*dTdt - ((1 - parameters.e_bed) * parameters.rho_s * qRef.*R.*TRef./PRef.*dq1dt) + ...
    ((parameters.y1_in.*parameters.F_in./volFluxRef-y1.*Fout./volFluxRef)./parameters.V_column));

% Pack ODEs to output vector
dXdt(1) = dy1dt;
dXdt(2) = dq1dt;
dXdt(3) = dq2dt;
dXdt(4) = dTdt;
end



