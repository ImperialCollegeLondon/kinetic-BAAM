%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Imperial College London, United Kingdom
% Multiphase Systems Laboratory
% Year:     2025
% MATLAB:   R2024a
% Authors:  Hassan Azzan (HA), Ayca Yilmaz (AY)
%
% Purpose:
% Evaluates the RHS of the dimensionless ODEs for the non-isothermal k-BAAM with pressure
% drop. Takes dimensionless time t, state vector X, parameters struct, cycle step name,
% and optional mode flag, and returns dXdt. Handles four PSA steps (ads, blo, evac, pres)
% via a switch statement. Supports both isothermal and non-isothermal models.
%
% The DAE system is written in semi-explicit mass-matrix form M*dX/dt = f. The function
% solves this directly by computing M\f at each call, using a 6x6 LU factorisation.
% This is efficient for the 6-state system here and avoids solver overhead from the
% odeset Mass option, which offers no benefit for small, strongly state-dependent M.
%
% State vector X (all dimensionless):
%   X(1) = y1     mole fraction of component 1 (CO2)        [-]
%   X(2) = q1/qRef  adsorbed amount of component 1           [-]
%   X(3) = q2/qRef  adsorbed amount of component 2 (N2)      [-]
%   X(4) = T/TRef   gas and solid temperature                 [-]
%   X(5) = Tw/TwRef wall temperature                          [-]
%   X(6) = P/PRef   total pressure                            [-]
%
% Equilibrium model: DSL (default) or SSLSTA (parameters.SSLSTA = 1)
% Kinetics: Linear Driving Force (LDF) via LDFCoefficient.m
% Flow: Darcy's law with linear pressure profile along bed length
% Pressure drop: Ergun-based Darcy permeability precomputed in Outputs
%
% Last modified:
% - 2026-04-23, HA: Add section headers, inline comments, precomputed constants,
%                   state clamping, dimensional caching, fewer LDF calls in ads step,
%                   remove buildMassMatrix / mode-switching infrastructure
% - 2025-10-09, HA: Add wall energy balance
% - 2025-09-21, HA: Initial creation
%
% Input arguments:
%   - t:          dimensionless time [-]
%   - X:          6x1 column vector of dimensionless state variables (see above)
%   - parameters: struct of adsorbent properties and process parameters
%   - stepName:   string identifying the current cycle step ('ads','blo','evac','pres')
%
% Output arguments:
%   - dXdt: 6x1 vector of time derivatives
%
% Local functions:
%   - getEquilibriumLoadings: wraps DSL/SSLSTA isotherm call
%   - buildMassMatrix:        assembles sparse 6x6 mass matrix M (called internally)
%
% Dependencies:
%   - DSL.m
%   - SSLSTA.m
%   - LDFCoefficient.m
%   - computeDSLHeatUnary.m
%   - computeSSLSTAHeatBinaryBT.m
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function dXdt = kBAAM_ODEs_nonIsothermal_ND_dP(t,X,parameters,stepName)

% Unpack state variables
y1 = max(1e-11,min(1,X(1)));   % mole fraction of component 1 [-]
q1 = X(2);   % dimensionless adsorbed amount of component 1 [-]
q2 = X(3);   % dimensionless adsorbed amount of component 2 [-]
T  = max(X(4),1e-11);   % dimensionless temperature [-]
Tw = X(5);   % dimensionless wall temperature [-]
P  = max(X(6),1e-11);   % dimensionless pressure [-]

R = 8.3145;  % universal gas constant [J/molK]

% Use precomputed constants (set in Outputs before cycle loop)
tRef  = parameters.timeRef;
qRef  = parameters.qRef;
TRef  = parameters.TRef;
TwRef = parameters.TwRef;
PRef  = parameters.PRef;
Ab    = parameters.Ab;                    % (1-e)/e
coeff_q  = parameters.Ab_rhos_qRef_tRef;  % Ab*rho_s*qRef/tRef
cp_a  = parameters.cp_a;
cp_g  = parameters.cp_g;
e   = parameters.e_bed;
V   = parameters.V_column;
L   = parameters.L;
A   = parameters.A_in;
tRef_qRef = parameters.tRef_qRef;         % tRef/qRef
darcyK    = parameters.darcyK;            % Darcy permeability factor
cpg_eV    = parameters.cpg_eV;            % cp_g/(V*e)
two_hin   = parameters.two_hin_rin_e;     % 2*h_in/(r_in*e)

% Cache dimensional pressure and temperature (used repeatedly)
P_dim = P .* PRef;
T_dim = T .* TRef;

Qheat = 0; % external heat input [W/m3] (only nonzero in evac with heating)

% ---- Step-specific: equilibrium, kinetics, and flow ----
switch stepName
    case 'ads'
        P_out = parameters.P_ads(t.*tRef) ./ PRef;

        if parameters.cCSTR
            % cCSTR mode: both components use outlet composition equilibrium
            [q1_start, q2_start] = getEquilibriumLoadings(P, y1, T, PRef, TRef, parameters);
            [k1t, k2t] = LDFCoefficient(P_dim, y1, T_dim, q1_start, q2_start, parameters);
            dq1dt = tRef_qRef .* k1t .* (q1_start - q1.*qRef);
            dq2dt = tRef_qRef .* k2t .* (q2_start - q2.*qRef);
        else
            % Default mode: CO2 uses feed-side equilibrium, N2 uses outlet equilibrium
            [q1_starIn, q2_starIn] = getEquilibriumLoadings(P, parameters.y1_in, T, PRef, TRef, parameters);
            [q1_start, q2_start]   = getEquilibriumLoadings(P, y1, T, PRef, TRef, parameters);
            [k1In, ~]  = LDFCoefficient(P_dim, y1, T_dim, q1_starIn, q2_starIn, parameters);
            [~, k2t]   = LDFCoefficient(P_dim, y1, T_dim, q1_start,  q2_start,  parameters);

            dq1dt = tRef_qRef .* k1In .* (q1_starIn - q1.*qRef);
            dq2dt = tRef_qRef .* k2t  .* (q2_start  - q2.*qRef);
        end

        % Inlet flow from feed specification: F_in = v_in*A*e * P_in/(R*T_feed)
        % where P_in = 2*P_avg - P_out (linear pressure profile)
        F_in  = parameters.volFlowin .* (2.*P - P_out) .* PRef ./ (R .*  TRef);
        y1_in = parameters.y1_in;

        % Outlet flow from Darcy's law at product end
        v_out = (2/L) .* darcyK .* (P - P_out) .* PRef;
        Fout  = P_out *PRef .* A .* e ./ (R .* T.*TRef) .* v_out; 
        % Fout  = max(0,Fout);

    case 'blo'
        P_out = parameters.P_blo(t.*tRef) ./ PRef;

        % Equilibrium at outlet composition
        [q1_star, q2_star] = getEquilibriumLoadings(P, y1, T, PRef, TRef, parameters);

        % LDF coefficient evaluated at outlet conditions
        [k1, k2] = LDFCoefficient(P_dim, y1, T_dim, q1_star, q2_star, parameters);
        dq1dt = tRef_qRef .* k1 .* (q1_star - q1.*qRef);
        dq2dt = tRef_qRef .* k2 .* (q2_star - q2.*qRef);

        F_in  = 0;
        y1_in = 0;

        % Outlet flow from Darcy (co-current, factor of 2 for linear profile)
        v_out = (2/L) .* darcyK .* (P - P_out) .* PRef;
        Fout  = P_dim .* A .* e ./ (R .* T_dim) .* v_out;

    case 'evac'
        P_out = parameters.P_evac(t.*tRef) ./ PRef;

        % Equilibrium at outlet composition
        [q1_star, q2_star] = getEquilibriumLoadings(P, y1, T, PRef, TRef, parameters);

        [k1, k2] = LDFCoefficient(P_dim, y1, T_dim, q1_star, q2_star, parameters);
        dq1dt = tRef_qRef .* k1 .* (q1_star - q1.*qRef);
        dq2dt = tRef_qRef .* k2 .* (q2_star - q2.*qRef);

        % External heating (temperature swing)
        if parameters.heating
            if T_dim < parameters.Theat
                Qheat = parameters.heatPowerDensity .* (parameters.Theat - T_dim) ...
                    ./ (parameters.Theat - TRef) ./ (parameters.r_out - parameters.r_in);
            end
        end

        F_in  = 0;
        y1_in = 0;

        % Outlet flow from Darcy (counter-current exit at feed end)
        v_out = (2/L) .* darcyK .* (P - P_out) .* PRef;
        Fout  = P_dim .* A .* e ./ (R .* T_dim) .* v_out;

    case 'pres'
        P_out = parameters.P_press(t.*tRef) ./ PRef;

        % Equilibrium at instantaneous composition
        [q1_star, q2_star] = getEquilibriumLoadings(P, y1, T, PRef, TRef, parameters);

        [k1, k2] = LDFCoefficient(P_dim, y1, T_dim, q1_star, q2_star, parameters);
        dq1dt = tRef_qRef .* k1 .* (q1_star - q1.*qRef);
        dq2dt = tRef_qRef .* k2 .* (q2_star - q2.*qRef);

        % Inlet flow from Darcy (gas enters column)
        v_in = (2/L) .* darcyK .* (P_out - P) .* PRef;
        F_in = P_out.*PRef .* A .* e ./ (R .* TRef) .* v_in;
        Fout = 0;

        if parameters.pressType == "LPP"
            y1_in = parameters.y1_LPP;
        else
            y1_in = parameters.y1_in;
        end
end

% ---- Assemble RHS (source terms only, no cross-derivative terms) ----
f = zeros(6, 1);

% Row 1: Species balance — flow source only
f(1) = R.*T_dim ./ P_dim .* (y1_in.*F_in - y1.*Fout) ./ (e.*V);

% Rows 2,3: Adsorbed phase (LDF kinetics)
f(2) = dq1dt;
f(3) = dq2dt;

% Row 4: Energy balance — flow enthalpy + wall heat transfer
if ~parameters.isIsothermal
    f(4) = cpg_eV .* (F_in.*TRef - Fout.*T_dim) ...
         - two_hin .* (T_dim - Tw.*TwRef);
end

% Row 5: Wall energy balance
if ~parameters.isIsothermal
    f(5) = parameters.wall_prefactor .* ...
        (+parameters.wall_coeff1 .* (T_dim - Tw.*TwRef) ...
         -parameters.wall_coeff2 .* (Tw.*TwRef - TRef) ...
         + Qheat);
end

% Row 6: Overall material balance — flow source only
f(6) = R.*T_dim .* (F_in - Fout) ./ (e.*V);

M = buildMassMatrix(y1,q1,q2,T,P,parameters,R,tRef,qRef,TRef,PRef,Ab,coeff_q,cp_a,cp_g,stepName);
dXdt = M\f;
end

function [q1_star, q2_star] = getEquilibriumLoadings(P, y1, T, PRef, TRef, parameters)
if parameters.SSLSTA
    [q1_star, q2_star] = SSLSTA(P.*PRef, y1, T.*TRef, parameters);
else
    [q1_star, q2_star] = DSL(P.*PRef, y1, T.*TRef, parameters.qsb_1, parameters.qsd_1, parameters.qsb_2, parameters.qsd_2, parameters.bo_1, parameters.do_1, parameters.bo_2, parameters.do_2, parameters.delUb_1, parameters.delUd_1, parameters.delUb_2, parameters.delUd_2);
end
end

function M = buildMassMatrix(y1,q1,q2,T,P,parameters,R,tRef,qRef,TRef,PRef,Ab,coeff_q,cp_a,cp_g,stepName)
P_dim = P .* PRef;
T_dim = T .* TRef;

% Pre-allocate triplets (row, col, val) — max 14 nonzeros
ii = zeros(14,1); jj = zeros(14,1); vv = zeros(14,1);
nz = 0;

RT_PP = R * T_dim / P_dim;  % common factor

% Row 1: Species balance
nz=nz+1; ii(nz)=1; jj(nz)=1; vv(nz) = 1/tRef;
nz=nz+1; ii(nz)=1; jj(nz)=2; vv(nz) = RT_PP * coeff_q;
nz=nz+1; ii(nz)=1; jj(nz)=4; vv(nz) = -y1 / (T * tRef);
nz=nz+1; ii(nz)=1; jj(nz)=6; vv(nz) =  y1 / (P * tRef);

% Row 2: LDF component 1
nz=nz+1; ii(nz)=2; jj(nz)=2; vv(nz) = 1;

% Row 3: LDF component 2
nz=nz+1; ii(nz)=3; jj(nz)=3; vv(nz) = 1;

if ~parameters.isIsothermal
    % Effective heat capacity [J/m3K]
    Ceff = Ab * (parameters.rho_s * parameters.cp_s + cp_a * parameters.rho_s * qRef * (q1 + q2));
    
    if string(stepName) == "ads" && ~parameters.cCSTR
        y1val = parameters.y1_in;
    else
        y1val = y1;
    end
    % Heat of adsorption [J/mol]
    if parameters.SSLSTA
        [delH1, delH2] = computeSSLSTAHeatBinaryBT(P_dim./1e5, y1, T_dim, [parameters.SSLSTA1'; parameters.SSLSTA2']);
    else
        [delH1, ~] = computeDSLHeatUnary(P, y1val, T, PRef, TRef, parameters);
        [~, delH2] = computeDSLHeatUnary(P, y1, T, PRef, TRef, parameters);

        % [delH1,delH2] = computeQHeatUnary(P, y1, T, q1, q2, PRef, TRef, qRef, parameters);
    end

    if parameters.isResin
        delH1 = -parameters.delUb_1;
        delH2 = -parameters.delUb_2;
    end

    % Row 4: Energy balance
    nz=nz+1; ii(nz)=4; jj(nz)=2; vv(nz) = -coeff_q * (delH1 - cp_a * T_dim);
    nz=nz+1; ii(nz)=4; jj(nz)=3; vv(nz) = -coeff_q * (delH2 - cp_a * T_dim);
    nz=nz+1; ii(nz)=4; jj(nz)=4; vv(nz) = Ceff * TRef / tRef;
    nz=nz+1; ii(nz)=4; jj(nz)=6; vv(nz) = cp_g / R * PRef / tRef;
end

% Row 5: Wall energy balance
nz=nz+1; ii(nz)=5; jj(nz)=5; vv(nz) = 1;

% Row 6: Overall material balance
nz=nz+1; ii(nz)=6; jj(nz)=2; vv(nz) = R * T_dim * coeff_q;
nz=nz+1; ii(nz)=6; jj(nz)=3; vv(nz) = R * T_dim * coeff_q;
nz=nz+1; ii(nz)=6; jj(nz)=4; vv(nz) = -P_dim / (T * tRef);
nz=nz+1; ii(nz)=6; jj(nz)=6; vv(nz) = PRef / tRef;

M = sparse(ii(1:nz), jj(1:nz), vv(1:nz), 6, 6);

% Isothermal simplification
if parameters.isIsothermal
    M(1,4) = 0;                        % remove T coupling from y1 equation
    M(4,:) = sparse([1],[4],[1],1,6);  % identity row: dT/dt = 0
    M(5,:) = sparse([1],[5],[1],1,6);  % identity row: dTw/dt = 0
    M(6,4) = 0;                        % remove T coupling from P equation
end
end