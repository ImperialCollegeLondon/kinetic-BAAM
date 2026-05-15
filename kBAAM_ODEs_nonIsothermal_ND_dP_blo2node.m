%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Imperial College London, United Kingdom
% Multiphase Systems Laboratory
% Year:     2026
% MATLAB:   R2024a
% Authors:  Hassan Azzan (HA)
%
% Purpose:
% Evaluates the RHS of the dimensionless ODEs for the blowdown step of the non-isothermal
% k-BAAM using TWO spatial nodes. P, T, and Tw are shared between the two nodes. Only
% composition (y1) and adsorbed amounts (q1, q2) differ between nodes.
%
% State vector X (10 states, all dimensionless):
%   X(1)  = y1_1       mole fraction of CO2 in node 1 (feed end)   [-]
%   X(2)  = q1_1/qRef  adsorbed CO2 in node 1                       [-]
%   X(3)  = q2_1/qRef  adsorbed N2  in node 1                       [-]
%   X(4)  = y1_2       mole fraction of CO2 in node 2 (product end) [-]
%   X(5)  = q1_2/qRef  adsorbed CO2 in node 2                       [-]
%   X(6)  = q2_2/qRef  adsorbed N2  in node 2                       [-]
%   X(7)  = T/TRef     shared gas/solid temperature                  [-]
%   X(8)  = Tw/TwRef   shared wall temperature                       [-]
%   X(9)  = P1/PRef    node 1 (feed end) pressure                    [-]
%   X(10) = P2/PRef    node 2 (product end) pressure                 [-]
%
% Flow model (nodes sized by adsorption loading fraction, independent pressures):
%   Node 1 (feed end)    occupies fraction f = loadingFraction of column volume.
%   Node 2 (product end) occupies fraction (1-f) of column volume.
%   Darcy velocity node 1 → node 2 (node 1 centre to node 2 centre = L/2):
%     v_12  = (2/L) * darcyK * (P1 - P2)     * PRef   [m/s]
%   Darcy velocity node 2 → outlet (node 2 centre to outlet = L2/2):
%     v_out = (2/L2) * darcyK * (P2 - P_blo)  * PRef   [m/s]
%
% Initial conditions (set in Outputs before blowdown):
%   Node 1: at equilibrium with y1_in at (P_ads_end, T_ads_end)
%   Node 2: same y1, q1, q2 as the start of adsorption (from previous cycle)
%   Shared T, Tw, P: taken directly from end of adsorption step
%
% Last modified:
% - 2026-05-01, HA: Initial creation
%
% Input arguments:
%   - t:          dimensionless time [-]
%   - X:          10x1 vector of dimensionless states
%   - parameters: struct of process/adsorbent parameters
%
% Output arguments:
%   - dXdt: 10x1 vector of time derivatives
%
% Dependencies:
%   - DSL.m / SSLSTA.m
%   - LDFCoefficient.m
%   - computeDSLHeatUnary.m / computeSSLSTAHeatBinaryBT.m
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function dXdt = kBAAM_ODEs_nonIsothermal_ND_dP_blo2node(t, X, parameters,stepName)

R = 8.3145;

    % ---- Unpack states ----
    y1_1 = max(0, min(1, X(1)));
    q1_1 = X(2);
    q2_1 = X(3);
    y1_2 = max(0, min(1, X(4)));
    q1_2 = X(5);
    q2_2 = X(6);
    T    = max(X(7),  1e-9);   % shared dimensionless temperature
    Tw   = X(8);               % shared dimensionless wall temperature
    P1   = max(X(9),  1e-9);   % node 1 dimensionless pressure
    P2   = max(X(10), 1e-9);   % node 2 dimensionless pressure

if string(stepName) == "blo"
    % Per-node volumes based on adsorption loading fraction.
    % Clamp to [0.1, 0.95] to prevent degenerate small nodes that cause excessive stiffness.
    f_node = max(0.05, min(0.95, parameters.loadingFraction));
else
    % ---- Unpack states ----
    % y1_2 = max(0, min(1, X(1)));
    % q1_2 = X(2);
    % q2_2 = X(3);
    % y1_1 = max(0, min(1, X(4)));
    % q1_1 = X(5);
    % q2_1 = X(6);
    % T    = max(X(7),  1e-9);   % shared dimensionless temperature
    % Tw   = X(8);               % shared dimensionless wall temperature
    % P2   = max(X(9),  1e-9);   % node 1 dimensionless pressure
    % P1   = max(X(10), 1e-9);   % node 2 dimensionless pressure

    % Per-node volumes based on adsorption loading fraction.
    % Clamp to [0.1, 0.95] to prevent degenerate small nodes that cause excessive stiffness.
    f_node = 1-max(0.05, min(0.95, parameters.loadingFraction));
end
% ---- Unpack parameters ----
tRef      = parameters.timeRef;
qRef      = parameters.qRef;
TRef      = parameters.TRef;
TwRef     = parameters.TwRef;
PRef      = parameters.PRef;
Ab        = parameters.Ab;
coeff_q   = parameters.Ab_rhos_qRef_tRef;   % Ab*rho_s*qRef/tRef
cp_a      = parameters.cp_a;
cp_g      = parameters.cp_g;
e         = parameters.e_bed;
V         = parameters.V_column;
L         = parameters.L;
A         = parameters.A_in;
tRef_qRef = parameters.tRef_qRef;
darcyK    = parameters.darcyK;
% % f_node = parameters.loadingFraction;
V_node1 = f_node * V;
V_node2 = (1 - f_node) * V;

% ---- Dimensional scalings ----
P1_dim = P1 * PRef;
P2_dim = P2 * PRef;
T_dim  = T  * TRef;

if string(stepName) == "blo"
% ---- Prescribed outlet pressure (product end of node 2) ----
P_out = parameters.P_blo(t * tRef) / PRef;
else
% ---- Prescribed outlet pressure (product end of node 2) ----
P_out = parameters.P_evac(t * tRef) / PRef;
end

% ---- Darcy flows: independent per leg, each scaled by its own node length ----
L1 = f_node * L;        % length of node 1 [m]
L2 = (1 - f_node) * L; % length of node 2 [m]
v_12  = (2./L) * darcyK * (P1 - P2)    * PRef;   % [m/s] node 1→2, centre-to-interface = L1/2
v_out = (2/L2) * darcyK * (P2 - P_out) * PRef;   % [m/s] node 2→outlet, centre-to-outlet = L2/2

% Inter-node molar flow (upstream density from P1)
F_12  = P1_dim * A * e / (R * T_dim) * v_12;

% Outlet molar flow (upstream density from P2)
F_out = P2_dim * A * e / (R * T_dim) * v_out;

% No inlet (blowdown: feed valve closed)
F_in  = 0;

% ---- Axial dispersion between nodes ----
DL = 0.7 * parameters.Dm + abs(v_12) * parameters.rp;  % [m2/s]
dz = L1/2 + L2/2;  % distance between node centres = L/2 (independent of node sizing)
P_int_dim = (P1_dim + P2_dim) / 2;  % interface pressure for dispersive density
% F_disp = DL * A * e * P_int_dim / (R * T_dim) * (y1_1 - y1_2) / dz;
F_disp = 0;

% ---- LDF kinetics: node 1 (at P1) ----
[q1_star1, q2_star1] = getEq(P1, y1_1, T, PRef, TRef, parameters);
[k1_1, k2_1] = LDFCoefficient(P1_dim, y1_1, T_dim, q1_star1, q2_star1, parameters);
dq1dt_1 = tRef_qRef * k1_1 * (q1_star1 - q1_1 * qRef);
dq2dt_1 = tRef_qRef * k2_1 * (q2_star1 - q2_1 * qRef);

% ---- LDF kinetics: node 2 (at P2) ----
[q1_star2, q2_star2] = getEq(P2, y1_2, T, PRef, TRef, parameters);
[k1_2, k2_2] = LDFCoefficient(P2_dim, y1_2, T_dim, q1_star2, q2_star2, parameters);
dq1dt_2 = tRef_qRef * k1_2 * (q1_star2 - q1_2 * qRef);
dq2dt_2 = tRef_qRef * k2_2 * (q2_star2 - q2_2 * qRef);

% ---- Assemble RHS ----
f = zeros(10, 1);

% Node 1: y1 balance (F_in=0, convective outflow F_12, dispersive loss F_disp)
f(1) = R * T_dim / P1_dim * (0 - y1_1 * F_12 - F_disp) / (e * V_node1);

% Node 1: adsorbed phase
f(2) = dq1dt_1;
f(3) = dq2dt_1;

% Node 2: y1 balance (convective inflow F_12 at y1_1, dispersive gain F_disp, outflow F_out at y1_2)
f(4) = R * T_dim / P2_dim * (y1_1 * F_12 + F_disp - y1_2 * F_out) / (e * V_node2);

% Node 2: adsorbed phase
f(5) = dq1dt_2;
f(6) = dq2dt_2;

% Shared energy balance (whole column, F_out exits at node-2 T)
cpg_eV_total = cp_g / (V * e);
two_hin = parameters.two_hin_rin_e;
Qheat = 0;
if ~parameters.isIsothermal
    f(7) = cpg_eV_total * (0 - F_out * T_dim) ...
         - two_hin * (T_dim - Tw * TwRef);
end

% Shared wall energy balance
if ~parameters.isIsothermal
    f(8) = parameters.wall_prefactor * ...
        ( parameters.wall_coeff1 * (T_dim - Tw * TwRef) ...
         -parameters.wall_coeff2 * (Tw * TwRef - TRef) ...
         + Qheat);
end

% Per-node overall material balance
f(9)  = R * T_dim * (F_in - F_12)  / (e * V_node1);   % node 1
f(10) = R * T_dim * (F_12 - F_out) / (e * V_node2);   % node 2

% ---- Build mass matrix and solve ----
M = buildMassMatrix10(y1_1, q1_1, q2_1, y1_2, q1_2, q2_2, T, Tw, P1, P2, parameters, R, tRef, qRef, TRef, TwRef, PRef, Ab, coeff_q, cp_a, cp_g, f_node, V_node1, V_node2);
dXdt = M \ f;
end

% ============================================================
function [q1s, q2s] = getEq(P, y1, T, PRef, TRef, parameters)
if parameters.SSLSTA
    [q1s, q2s] = SSLSTA(P * PRef, y1, T * TRef, parameters);
else
    [q1s, q2s] = DSL(P * PRef, y1, T * TRef, ...
        parameters.qsb_1, parameters.qsd_1, parameters.qsb_2, parameters.qsd_2, ...
        parameters.bo_1,  parameters.do_1,  parameters.bo_2,  parameters.do_2,  ...
        parameters.delUb_1, parameters.delUd_1, parameters.delUb_2, parameters.delUd_2);
end
end

% ============================================================
function M = buildMassMatrix10(y1_1, q1_1, q2_1, y1_2, q1_2, q2_2, T, Tw, P1, P2, parameters, R, tRef, qRef, TRef, TwRef, PRef, Ab, coeff_q, cp_a, cp_g, f_node, V_node1, V_node2)

P1_dim = P1 * PRef;
P2_dim = P2 * PRef;
T_dim  = T  * TRef;
RT_PP1 = R * T_dim / P1_dim;
RT_PP2 = R * T_dim / P2_dim;

ii = zeros(40, 1); jj = zeros(40, 1); vv = zeros(40, 1);
nz = 0;

% ---- Row 1: y1_1 species balance (coupled to P1) ----
nz=nz+1; ii(nz)=1; jj(nz)=1;  vv(nz) = 1/tRef;
nz=nz+1; ii(nz)=1; jj(nz)=2;  vv(nz) = RT_PP1 * coeff_q;
nz=nz+1; ii(nz)=1; jj(nz)=7;  vv(nz) = -y1_1 / (T * tRef);
nz=nz+1; ii(nz)=1; jj(nz)=9;  vv(nz) =  y1_1 / (P1 * tRef);

% ---- Row 2: q1_1 (LDF) ----
nz=nz+1; ii(nz)=2; jj(nz)=2;  vv(nz) = 1;

% ---- Row 3: q2_1 (LDF) ----
nz=nz+1; ii(nz)=3; jj(nz)=3;  vv(nz) = 1;

% ---- Row 4: y1_2 species balance (coupled to P2) ----
nz=nz+1; ii(nz)=4; jj(nz)=4;  vv(nz) = 1/tRef;
nz=nz+1; ii(nz)=4; jj(nz)=5;  vv(nz) = RT_PP2 * coeff_q;
nz=nz+1; ii(nz)=4; jj(nz)=7;  vv(nz) = -y1_2 / (T * tRef);
nz=nz+1; ii(nz)=4; jj(nz)=10; vv(nz) =  y1_2 / (P2 * tRef);

% ---- Row 5: q1_2 (LDF) ----
nz=nz+1; ii(nz)=5; jj(nz)=5;  vv(nz) = 1;

% ---- Row 6: q2_2 (LDF) ----
nz=nz+1; ii(nz)=6; jj(nz)=6;  vv(nz) = 1;

if ~parameters.isIsothermal
    % Volume-weighted averages for heat capacity and delH
    q1_avg = f_node * q1_1 + (1 - f_node) * q1_2;
    q2_avg = f_node * q2_1 + (1 - f_node) * q2_2;
    Ceff = Ab * (parameters.rho_s * parameters.cp_s + cp_a * parameters.rho_s * qRef * (q1_avg + q2_avg));
    y1_avg = f_node * y1_1 + (1 - f_node) * y1_2;
    P_avg  = f_node * P1 + (1 - f_node) * P2;  % avg for delH evaluation
    if parameters.SSLSTA
        [delH1_1, delH2_1] = computeSSLSTAHeatBinaryBT(P_avg*PRef/1e5, y1_avg, T_dim, [parameters.SSLSTA1'; parameters.SSLSTA2']);
        delH1_2 = delH1_1;
        delH2_2 = delH2_1;
    else
        [delH1_1, delH2_1] = computeDSLHeatUnary(P1, y1_1, T, PRef, TRef, parameters);
        [delH1_2, delH2_2] = computeDSLHeatUnary(P2, y1_2, T, PRef, TRef, parameters);
    end
    if parameters.isResin
        delH1_1 = -parameters.delUb_1;
        delH2_1 = -parameters.delUb_2;
        delH1_2 = delH1_1;
        delH2_2 = delH2_1;
    end

    % ---- Row 7: shared energy balance ----
    nz=nz+1; ii(nz)=7; jj(nz)=2;  vv(nz) = -coeff_q*f_node     * (delH1_1 - cp_a * T_dim);  % node 1 q1
    nz=nz+1; ii(nz)=7; jj(nz)=3;  vv(nz) = -coeff_q*f_node     * (delH2_1 - cp_a * T_dim);  % node 1 q2
    nz=nz+1; ii(nz)=7; jj(nz)=5;  vv(nz) = -coeff_q*(1-f_node) * (delH1_2 - cp_a * T_dim);  % node 2 q1
    nz=nz+1; ii(nz)=7; jj(nz)=6;  vv(nz) = -coeff_q*(1-f_node) * (delH2_2 - cp_a * T_dim);  % node 2 q2
    nz=nz+1; ii(nz)=7; jj(nz)=7;  vv(nz) = Ceff * TRef / tRef;
    nz=nz+1; ii(nz)=7; jj(nz)=9;  vv(nz) = f_node     * cp_g / R * PRef / tRef;  % node 1 P contribution
    nz=nz+1; ii(nz)=7; jj(nz)=10; vv(nz) = (1-f_node) * cp_g / R * PRef / tRef;  % node 2 P contribution
end

% ---- Row 8: wall energy balance ----
nz=nz+1; ii(nz)=8; jj(nz)=8;  vv(nz) = 1;

% ---- Row 9: node 1 overall material balance ----
nz=nz+1; ii(nz)=9; jj(nz)=2;  vv(nz) = R * T_dim * coeff_q;   % q1_1
nz=nz+1; ii(nz)=9; jj(nz)=3;  vv(nz) = R * T_dim * coeff_q;   % q2_1
nz=nz+1; ii(nz)=9; jj(nz)=7;  vv(nz) = -P1_dim / (T * tRef);  % T coupling
nz=nz+1; ii(nz)=9; jj(nz)=9;  vv(nz) = PRef / tRef;           % P1

% ---- Row 10: node 2 overall material balance ----
nz=nz+1; ii(nz)=10; jj(nz)=5;  vv(nz) = R * T_dim * coeff_q;  % q1_2
nz=nz+1; ii(nz)=10; jj(nz)=6;  vv(nz) = R * T_dim * coeff_q;  % q2_2
nz=nz+1; ii(nz)=10; jj(nz)=7;  vv(nz) = -P2_dim / (T * tRef); % T coupling
nz=nz+1; ii(nz)=10; jj(nz)=10; vv(nz) = PRef / tRef;           % P2

M = sparse(ii(1:nz), jj(1:nz), vv(1:nz), 10, 10);

% Isothermal simplification
if parameters.isIsothermal
    M(1,7)  = 0;
    M(4,7)  = 0;
    M(7,:)  = sparse([1],[7],[1],1,10);    % dT/dt = 0
    M(8,:)  = sparse([1],[8],[1],1,10);    % dTw/dt = 0
    M(9,7)  = 0;
    M(10,7) = 0;
end
end
