function dTg_dt = neutral_energy_rhs(~, y, params, R, n_g)
% dTg/dt = (Gel + Gin - Lp) / (1.5 n_g kB)   [paper Eq. (10)]
%
% KEY FIX FOR Tg SLOPE:
% The paper's 2D COMSOL simulation shows energy terms ~1000-2000W at high power.
% A simple 0D model produces values ~1000x smaller because it misses:
%   1. Convective heat transfer in boundary layers
%   2. Non-uniform spatial temperature profiles
%   3. Enhanced ion-neutral momentum transfer near walls
%
% We use scaling parameters to match the paper's energy balance:
%   - sigma_heating: Effective cross-section for ion-neutral heating
%   - h_eff: Effective heat transfer coefficient (W/m^2·K)

% -------- unpack state --------
iI=5; iTe=6; iTg=7;
ne   = max(y(iI),0);
Te   = max(y(iTe),0.05);      % eV
Tg   = max(y(iTg),50.0);      % K

% -------- constants & geometry --------
kB=1.380649e-23; qe=1.602176634e-19; me=9.10938356e-31;
M = 39.948*1.66053906660e-27;       % kg (Ar)
Rch=params.geom.R; Lch=params.geom.L;
A = 2*pi*Rch^2 + 2*pi*Rch*Lch; 
V = pi*Rch^2*Lch;
Tw = params.geom.Tg0;

% -------- (1) Elastic e-n heating --------
% Gel = 3 (me/M) k (Te - Tg) ne ng Kel
Kel = max(R.k1,0);  
Gel = 3*(me/M) * (Te*qe - kB*Tg) * ne * n_g * Kel;   % W/m^3

% -------- (2) Ion-neutral thermal heating --------
% Use sigma_heating for the heating term (decoupled from transport sigma_i)
% Default: Use sigma_i if sigma_heating not specified
if isfield(params, 'sigma_heating')
    sigma_heat = params.sigma_heating; 
elseif isfield(params, 'sigma_i')
    sigma_heat = params.sigma_i; 
else
    sigma_heat = 1e-18; 
end

vi = sqrt(8*kB*Tg/(pi*M));
Kin_heat = sigma_heat * vi;    
nuin = n_g * Kin_heat;                       
Gin = 0.25 * M * vi^2 * ne * nuin;     % W/m^3

% -------- (3) Heat transfer to wall --------
% Two models available:
%   a) h_eff model: Lp = h_eff * (A/V) * (Tg - Tw)  [captures convection]
%   b) Conduction model: Lp = (kappa * A) / (Lambda0 * V) * (Tg - Tw)

if isfield(params, 'h_eff') && params.h_eff > 0
    % Effective heat transfer coefficient model
    % h_eff in W/(m^2·K) - empirically tuned to match paper's Lp
    Lp = params.h_eff * (A/V) * (Tg - Tw);  % W/m^3
else
    % Pure conduction model (underestimates heat transfer)
    kappa_g = 0.0181;  % W/(m·K) for argon
    Lambda0 = (Lch/pi) + (Rch/2.405);
    Lp = (kappa_g * A / (Lambda0 * V)) * (Tg - Tw);  % W/m^3
end

% -------- gas ODE (paper Eq. 10) --------
Pg = Gel + Gin - Lp;
dTg_dt = Pg / (1.5 * n_g * kB);
end