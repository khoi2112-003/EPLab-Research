function dTg_dt = neutral_energy_rhs(~, y, params, R)
% dTg/dt = (Gel + Gin - Lp) / (1.5 n_g kB)   [paper Eq. (10)]
% Also computes Li (neutral enthalpy loss due to ionization) for diagnostics.

% -------- unpack state --------
iAr=1; iArm=2; iArr=3; iArp=4; iI=5; iTe=6; iTg=7;
nAr  = max(y(iAr),0);
nArm = max(y(iArm),0);
nArr = max(y(iArr),0);
nArp = max(y(iArp),0);
ne   = max(y(iI),0);
Te   = max(y(iTe),0.05);      % eV
Tg   = max(y(iTg),50.0);      % K

% -------- constants & geometry --------
kB=1.380649e-23; qe=1.602176634e-19; me=9.10938356e-31;
M = 39.948*1.66053906660e-27;       % kg (Ar)
Rch=params.geom.R; Lch=params.geom.L;
A = 2*pi*Rch^2 + 2*pi*Rch*Lch; 
V = pi*Rch^2*Lch;

% FIXED: Use actual gas temperature Tg for neutral density
% At constant pressure: n_g = p0/(kB*Tg)
%n_g = params.p0/(kB*Tg);   % CORRECTED: Using actual Tg, not Tg0
n_g = nAr + nArp + nArm + nArr;
% (1) elastic e–n heating  (equal & opposite to electron elastic cooling)
Kel = max(R.k1,0);                      % m^3/s  (k1 in Table)
Gel = 3*(me/M) * (Te*qe - kB*Tg) * ne * n_g * Kel;   % W/m^3

% (2) ion–neutral thermal heating (Ti = Tg)
vi   = sqrt(8*kB*Tg/(pi*M));
Kin  = params.sigma_i * vi;             % m^3/s (Ar+–Ar)
nuin = n_g * Kin;                       % 1/s
Gin  = 0.25 * M * vi^2 * ne * nuin;     % W/m^3

% (3) heat conduction to wall
kappa_g = 0.0181;                        % W m^-1 K^-1
Lambda0  = (Lch/pi) + (Rch/2.405);
Tw       = params.geom.Tg0;
Lp       = (kappa_g * A / (Lambda0 * V)) * (Tg - Tw);  % W/m^3
Li = 1.5 * kB * Tg * ne * n_g * Kin;

% -------- gas ODE (paper Eq. 10) --------
Pg     = Gel + Gin - Lp;% - Li;                 % do NOT subtract Li here

% Add before the final line of neutral_energy_rhs.m
if Tg > 340 && Tg < 342  % Only print near final state
    fprintf('\n=== INSIDE neutral_energy_rhs at Tg=%.1fK ===\n', Tg);
    fprintf('  ne = %.3e m^-3\n', ne);
    fprintf('  n_g = %.3e m^-3\n', n_g);
    fprintf('  Gel = %.3f W/m^3 (total: %.3f W)\n', Gel, Gel*V);
    fprintf('  Gin = %.3f W/m^3 (total: %.3f W)\n', Gin, Gin*V);
    fprintf('  Lp  = %.3f W/m^3 (total: %.3f W)\n', Lp, Lp*V);
    fprintf('  Net = %.3f W/m^3 (total: %.3f W)\n', Gel+Gin-Lp, (Gel+Gin-Lp)*V);
    fprintf('  dTg/dt = %.3f K/s\n', (Gel+Gin-Lp)/(1.5*n_g*kB));
end

dTg_dt = Pg / (1.5 * n_g * kB);
end