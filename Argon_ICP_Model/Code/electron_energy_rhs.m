function dTe_dt = electron_energy_rhs(~, y, params, R)
% y = [nAr nArm nArr nArp nI Te Tg]; Te in eV -> dTe/dt (eV/s)

% unpack / floors
nAr  = max(y(1),0);  
nArm = max(y(2),0);  
nArr = max(y(3),0);  
nArp = max(y(4),0);
ne   = max(y(5),1e10);  
Te = max(y(6),0.05);  
Tg = max(y(7),50);
nI   = ne;  % ion density (quasi-neutral)

% constants & geom
qe=1.602176634e-19; 
kB=1.380649e-23; 
me=9.10938356e-31; 
MAr=39.948*1.66053906660e-27;
Rch=params.geom.R; 
Lch=params.geom.L; 
A_ends=2*pi*Rch^2; 
A_side=2*pi*Rch*Lch; 
A=A_ends+A_side;
V=pi*Rch^2*Lch;

% CORRECTED: Calculate absorbed power properly (PRF - Pcoil)
% First calculate Rind and Rcoil to find Pcoil
n_g = params.p0/(kB*Tg);  % Use Tg0 (300K) not current Tg
%n_g = nAr + nArp + nArm + nArr;
Kel = max(R.k1, 0);
nu_m = n_g * Kel;  % electron-neutral collision frequency

% Calculate Rind using the separate function
Rind = calc_rind(params, ne, nu_m);

% Get Rcoil from parameters
Rcoil = 2;  % From Table 2 in the paper
if isfield(params, 'Rcoil')
    Rcoil = params.Rcoil;
end



% Calculate current from PRF = 0.5*(Rind + Rcoil)*I^2
PRF = max(params.PRF, 0);
Rtotal = max(Rind + Rcoil, 1e-9);
I_squared = 2 * PRF / Rtotal;

% Calculate Pcoil = 0.5 * Rcoil * I^2
Pcoil = 0.5 * Rcoil * I_squared;

% CORRECTED: Pabs = (PRF - Pcoil) / V
Pabs_total = max(PRF - Pcoil, 0);  % Total absorbed power in W
Pabs = Pabs_total / max(V, 1e-12);  % Power density in W/m^3

% Ion mean free path and effective area factors
lambda_i = 1 / max(params.sigma_i*n_g, 1e-12);

% Effective area factors (hL and hR from paper)
hL = 0.86*(3 + Lch/(2*lambda_i))^(-0.5);  
hR = 0.80*(4 + Rch/lambda_i)^(-0.5);
Aeff = hL*A_ends + hR*A_side;

% Bohm velocity and sheath potential
uB = sqrt(qe*Te/MAr);  % Bohm velocity (same as cs)
Vs_eV = 0.5*Te*log(MAr/(2*pi*me));  % sheath potential (eV)

% Pew uses electron density at sheath (n_es)  
n_es = nI * exp(-Vs_eV/Te);          % Boltzmann relation
n_is = n_es;  % For electropositive plasma
epsilon_ew = 2.0*Te;  % energy per electron lost (eV)
Pew = (A/V)*n_es*uB*(epsilon_ew*qe);  % electron wall power density

% Ion wall power loss
epsilon_iw = Vs_eV + 0.5*Te;  % energy per ion lost (eV)
Piw = (A/V)*n_is*uB*(epsilon_iw*qe);  % ion wall power density


% Energy levels
eps_2 = 15.76; 
eps_3 = 4.43;
eps_4 = 3.96;
eps_5 = 2.26;
eps_6 = 11.5; 
eps_7 = 11.6; 
eps_8 = 12.9;
eps_12 = 1.6; 
eps_13 = 1.2; 
eps_14 = 1.1;

% CORRECTED: Calculate epsilon_c using equation (13) from the paper
% epsilon_c^(X) = eps_iz + sum_l(k_exl^(X)/k_iz) + (3*me*kel*Te)/(M*k_iz)
% For each ionizing species X, calculate its epsilon_c

% Ground state Ar ionization with associated excitations
if R.k2 > 1e-30  % Avoid division by zero
    % Excitation to loss ratio for ground state
    excit_ratio_g = (R.k6*eps_6 + R.k7*eps_7 + R.k8*eps_8)/R.k2;
    % Elastic energy transfer per ionization
    elastic_ratio_g = (3*me/MAr)*Kel*Te/R.k2;
    epsilon_c_g = eps_2 + excit_ratio_g + elastic_ratio_g;
else
    epsilon_c_g = eps_2;
end

% Metastable Ar* ionization
if R.k3 > 1e-30
    % For metastable, consider excitations from metastable state
    excit_ratio_m = (R.k12*eps_12 + R.k13*eps_13)/R.k3;
    elastic_ratio_m = (3*me/MAr)*Kel*Te/R.k3;
    epsilon_c_m = eps_3 + excit_ratio_m + elastic_ratio_m;
else
    epsilon_c_m = eps_3;
end

% Resonant Ar^r ionization
if R.k4 > 1e-30
    excit_ratio_r = (R.k14*eps_14)/R.k4;
    elastic_ratio_r = (3*me/MAr)*Kel*Te/R.k4;
    epsilon_c_r = eps_4 + excit_ratio_r + elastic_ratio_r;
else
    epsilon_c_r = eps_4;
end

% 4p state Ar^p ionization (no further excitations from this level)
if R.k5 > 1e-30
    elastic_ratio_p = (3*me/MAr)*Kel*Te/R.k5;
    epsilon_c_p = eps_5 + elastic_ratio_p;
else
    epsilon_c_p = eps_5;
end

% Power lost due to ionization from each state using its specific epsilon_c
P_ion = ne * qe * (R.k2*nAr*epsilon_c_g + R.k3*nArm*epsilon_c_m + ...
                   R.k4*nArr*epsilon_c_r + R.k5*nArp*epsilon_c_p);
%P_ion = ne * qe * R.k2*nAr*epsilon_c_g
TeJ = Te * qe;
% CORRECTED: Pev now only includes P_ion (which contains all losses via epsilon_c formula)
% and de-excitation returns (energy gains)
% Pev = P_ion + P_dex_g + P_dex_x;
Pev = P_ion;
% Continuous elastic cooling (required for energy conservation with gas equation)
%P_el = 3*(me/MAr) * ne * n_g * Kel * (TeJ - kB*Tg);

% Total 
%Pev = P_ion + P_el;  % Even though paper shows only P_ion

% total electron power balance -> Te ODE
Ploss = Piw + Pew + Pev;
denom = 1.5*ne*qe;  % J m^-3 per eV
dTe_dt = (Pabs - Ploss)/max(denom,1e-30);
end