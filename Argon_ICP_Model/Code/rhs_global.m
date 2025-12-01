function dydt = rhs_global(t, y, params, ~)
% State: y = [ nAr  nArm  nArr  nArp  nI  Te  Tg ]
iAr=1; iArm=2; iArr=3; iArp=4; iI=5; iTe=6; iTg=7;

% ----------------- unpack & floors -----------------
nAr  = max(y(iAr),0);
nArm = max(y(iArm),0);
nArr = max(y(iArr),0);
nArp = max(y(iArp),0);
nI   = max(y(iI),0);
Te   = max(y(iTe),0.05);   % eV
Tg   = max(y(iTg),50);     % K
ne   = nI;                  % quasi-neutral

% ----------------- rate coefficients -----------------
R = RateCoefficients_Ar(params.geom.Tg0, params.p0, Te);
R = adapt_rates(R, params.ps);     % pass-through unless you set scales

% ----------------- geometry & constants --------------
qe  = 1.602176634e-19;
me  = 9.10938356e-31;
MAr = 39.948*1.66053906660e-27;
kB  = 1.380649e-23;

Rch   = params.geom.R;
Lch   = params.geom.L;
A_ends= 2*pi*Rch^2;
A_side= 2*pi*Rch*Lch;
A     = A_ends + A_side;
V     = pi*Rch^2*Lch;

% FIXED: background neutral density from gas law (at actual Tg, not Tg0)
%n_g = params.p0/(kB*Tg);   % CORRECTED: Using actual Tg
n_g = nAr + nArp + nArm + nArr;
% ----------------- transport to walls ----------------
% Ion Bohm speed and effective area via ion mfp with background gas
cs       = sqrt(qe*Te/MAr);
lambda_i = 1 / max(params.sigma_i * n_g, 1e-12);  % Ar+-Ar mfp
hL = 0.86*(3 + Lch/(2*lambda_i))^(-0.5);
hR = 0.80*(4 + Rch/(   lambda_i))^(-0.5);
Aeff = hL*A_ends + hR*A_side;

k_iwall    = (cs*Aeff)/V;                      % 1/s
Vs_eV = 0.5*Te*log(MAr/(2*pi*me));  % Sheath potential
n_s = nI * exp(-Vs_eV/Te);          % Boltzmann relation
L_ion_wall = k_iwall * nI;

% Neutral/excited wall-loss rates (paper Eqs. 8,9)
k_gwall = kn_wall_paper(Tg, params.geom, params.p0, params.gamma_g, params.sigma_nn);
k_mwall = kn_wall_paper(Tg, params.geom, params.p0, params.gamma_m, params.sigma_nn);
k_rwall = kn_wall_paper(Tg, params.geom, params.p0, params.gamma_r, params.sigma_nn);
k_pwall = kn_wall_paper(Tg, params.geom, params.p0, params.gamma_p, params.sigma_nn);

L_g_wall = k_gwall * nAr;
L_m_wall = k_mwall * nArm;
L_r_wall = k_rwall * nArr;
L_p_wall = k_pwall * nArp;

% ----------------- gas-flow closure ------------------
% Use Q0 (particles/s) from the paper's table
n_total_actual = nAr + nArm + nArr + nArp;  % Current total neutrals
%n_g = n_total_actual
n_ref = params.p0/(kB * Tg);
S_flow = (params.Q0 / V) * (1 - n_total_actual / n_ref);

% ----------------- volume reaction sources -----------
% Ionization (k2-k5)
S_ion_Ar  = ne*R.k2*nAr;% R2  Ar + e -> Ar+ + 2e  
S_ion_Arm = ne*R.k3*nArm; % R3 Arm + e -> Ar+ + 2e
S_ion_Arr = ne*R.k4*nArr; % R4 Arr + e -> Ar+ + 2e
S_ion_Arp = ne*R.k5*nArp; % R5 Arp + e -> Ar+ + 2e

% Excitation from ground (k6-k8)
S_exc_gm = ne*R.k6*nAr; % Ar + e -> Arm + e    
S_exc_gr  = ne*R.k7*nAr; % Ar + e -> Arr + e  
S_exc_gp = ne*R.k8*nAr; % Ar + e -> Arp + e

% De-excitation to ground (k9-k11)
S_dex_m  = ne*R.k9*nArm; % Arm + e -> Ar + e    
S_dex_r   = ne*R.k10*nArr; % Arr + e -> Ar + e
S_dex_p  = ne*R.k11*nArp; % Arp + e -> Ar + e

% Excited-to-excited transfers (k12-k17)
S_m_to_r = ne*R.k12*nArm; % Arm + e -> Arr + e   
S_m_to_p  = ne*R.k13*nArm; % Arm + e -> Arp + e
S_r_to_p = ne*R.k14*nArr;  % Arm + e -> Arp + e 
S_r_to_m  = ne*R.k15*nArr; % Arr + e -> Arm + e
S_p_to_m = ne*R.k16*nArp;   % Arp + e -> Arm + e
S_p_to_r  = ne*R.k17*nArp; % Arp + e -> Arr + e

% Radiative (k18-k20)
S_r_rad     = R.k18*nArr;   % Arr -> Ar + hv
S_p_to_m_hv = R.k19*nArp;   % Arp -> Arm + hv  
S_p_to_r_hv = R.k20*nArp;   % Arp -> Arr + hv

% ----------------- species ODEs ----------------------
dnAr_dt  = +S_dex_m +S_dex_r +S_dex_p +S_r_rad +L_ion_wall ...
           -(S_exc_gm +S_exc_gr +S_exc_gp) -S_ion_Ar -L_g_wall... 
            + S_flow;

%dnAr_dt  = +S_dex_m +S_dex_r +S_dex_p +S_r_rad +L_ion_wall ...
 %          -(S_exc_gm +S_exc_gr +S_exc_gp) -S_ion_Ar -L_g_wall

dnArm_dt = +S_exc_gm +S_r_to_m +S_p_to_m +S_p_to_m_hv ...
           -(S_m_to_r +S_m_to_p) -S_dex_m -S_ion_Arm -L_m_wall;

dnArr_dt = +S_exc_gr +S_m_to_r +S_p_to_r +S_p_to_r_hv ...
           -(S_r_to_p +S_r_to_m) -S_dex_r -S_ion_Arr -S_r_rad -L_r_wall;

dnArp_dt = +S_exc_gp +S_m_to_p +S_r_to_p ...
           -(S_p_to_m +S_p_to_r) -S_dex_p -S_ion_Arp ...
           -(S_p_to_m_hv +S_p_to_r_hv) -L_p_wall;

dnI_dt   = +(S_ion_Ar +S_ion_Arm +S_ion_Arr +S_ion_Arp) - L_ion_wall;

% ----------------- energy ODEs -----------------------
dTe_dt_raw = electron_energy_rhs(t, y, params, R);
dTe_dt     = dTe_dt_raw - (Te/max(ne,1e-30))*dnI_dt;   % ∂(3/2 ne Te)/∂t coupling

dTg_dt = neutral_energy_rhs(t, y, params, R) - (Tg/n_g) * (dnAr_dt + dnArm_dt + dnArr_dt + dnArp_dt);

% ----------------- pack ------------------------------
dydt        = zeros(size(y));
dydt(iAr)   = dnAr_dt;
dydt(iArm)  = dnArm_dt;
dydt(iArr)  = dnArr_dt;
dydt(iArp)  = dnArp_dt;
dydt(iI)    = dnI_dt;
dydt(iTe)   = dTe_dt;
dydt(iTg)   = dTg_dt;
end