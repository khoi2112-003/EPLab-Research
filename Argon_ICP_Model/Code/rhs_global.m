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
ne   = nI;                 % quasi-neutral

% Total neutral density (sum of all neutral states)
n_g = nAr + nArm + nArr + nArp;

% ----------------- rate coefficients -----------------
R = RateCoefficients_Ar(params.geom.Tg0, params.p0, Te);
R = adapt_rates(R, params.ps); 

% ----------------- geometry & constants --------------
qe  = 1.602176634e-19;
me  = 9.10938356e-31;
MAr = 39.948*1.66053906660e-27;
kB  = 1.380649e-23;

Rch   = params.geom.R;
Lch   = params.geom.L;
A_ends= 2*pi*Rch^2;
A_side= 2*pi*Rch*Lch;
V     = pi*Rch^2*Lch;
A_grid = params.beta_g * pi * Rch^2; % Open area for neutral loss

% ----------------- Gas Flow Balance ------------------
% CRITICAL FIX: Calculate the Q0 required to sustain p0 at T_g0.
% This overrides the hardcoded Q0 which may be inconsistent with the pressure sweep.
% At steady state (no plasma): Q0 = (1/4) * n_g0 * v_g0 * A_grid
n_ref = params.p0 / (kB * params.geom.Tg0);
v_ref = sqrt(8 * kB * params.geom.Tg0 / (pi * MAr));
S_in  = (0.25 * n_ref * v_ref * A_grid) / V;
%S_in = params.Q0 /V;
% Outflow (Effusion): scales with current n_g and sqrt(Tg)
vg = sqrt(8*kB*Tg/(pi*MAr)); 
S_out_factor = (0.25 * vg * A_grid) / V;

% ----------------- transport to walls ----------------
% Ion Bohm speed and effective area 
cs       = sqrt(qe*Te/MAr);
lambda_i = 1 / max(params.sigma_i * n_g, 1e-12);  
hL = 0.86*(3 + Lch/(2*lambda_i))^(-0.5);
hR = 0.80*(4 + Rch/(   lambda_i))^(-0.5);
Aeff = hL*A_ends + hR*A_side;

k_iwall    = (cs*Aeff)/V;                 
Vs_eV = 0.5*Te*log(MAr/(2*pi*me));  
L_ion_wall = k_iwall * nI;

% Neutral/excited wall-loss rates (Diffusive)
k_gwall = kn_wall_paper(Tg, params.geom, n_g, params.gamma_g, params.sigma_nn);
k_mwall = kn_wall_paper(Tg, params.geom, n_g, params.gamma_m, params.sigma_nn);
k_rwall = kn_wall_paper(Tg, params.geom, n_g, params.gamma_r, params.sigma_nn);
k_pwall = kn_wall_paper(Tg, params.geom, n_g, params.gamma_p, params.sigma_nn);

L_g_wall = k_gwall * nAr;
L_m_wall = k_mwall * nArm;
L_r_wall = k_rwall * nArr;
L_p_wall = k_pwall * nArp;

% ----------------- volume reaction sources -----------
% Ionization
S_ion_Ar  = ne*R.k2*nAr;
S_ion_Arm = ne*R.k3*nArm;
S_ion_Arr = ne*R.k4*nArr;
S_ion_Arp = ne*R.k5*nArp;

% Excitation 
S_exc_gm = ne*R.k6*nAr;     
S_exc_gr = ne*R.k7*nAr;   
S_exc_gp = ne*R.k8*nAr; 

% De-excitation 
S_dex_m  = ne*R.k9*nArm;    
S_dex_r  = ne*R.k10*nArr; 
S_dex_p  = ne*R.k11*nArp; 

% Excited transfers
S_m_to_r = ne*R.k12*nArm;    
S_m_to_p = ne*R.k13*nArm; 
S_r_to_p = ne*R.k14*nArr;   
S_r_to_m = ne*R.k15*nArr; 
S_p_to_m = ne*R.k16*nArp;   
S_p_to_r = ne*R.k17*nArp; 

% Radiative
S_r_rad     = R.k18*nArr;   
S_p_to_m_hv = R.k19*nArp;    
S_p_to_r_hv = R.k20*nArp;   

% ----------------- species ODEs ----------------------
% Note: S_in is added entirely to the Ground State equation.
dnAr_dt  = +S_dex_m +S_dex_r +S_dex_p +S_r_rad +L_ion_wall ...
           -(S_exc_gm +S_exc_gr +S_exc_gp) -S_ion_Ar -L_g_wall ...
           + S_in - (nAr * S_out_factor);

dnArm_dt = +S_exc_gm +S_r_to_m +S_p_to_m +S_p_to_m_hv ...
           -(S_m_to_r +S_m_to_p) -S_dex_m -S_ion_Arm -L_m_wall ...
           - (nArm * S_out_factor);

dnArr_dt = +S_exc_gr +S_m_to_r +S_p_to_r +S_p_to_r_hv ...
           -(S_r_to_p +S_r_to_m) -S_dex_r -S_ion_Arr -S_r_rad -L_r_wall ...
           - (nArr * S_out_factor);

dnArp_dt = +S_exc_gp +S_m_to_p +S_r_to_p ...
           -(S_p_to_m +S_p_to_r) -S_dex_p -S_ion_Arp ...
           -(S_p_to_m_hv +S_p_to_r_hv) -L_p_wall ...
           - (nArp * S_out_factor);

dnI_dt   = +(S_ion_Ar +S_ion_Arm +S_ion_Arr +S_ion_Arp) - L_ion_wall;

% ----------------- energy ODEs -----------------------
dTe_dt_raw = electron_energy_rhs(t, y, params, R, n_g);
dTe_dt     = dTe_dt_raw - (Te/max(ne,1e-30))*dnI_dt; 

dTg_dt = neutral_energy_rhs(t, y, params, R, n_g);

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