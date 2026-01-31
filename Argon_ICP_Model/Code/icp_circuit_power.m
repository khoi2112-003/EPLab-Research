function Pabs = icp_circuit_power_improved(params, ne, Te_eV, Tg_K, Rrates, geom)
% IMPROVED power absorption with power-dependent coupling efficiency
%
% FIXES: At low power, Ar+ and Arp are overestimated because the model
% assumes 100% absorption. In reality, coupling efficiency varies with:
%   - Plasma density (low at low power)
%   - Impedance matching
%   - Skin depth effects
%
% This adds empirical power-dependent absorption: η(P)
%   Low power (< 300W):  40-80% absorption
%   High power (> 300W): 80-95% absorption
%
% This matches the observed trend: Ar+/Arp errors decrease with power!

% ------------------ common constants ------------------
eps0 = 8.8541878128e-12;
mu0  = 4*pi*1e-7;
c0   = 299792458;
me   = 9.10938356e-31;
qe   = 1.602176634e-19;
kB   = 1.380649e-23;

omega = params.omega;

% ------------------ determine mode --------------
mode  = 'absorbed';
if isfield(params,'power_model') && ~isempty(params.power_model)
    mode = lower(string(params.power_model));
end

% ------------------ IMPROVED absorbed mode --------------
if mode == "absorbed"
    PRF = max(params.PRF, 0);
    
    % ========== POWER-DEPENDENT ABSORPTION (FIX) ==========
    % Physical basis: 
    %   - Low density → poor coupling → some power lost in coil/matching network
    %   - High density → good coupling → nearly full absorption
    
    % Threshold power (below this, coupling degrades)
    P_threshold = 300;  % W (tunable based on your data)
    
    if PRF < P_threshold
        % Low power regime: poor impedance matching
        % Linear increase from 40% to 80%
        eta_abs = 0.40 + 0.40 * (PRF / P_threshold);
    else
        % High power regime: good impedance matching  
        % Asymptotic approach to 95%
        excess_power = min((PRF - P_threshold) / 1000, 1.0);
        eta_abs = 0.80 + 0.15 * excess_power;
    end
    
    % Allow user override with fixed absorption fraction
    if isfield(params,'absorb_frac') && ~isempty(params.absorb_frac)
        if params.absorb_frac >= 0 && params.absorb_frac <= 1
            eta_abs = params.absorb_frac;  % Use fixed value if provided
        end
    end
    
    % Calculate absorbed power
    Pabs = eta_abs * PRF;
    
    % Store efficiency for diagnostics/plotting
    if isfield(params, 'store_diagnostics') && params.store_diagnostics
        params.eta_absorption = eta_abs;
        params.PRF_input = PRF;
        params.Pabs_actual = Pabs;
        params.Plost_coil = PRF - Pabs;
    end
    
    return
end

% ------------------ circuit mode (unchanged) -----------------------
% geometry
R = geom.R; L = geom.L;

% gas density & electron momentum collision frequency
Tg = max(Tg_K, 50);
n_g = params.p0/(kB*Tg);

% Use the elastic momentum-transfer coefficient from rates (m^3/s)
Kel = max(Rrates.k1, 0);
nu_m = n_g * Kel;                 % 1/s

% complex plasma permittivity
wpe2  = ne*qe^2/(me*eps0);
eps_p = 1 - wpe2 ./ (omega*(omega - 1i*nu_m));

% axial wave number in plasma & Bessel boundary factor
k0  = omega/c0;
kpl = k0 * sqrt(eps_p);

% Avoid singularities for very small |kpl*R|
x = kpl*R;
J0 = besselj(0, x);
J1 = besselj(1, x);
bfr = 1i*x .* J1 ./ (eps_p .* max(J0, 1e-30));   % dimensionless load factor

% plasma "resistance" seen by the coil (empirical ICP form)
Nturns = max(getfield_or(params,'Nturns',5), 1);
Rind   = real( (2*pi*Nturns^2) / (L*omega*eps0) * bfr );  % ohms
Rind   = max(Rind, 0);                                    % clamp

% series coil resistance
Rcoil  = max(getfield_or(params,'Rcoil',0), 0);

% current from generator power PRF into series (Rind + Rcoil)
PRF = max(getfield_or(params,'PRF',0), 0);
den = max(Rind + Rcoil, 1e-9);
I2  = 2*PRF / den;               % since PRF = 0.5*(Rind+Rcoil)*I^2

% power absorbed in plasma is the portion in Rind
Pabs = 0.5 * Rind * I2;          % [W]

% optional overall coupling fudge factor (kept at 1 unless provided)
eta = getfield_or(params,'coupling_eta', 1.0);
Pabs = max(real(Pabs)*eta, 0);
end

% ---------- helper ----------
function v = getfield_or(S, name, default)
if isstruct(S) && isfield(S, name) && ~isempty(S.(name))
    v = S.(name);
else
    v = default;
end
end