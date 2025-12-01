function Pabs = icp_circuit_power(params, ne, Te_eV, Tg_K, Rrates, geom)
% Return TOTAL power absorbed by the plasma [W].
% Modes:
%   params.power_model = 'absorbed'  ->  Pabs = absorb_frac * PRF
%   params.power_model = 'circuit'   ->  simple ICP circuit using coil+plasma
%
% Notes:
% - 'absorbed' mode corresponds to the paper's use of "RF power" as the
%   power deposited in the plasma volume in the energy balance.
% - 'circuit' mode is provided only as an optional exploratory model.

% ------------------ common constants ------------------
eps0 = 8.8541878128e-12;
mu0  = 4*pi*1e-7;
c0   = 299792458;
me   = 9.10938356e-31;
qe   = 1.602176634e-19;
kB   = 1.380649e-23;

omega = params.omega;

% ------------------ paper (absorbed) mode --------------
mode  = 'absorbed';
if isfield(params,'power_model') && ~isempty(params.power_model)
    mode = lower(string(params.power_model));
end

if mode == "absorbed"
    f = 1.0;
    if isfield(params,'absorb_frac') && ~isempty(params.absorb_frac)
        f = max(params.absorb_frac,0);
    end
    Pabs = f * max(params.PRF,0);   % [W], directly used in electron energy eq.
    return
end

% ------------------ circuit mode -----------------------
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
