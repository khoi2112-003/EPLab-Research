function Rind = calc_rind(params, ne, nu_m)
% Calculate the induced resistance Rind from equation (18) of the paper
% Inputs:
%   params - parameter structure with omega, Nturns, geom
%   ne - electron density (m^-3)
%   nu_m - electron-neutral collision frequency (s^-1)
% Output:
%   Rind - induced resistance (ohms)

% Constants
eps0 = 8.8541878128e-12;  % permittivity of free space
me   = 9.10938356e-31;     % electron mass
qe   = 1.602176634e-19;    % elementary charge
c0   = 299792458;          % speed of light

% Extract parameters
omega = params.omega;
N = params.Nturns;
R = params.geom.R;  % chamber radius
L = params.geom.L;  % chamber length

% Plasma frequency squared
wpe2 = ne * qe^2 / (me * eps0);

% Complex plasma permittivity (equation 20)
eps_p = 1 - wpe2 / (omega * (omega - 1i*nu_m));

% Wave number in vacuum and plasma
k0 = omega / c0;
k = k0 * sqrt(eps_p);  % wave number in plasma

% Calculate Bessel functions for kR
x = k * R;
J0_kR = besselj(0, x);
J1_kR = besselj(1, x);

% Avoid division by zero
if abs(J0_kR) < 1e-30
    J0_kR = 1e-30;
end

% Calculate Rind from equation (18)
% Rind = (2π N²)/(Lωε₀) * Re[(ikRJ₁(kR))/(εₚJ₀(kR))]
factor = (2 * pi * N^2) / (L * omega * eps0);
bessel_ratio = (1i * k * R * J1_kR) / (eps_p * J0_kR);
Rind = factor * real(bessel_ratio);

% Ensure non-negative
Rind = max(Rind, 0);

end