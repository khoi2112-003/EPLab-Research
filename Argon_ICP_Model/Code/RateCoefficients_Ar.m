function R = RateCoefficients_Ar(~, ~, Te)
% Argon reaction set (m^3/s for two-body; s^-1 for k18..k20)
Te = max(Te, 0.05);

% Elastic
R.k1  = 2.3e-14 .* Te.^1.61 .* exp( 0.06*(log(Te)).^2 - 0.12*(log(Te)).^3 );

% Direct Ionization
R.k2  = 2.39e-14 .* Te.^0.57 .* exp(-17.43./Te);  % Ar  -> Ar+ + 2e

% Multistep Ionization
R.k3  = 2.71e-13 .* Te.^0.26 .* exp( -4.59./Te);  % Arm -> Ar+ + 2e
R.k4  = 2.71e-13 .* Te.^0.28 .* exp( -4.24./Te);  % Arr -> Ar+ + 2e
R.k5  = 1.09e-12 .* Te.^0.29 .* exp( -3.42./Te);  % Arp -> Ar+ + 2e

% Excitation (to 4s m/r; to 4p)
R.k6  = 2.5e-15  .* Te.^0.74  .* exp(-11.56./Te); % Ar -> Arm
R.k7  = 3.93e-15 .* Te.^0.46  .* exp(-12.09./Te); % Ar -> Arr
R.k8  = 8.91e-15 .* Te.^(-0.04).* exp(-14.18./Te);% Ar -> Arp

% De-excitation to ground
R.k9  = 2.25e-16 .* Te.^(-0.17).* exp( -1.65./Te); % Arm -> Ar
R.k10 = 6.82e-16 .* Te.^( 0.44).* exp( -0.43./Te); % Arr -> Ar
R.k11 = 2.97e-16 .* Te.^(-0.11).* exp( -1.38./Te); % Arp -> Ar

% Excited→Excited transfers
R.k12 = 3.7e-13;                                   % Arm -> Arr
R.k13 = 2.39e-12 .* Te.^(-0.15).* exp( -1.82./Te); % Arm -> Arp
R.k14 = 2.48e-12 .* Te.^(-0.16).* exp( -1.79./Te); % Arr -> Arp
R.k15 = 9.1e-13;                                   % Arr -> Arm
R.k16 = 1.5e-13 .* Te.^( 0.51);                    % Arp -> Arm
R.k17 = 1.5e-13 .* Te.^( 0.51);                    % Arp -> Arr

% Radiative (first-order)
R.k18 = 1.0e5;  % Arr -> Ar + hv
R.k19 = 3.0e7;  % Arp -> Arm + hv
R.k20 = 3.0e7;  % Arp -> Arr + hv
end