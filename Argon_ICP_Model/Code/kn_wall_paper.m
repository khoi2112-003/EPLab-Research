function k_wall = kn_wall_paper(Tg, geom, p0, gamma_X, sigma_nn)
% Paper Eq. (8): k_wall = [ Λ_n^2/D_n  +  2V(2-γ)/(A v_g γ) ]^{-1}
% Λ_n (Eq. 9), D_n ≈ (1/3) v_th λ_n  (collisional regime)

kB  = 1.380649e-23;  MAr = 39.948*1.66053906660e-27;

R = geom.R;  L = geom.L;
A = 2*pi*R^2 + 2*pi*R*L;     % m^2
V = pi*R^2*L;                % m^3

% neutral number density (background) from gas law at Tg
TgK = max(Tg,50);
n_g = p0 / (kB*Tg);                                  % m^-3

% neutral thermal speed & mean free path
v_th    = sqrt(8*kB*Tg/(pi*MAr));                    % m/s
lambda_n = 1 / max(n_g*max(sigma_nn,1e-25),1e-6);  % m

% diffusion coefficient (collisional)
D_n = (1/3) * v_th * lambda_n;                        % m^2/s

% effective diffusion length Λ_n (paper Eq. 9)
Lambda_n = 1 / sqrt( (pi/L)^2 + (2.405/R)^2 );        % m

% mean molecular speed used in the surface term
v_g = v_th;                                           % m/s

% assemble k_wall
term_diff = (Lambda_n^2) / max(D_n,1e-30);            % s
term_surf = (2*V*(2 - gamma_X)) / max(A*v_g*gamma_X,1e-30); % s
k_wall    = 1 / (term_diff + term_surf);              % 1/s
end
