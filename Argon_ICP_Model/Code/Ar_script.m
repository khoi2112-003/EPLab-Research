%% Plot Densities vs Absorbed Power (Pabs)

clear; clc; close all

% --- geometry & constants
geom.R   = 0.06;
geom.L   = 0.10;
geom.Tg0 = 300;
p0       = 0.266645;

% --- params
params.geom    = geom;
params.p0      = p0;
params.Vgrid   = 1000;
params.beta_i  = 0.7;
params.beta_g  = 0.3;
params.sigma_i = 5e-18;   % m^2 (optimized value)
params.Tg0     = geom.Tg0;
params.gamma_m = 0.2;
params.gamma_r = 0.5;
params.gamma_p = 1;
params.gamma_g = 1e-5;
params.Q0 = 1.13e18;
params.s_sheath = 1.5e-3;
params.sigma_nn = 4e-19;
params.ps = struct();
params.f_RF   = 13.56e6;
params.omega  = 2*pi*params.f_RF;
params.Nturns = 5;
params.Rcoil  = 2;
params.Rc     = 0.07;
params.lc     = geom.L;
params.f_sheath = 0.55;

% --- Constants
kB = 1.380649e-23; 
qe = 1.602176634e-19; 
MAr = 39.948*1.66053906660e-27;
V = pi*geom.R^2*geom.L;

% --- Initial conditions at low power
iAr=1; iArm=2; iArr=3; iArp=4; iI=5; iTe=6; iTg=7;

nAr_0  = 6.4e19;
nArm_0 = 1e17;
nArr_0 = 7e16;
nArp_0 = 1e15;
ne_0   = 3e17;
Te_0   = 3.6;
Tg_0   = 310;

y0 = [nAr_0; nArm_0; nArr_0; nArp_0; ne_0; Te_0; Tg_0];

% Power sweep - adjust to cover Pabs = 0-600W
% Since Pabs ≈ 0.3-0.5 * PRF, use PRF from 0 to ~1800W
PRF_vec = linspace(10, 100, 30).';  % Start from 50W to avoid convergence issues
N = length(PRF_vec);

% Storage
Y = nan(N, 7);
Pabs_vec = nan(N, 1);
PRF_actual = nan(N, 1);

% ODE options
opts = odeset('RelTol', 1e-7, ...
              'AbsTol', [1e14 1e13 1e13 1e12 1e13 1e-5 0.01], ...
              'NonNegative', 1:7);

fprintf('Running power sweep for Pabs vs densities...\n');
fprintf('PRF (W)  Pabs (W)   ne (m^-3)    Te (eV)    Tg (K)\n');
fprintf('--------------------------------------------------------\n');

% First point
params.PRF = PRF_vec(1);
tspan = [0, 0.1];
sol = ode15s(@(t,z) rhs_global(t,z,params,[]), tspan, y0, opts);
y = sol.y(:,end);
Y(1,:) = y.';
PRF_actual(1) = params.PRF;

% Calculate Pabs for first point
ne_temp = y(iI);
Tg_temp = y(iTg);
n_g_temp = y(iAr) + y(iArm) + y(iArr) + y(iArp);
R_temp = RateCoefficients_Ar(params.Tg0, params.p0, y(iTe));
Kel_temp = max(R_temp.k1, 0);
nu_m_temp = n_g_temp * Kel_temp;
Rind_temp = calc_rind(params, ne_temp, nu_m_temp);
I2_temp = 2*params.PRF / (Rind_temp + params.Rcoil);
Pcoil_temp = 0.5 * params.Rcoil * I2_temp;
Pabs_vec(1) = (params.PRF - Pcoil_temp);

fprintf('%6.0f   %6.1f   %.2e   %.3f   %.1f\n', ...
        PRF_actual(1), Pabs_vec(1), y(iI), y(iTe), y(iTg));

% Continue sweep
for i = 2:N
    params.PRF = PRF_vec(i);
    tspan = [0, 0.02];
    
    sol = ode15s(@(t,z) rhs_global(t,z,params,[]), tspan, y, opts);
    y = sol.y(:,end);
    Y(i,:) = y.';
    PRF_actual(i) = params.PRF;
    
    % Calculate Pabs
    ne_temp = y(iI);
    Tg_temp = y(iTg);
    n_g_temp = y(iAr) + y(iArm) + y(iArr) + y(iArp);
    R_temp = RateCoefficients_Ar(params.Tg0, params.p0, y(iTe));
    Kel_temp = max(R_temp.k1, 0);
    nu_m_temp = n_g_temp * Kel_temp;
    Rind_temp = calc_rind(params, ne_temp, nu_m_temp);
    I2_temp = 2*params.PRF / (Rind_temp + params.Rcoil);
    Pcoil_temp = 0.5 * params.Rcoil * I2_temp;
    Pabs_vec(i) = (params.PRF - Pcoil_temp);
    
    fprintf('%6.0f   %6.1f   %.2e   %.3f   %.1f\n', ...
            PRF_actual(i), Pabs_vec(i), y(iI), y(iTe), y(iTg));
end

%% Print densities at 15 W (add after the power sweep loop, before plotting)

%% Extract exact densities at Pabs = 15 W (add after the power sweep loop)

% Interpolate all variables at exactly Pabs = 15W
Pabs_target = 15.0;

nAr_15W  = interp1(Pabs_vec, Y(:,iAr), Pabs_target, 'linear');
nArm_15W = interp1(Pabs_vec, Y(:,iArm), Pabs_target, 'linear');
nArr_15W = interp1(Pabs_vec, Y(:,iArr), Pabs_target, 'linear');
nArp_15W = interp1(Pabs_vec, Y(:,iArp), Pabs_target, 'linear');
ne_15W   = interp1(Pabs_vec, Y(:,iI), Pabs_target, 'linear');
Te_15W   = interp1(Pabs_vec, Y(:,iTe), Pabs_target, 'linear');
Tg_15W   = interp1(Pabs_vec, Y(:,iTg), Pabs_target, 'linear');

fprintf('\n=== DENSITIES AT Pabs = %.2f W (interpolated) ===\n', Pabs_target);
fprintf('  nAr  = %.4e m^-3\n', nAr_15W);
fprintf('  nArm = %.4e m^-3\n', nArm_15W);
fprintf('  nArr = %.4e m^-3\n', nArr_15W);
fprintf('  nArp = %.4e m^-3\n', nArp_15W);
fprintf('  ne   = %.4e m^-3\n', ne_15W);
fprintf('  Te   = %.3f eV\n', Te_15W);
fprintf('  Tg   = %.1f K\n\n', Tg_15W);



%% Create Plots

% Extract densities
nAr = Y(:,iAr);
nArm = Y(:,iArm);
nArr = Y(:,iArr);
nArp = Y(:,iArp);
ne = Y(:,iI);
Te = Y(:,iTe);
Tg = Y(:,iTg);

% Filter to Pabs = 0-600W range
idx = Pabs_vec <= 600;
Pabs_plot = Pabs_vec(idx);
nAr_plot = nAr(idx);
nArm_plot = nArm(idx);
nArr_plot = nArr(idx);
nArp_plot = nArp(idx);
ne_plot = ne(idx);
Te_plot = Te(idx);
Tg_plot = Tg(idx);

figure('Position', [100 100 1200 800]);

% Plot 1: All neutral species
subplot(2,2,1)
hold on
plot(Pabs_plot, nAr_plot, 'b-o', 'LineWidth', 2, 'MarkerSize', 6)
plot(Pabs_plot, nArm_plot, 'r-s', 'LineWidth', 2, 'MarkerSize', 6)
plot(Pabs_plot, nArr_plot, 'g-d', 'LineWidth', 2, 'MarkerSize', 6)
plot(Pabs_plot, nArp_plot, 'm-^', 'LineWidth', 2, 'MarkerSize', 6)
hold off
xlabel('Absorbed Power, P_{abs} (W)', 'FontSize', 12)
xlim = [0 100];
xline(15, 'k--', 'LineWidth', 1.5, 'Label', '15W');
ylabel('Density (m^{-3})', 'FontSize', 12)
title('Neutral Species Densities vs Absorbed Power', 'FontSize', 14)
legend('Ar (ground)', 'Ar^m (metastable)', 'Ar^r (resonant)', 'Ar^p (4p)', ...
       'Location', 'best', 'FontSize', 10)
grid on
set(gca, 'FontSize', 11)

% Plot 2: Electron density
subplot(2,2,2)
plot(Pabs_plot, ne_plot, 'b-o', 'LineWidth', 2, 'MarkerSize', 6)
xlabel('Absorbed Power, P_{abs} (W)', 'FontSize', 12)
xline(15, 'k--', 'LineWidth', 1.5, 'Label', '15W');
ylabel('Electron Density (m^{-3})', 'FontSize', 12)
title('Electron Density vs Absorbed Power', 'FontSize', 14)
grid on
set(gca, 'FontSize', 11)

% Plot 3: Excited species (log scale)
subplot(2,2,3)
semilogy(Pabs_plot, nArm_plot, 'r-s', 'LineWidth', 2, 'MarkerSize', 6)
hold on
semilogy(Pabs_plot, nArr_plot, 'g-d', 'LineWidth', 2, 'MarkerSize', 6)
semilogy(Pabs_plot, nArp_plot, 'm-^', 'LineWidth', 2, 'MarkerSize', 6)
hold off
xlabel('Absorbed Power, P_{abs} (W)', 'FontSize', 12)
xline(15, 'k--', 'LineWidth', 1.5, 'Label', '15W');
ylabel('Density (m^{-3})', 'FontSize', 12)
title('Excited State Densities vs Absorbed Power (log scale)', 'FontSize', 14)
legend('Ar^m', 'Ar^r', 'Ar^p', 'Location', 'best', 'FontSize', 10)
grid on
set(gca, 'FontSize', 11)

% Plot 4: Temperatures
subplot(2,2,4)
yyaxis left
plot(Pabs_plot, Te_plot, 'b-o', 'LineWidth', 2, 'MarkerSize', 6)
ylabel('Electron Temperature, T_e (eV)', 'FontSize', 12, 'Color', 'b')
set(gca, 'YColor', 'b')

yyaxis right
plot(Pabs_plot, Tg_plot, 'r-s', 'LineWidth', 2, 'MarkerSize', 6)
ylabel('Gas Temperature, T_g (K)', 'FontSize', 12, 'Color', 'r')
set(gca, 'YColor', 'r')

xlabel('Absorbed Power, P_{abs} (W)', 'FontSize', 12)
title('Temperatures vs Absorbed Power', 'FontSize', 14)
legend('T_e', 'T_g', 'Location', 'best', 'FontSize', 10)
grid on
set(gca, 'FontSize', 11)

sgtitle('Species Densities and Temperatures vs Absorbed Power', 'FontSize', 16, 'FontWeight', 'bold')

%% Additional plot: All densities on one log plot
figure('Position', [150 150 800 600]);
semilogy(Pabs_plot, nAr_plot, 'b-o', 'LineWidth', 2, 'MarkerSize', 6)
hold on
semilogy(Pabs_plot, nArm_plot, 'r-s', 'LineWidth', 2, 'MarkerSize', 6)
semilogy(Pabs_plot, nArr_plot, 'g-d', 'LineWidth', 2, 'MarkerSize', 6)
semilogy(Pabs_plot, nArp_plot, 'm-^', 'LineWidth', 2, 'MarkerSize', 6)
semilogy(Pabs_plot, ne_plot, 'k-*', 'LineWidth', 2, 'MarkerSize', 8)
hold off
xline(15, 'k--', 'LineWidth', 2, 'Label', 'P_{abs} = 15W')
xlabel('Absorbed Power, P_{abs} (W)', 'FontSize', 14)
ylabel('Density (m^{-3})', 'FontSize', 14)
title('All Species Densities vs Absorbed Power', 'FontSize', 16, 'FontWeight', 'bold')
legend('Ar', 'Ar^m', 'Ar^r', 'Ar^p', 'n_e', 'Location', 'best', 'FontSize', 12)
grid on
set(gca, 'FontSize', 12)

%% Print summary statistics
fprintf('\n=== SUMMARY ===\n');
fprintf('Pabs range: %.1f - %.1f W\n', min(Pabs_plot), max(Pabs_plot));
fprintf('At Pabs = %.1f W:\n', Pabs_plot(1));
fprintf('  ne = %.2e m^-3, Te = %.2f eV, Tg = %.1f K\n', ne_plot(1), Te_plot(1), Tg_plot(1));
fprintf('At Pabs = %.1f W:\n', Pabs_plot(end));
fprintf('  ne = %.2e m^-3, Te = %.2f eV, Tg = %.1f K\n', ne_plot(end), Te_plot(end), Tg_plot(end));
fprintf('Electron density increase: %.1fx\n', ne_plot(end)/ne_plot(1));
fprintf('Te change: %.3f eV\n', Te_plot(end) - Te_plot(1));
fprintf('Tg change: %.1f K\n', Tg_plot(end) - Tg_plot(1));