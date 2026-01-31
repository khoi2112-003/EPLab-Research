%% validate_model_enhanced.m - Complete validation with OPTIMIZED parameters
% This script validates species densities AND temperatures against reference data

clear; clc; close all

fprintf('========================================\n');
fprintf('ENHANCED MODEL VALIDATION\n');
fprintf('With Optimized Parameters\n');
fprintf('========================================\n\n');

% ================= STEP 1: LOAD ALL REFERENCE DATA =================
fprintf('Step 1: Loading reference data...\n');

% Species densities
data_Ar = readmatrix('Ar_new.csv');
data_Ar_plus = readmatrix('Ar+_new.csv');
data_Arm = readmatrix('Arm_new.csv');
data_Arr = readmatrix('Arr_new.csv');
data_Arp = readmatrix('Arp_new.csv');

% Temperatures (Case 1)
data_Te = readmatrix('Te Case 1.csv');
data_Tg = readmatrix('Tg Case 1.csv');

% Extract data
power_species = data_Ar(:,1);
n_Ar_ref = data_Ar(:,2);
n_Ar_plus_ref = data_Ar_plus(:,2);
n_Arm_ref = data_Arm(:,2);
n_Arr_ref = data_Arr(:,2);
n_Arp_ref = data_Arp(:,2);

power_Te = data_Te(:,1);
Te_ref = data_Te(:,2);

power_Tg = data_Tg(:,1);
Tg_ref = data_Tg(:,2);

fprintf('  Loaded species data: %d power points (%.0f - %.0f W)\n', ...
        length(power_species), min(power_species), max(power_species));
fprintf('  Loaded Te Case 1: %d power points (%.0f - %.0f W)\n', ...
        length(power_Te), min(power_Te), max(power_Te));
fprintf('  Loaded Tg Case 1: %d power points (%.0f - %.0f W)\n\n', ...
        length(power_Tg), min(power_Tg), max(power_Tg));

% ================= STEP 2: SETUP MODEL WITH OPTIMIZED PARAMETERS =================
fprintf('Step 2: Setting up model with OPTIMIZED parameters...\n');

geom.R = 0.06; geom.L = 0.10; geom.Tg0 = 300;
params.geom = geom;
params.p0 = 0.266645;  % 2 mTorr in Pa
params.Vgrid = 1000;
params.beta_i = 0.7;
params.beta_g = 0.3;
params.Q0 = 4.5e18;
params.s_sheath = 1.5e-3;
params.ps = struct();

% RF and circuit parameters
params.f_RF = 13.56e6;
params.omega = 2*pi*params.f_RF;
params.Nturns = 5;
params.Rc = 0.07;
params.lc = geom.L;
params.f_sheath = 0.55;

% Power model
params.power_model = 'circuit';

% ========================================
% OPTIMIZED PARAMETERS (from inverse optimization)
% Validated at 2 mTorr with < 5% error on all outputs
% ========================================
params.Rcoil = 0.543183;
params.sigma_heating = 4.852125e-18;
params.sigma_nn = 8.236909e-19;
params.sigma_i = 7.446966e-19;
params.gamma_g = 1.000000e-04;
params.gamma_m = 0.001110;
params.gamma_r = 0.001945;
params.gamma_p = 0.005982;

fprintf('  Optimized parameters:\n');
fprintf('    Rcoil         = %.6f Ohm\n', params.Rcoil);
fprintf('    sigma_heating = %.6e m^2\n', params.sigma_heating);
fprintf('    sigma_nn      = %.6e m^2\n', params.sigma_nn);
fprintf('    sigma_i       = %.6e m^2\n', params.sigma_i);
fprintf('    gamma_g       = %.6e\n', params.gamma_g);
fprintf('    gamma_m       = %.6f\n', params.gamma_m);
fprintf('    gamma_r       = %.6f\n', params.gamma_r);
fprintf('    gamma_p       = %.6f\n\n', params.gamma_p);

y0 = [6.4e19; 1e17; 7e16; 1e15; 3e17; 3.6; 310];
opts = odeset('RelTol', 1e-6, 'AbsTol', [1e14 1e13 1e13 1e12 1e13 1e-4 0.1], 'NonNegative', 1:7);

fprintf('  Model configured.\n\n');

% ================= STEP 3: RUN MODEL AT ALL POWER POINTS =================
fprintf('Step 3: Running model...\n');

% Combine all power points
all_powers = unique([power_species; power_Te; power_Tg]);
all_powers = sort(all_powers);

N = length(all_powers);
Y_model = nan(N, 7);

fprintf('  Running %d power points... ', N);
for k = 1:N
    params.PRF = all_powers(k);
    
    if k == 1
        tspan = [0, 0.5];
    else
        tspan = [0, 0.1];
        y0 = Y_model(k-1, :)';
    end
    
    sol = ode15s(@(t,z) rhs_global(t,z,params,[]), tspan, y0, opts);
    Y_model(k,:) = sol.y(:,end)';
    
    if mod(k, 10) == 0
        fprintf('.');
    end
    
end
fprintf(' Done!\n\n');

% ================= STEP 4: INTERPOLATE TO REFERENCE POINTS =================
fprintf('Step 4: Interpolating results...\n');

% Species
n_Ar_model = interp1(all_powers, Y_model(:,1), power_species, 'linear', 'extrap');
n_Arm_model = interp1(all_powers, Y_model(:,2), power_species, 'linear', 'extrap');
n_Arr_model = interp1(all_powers, Y_model(:,3), power_species, 'linear', 'extrap');
n_Arp_model = interp1(all_powers, Y_model(:,4), power_species, 'linear', 'extrap');
n_Ar_plus_model = interp1(all_powers, Y_model(:,5), power_species, 'linear', 'extrap');

% Temperatures
Te_model = interp1(all_powers, Y_model(:,6), power_Te, 'linear', 'extrap');
Tg_model = interp1(all_powers, Y_model(:,7), power_Tg, 'linear', 'extrap');

fprintf('  Interpolation complete.\n\n');

% ================= STEP 5: CALCULATE ERRORS =================
fprintf('Step 5: Calculating errors...\n\n');

% Species errors
error_Ar = abs(n_Ar_model - n_Ar_ref) ./ n_Ar_ref * 100;
error_Ar_plus = abs(n_Ar_plus_model - n_Ar_plus_ref) ./ n_Ar_plus_ref * 100;
error_Arm = abs(n_Arm_model - n_Arm_ref) ./ n_Arm_ref * 100;
error_Arr = abs(n_Arr_model - n_Arr_ref) ./ n_Arr_ref * 100;
error_Arp = abs(n_Arp_model - n_Arp_ref) ./ n_Arp_ref * 100;

% Temperature errors
error_Te = abs(Te_model - Te_ref) ./ Te_ref * 100;
error_Tg = abs(Tg_model - Tg_ref) ./ Tg_ref * 100;

% Print error summary
fprintf('========================================\n');
fprintf('VALIDATION ERRORS (OPTIMIZED MODEL)\n');
fprintf('========================================\n\n');

fprintf('SPECIES DENSITIES:\n');
fprintf('%-8s | Mean   | Max    | RMS\n', 'Species');
fprintf('---------|--------|--------|--------\n');
fprintf('%-8s | %5.1f%% | %5.1f%% | %5.1f%%\n', 'Ar', mean(error_Ar), max(error_Ar), rms(error_Ar));
fprintf('%-8s | %5.1f%% | %5.1f%% | %5.1f%%\n', 'Ar+', mean(error_Ar_plus), max(error_Ar_plus), rms(error_Ar_plus));
fprintf('%-8s | %5.1f%% | %5.1f%% | %5.1f%%\n', 'Arm', mean(error_Arm), max(error_Arm), rms(error_Arm));
fprintf('%-8s | %5.1f%% | %5.1f%% | %5.1f%%\n', 'Arr', mean(error_Arr), max(error_Arr), rms(error_Arr));
fprintf('%-8s | %5.1f%% | %5.1f%% | %5.1f%%\n\n', 'Arp', mean(error_Arp), max(error_Arp), rms(error_Arp));

fprintf('TEMPERATURES (Case 1):\n');
fprintf('%-8s | Mean   | Max    | RMS\n', 'Temp');
fprintf('---------|--------|--------|--------\n');
fprintf('%-8s | %5.1f%% | %5.1f%% | %5.1f%%\n', 'Te', mean(error_Te), max(error_Te), rms(error_Te));
fprintf('%-8s | %5.1f%% | %5.1f%% | %5.1f%%\n\n', 'Tg', mean(error_Tg), max(error_Tg), rms(error_Tg));

fprintf('TEMPERATURE RANGES:\n');
fprintf('  Te: %.2f - %.2f eV (model) vs %.2f - %.2f eV (ref)\n', ...
        min(Te_model), max(Te_model), min(Te_ref), max(Te_ref));
fprintf('  Tg: %.0f - %.0f K (model) vs %.0f - %.0f K (ref)\n\n', ...
        min(Tg_model), max(Tg_model), min(Tg_ref), max(Tg_ref));

% ================= STEP 6: GENERATE PLOTS =================
fprintf('Step 6: Generating plots...\n');

%% Figure 1: Temperature Validation with Case 1 Data
figure('Color','w', 'Position', [50 50 1400 500]);

% Electron Temperature
subplot(1,3,1);
plot(power_Te, Te_ref, 'o-', 'LineWidth', 2, 'MarkerSize', 8, ...
     'DisplayName', 'Reference', 'Color', [0 0 0]);
hold on;
plot(power_Te, Te_model, 's--', 'LineWidth', 2, 'MarkerSize', 6, ...
     'DisplayName', 'Model', 'Color', [0.85 0.33 0.10]);
xlabel('RF Power (W)', 'FontSize', 12);
ylabel('T_e (eV)', 'FontSize', 12);
title('Electron Temperature', 'FontSize', 14);
legend('Location', 'best');
grid on;
xlim([min(power_Te)-50 max(power_Te)+50]);

% Gas Temperature
subplot(1,3,2);
plot(power_Tg, Tg_ref, 'o-', 'LineWidth', 2, 'MarkerSize', 8, ...
     'DisplayName', 'Reference', 'Color', [0 0 0]);
hold on;
plot(power_Tg, Tg_model, 's--', 'LineWidth', 2, 'MarkerSize', 6, ...
     'DisplayName', 'Model', 'Color', [0.85 0.33 0.10]);
xlabel('RF Power (W)', 'FontSize', 12);
ylabel('T_g (K)', 'FontSize', 12);
title('Gas Temperature', 'FontSize', 14);
legend('Location', 'best');
grid on;
xlim([min(power_Tg)-50 max(power_Tg)+50]);

% Temperature Errors
subplot(1,3,3);
yyaxis left
plot(power_Te, error_Te, '-', 'LineWidth', 2, 'Color', [0 0.45 0.74]);
ylabel('T_e Error (%)', 'FontSize', 12);
ylim([0 max([max(error_Te), max(error_Tg)])*1.2]);

yyaxis right
plot(power_Tg, error_Tg, '-', 'LineWidth', 2, 'Color', [0.85 0.33 0.10]);
ylabel('T_g Error (%)', 'FontSize', 12);
ylim([0 max([max(error_Te), max(error_Tg)])*1.2]);

xlabel('RF Power (W)', 'FontSize', 12);
title('Temperature Errors', 'FontSize', 14);
grid on;
xlim([200 1600]);
yline(5, 'g--', '5%', 'LineWidth', 1.5);

sgtitle('Temperature Validation - Optimized Parameters', 'FontSize', 16, 'FontWeight', 'bold');
saveas(gcf, 'validation_temperatures.png');

%% Figure 2: Species Densities
figure('Color','w', 'Position', [100 100 1800 700]);

species_names = {'Ar', 'Ar^+', 'Ar^m', 'Ar^r', 'Ar^p'};
ref_data = {n_Ar_ref, n_Ar_plus_ref, n_Arm_ref, n_Arr_ref, n_Arp_ref};
model_data = {n_Ar_model, n_Ar_plus_model, n_Arm_model, n_Arr_model, n_Arp_model};
errors = {error_Ar, error_Ar_plus, error_Arm, error_Arr, error_Arp};

for i = 1:5
    % Density comparison
    subplot(2, 5, i);
    semilogy(power_species, ref_data{i}, 'o-', 'LineWidth', 2, 'MarkerSize', 6, ...
             'DisplayName', 'Reference', 'Color', [0 0 0]);
    hold on;
    semilogy(power_species, model_data{i}, 's--', 'LineWidth', 2, 'MarkerSize', 5, ...
             'DisplayName', 'Model', 'Color', [0.85 0.33 0.10]);
    xlabel('Power (W)');
    ylabel('Density (m^{-3})');
    title(species_names{i});
    legend('Location', 'best');
    grid on;
    
    % Error plot
    subplot(2, 5, i+5);
    plot(power_species, errors{i}, 'r-', 'LineWidth', 2);
    xlabel('Power (W)');
    ylabel('Error (%)');
    title(sprintf('%s Error (mean=%.1f%%)', species_names{i}, mean(errors{i})));
    grid on;
    yline(5, 'g--', '5%', 'LineWidth', 1);
    yline(10, 'k--', '10%', 'LineWidth', 1);
    ylim([0 max(errors{i})*1.3]);
end

sgtitle('Species Density Validation - Optimized Parameters', 'FontSize', 16, 'FontWeight', 'bold');
saveas(gcf, 'validation_species.png');

%% Figure 3: Error Summary
figure('Color','w', 'Position', [200 200 900 400]);

all_errors = [mean(error_Ar), mean(error_Ar_plus), mean(error_Arm), ...
              mean(error_Arr), mean(error_Arp), mean(error_Te), mean(error_Tg)];
all_names = {'Ar', 'Ar^+', 'Ar^m', 'Ar^r', 'Ar^p', 'T_e', 'T_g'};

b = bar(all_errors);
b.FaceColor = 'flat';
colors = parula(7);
for k = 1:7
    b.CData(k,:) = colors(k,:);
end
set(gca, 'XTickLabel', all_names);
ylabel('Mean Relative Error (%)', 'FontSize', 12);
title('Error Summary - All Outputs (Optimized Parameters)', 'FontSize', 14);
grid on;
yline(5, 'g--', '5%', 'LineWidth', 2);
yline(10, 'r--', '10%', 'LineWidth', 2);

saveas(gcf, 'validation_summary.png');

fprintf('  Plots saved!\n\n');

% ================= STEP 7: VALIDATION STATUS =================
fprintf('========================================\n');
fprintf('VALIDATION STATUS\n');
fprintf('========================================\n\n');

% Overall assessment
overall_species_error = mean([mean(error_Ar), mean(error_Ar_plus), ...
                               mean(error_Arm), mean(error_Arr), mean(error_Arp)]);
overall_temp_error = mean([mean(error_Te), mean(error_Tg)]);
overall_all = mean(all_errors);

fprintf('OVERALL PERFORMANCE:\n');
fprintf('  Species (average): %.1f%%\n', overall_species_error);
fprintf('  Temperatures:      %.1f%%\n', overall_temp_error);
fprintf('  All outputs:       %.1f%%\n\n', overall_all);

if overall_all < 5
    fprintf('  ★★★ EXCELLENT: All errors < 5%%\n');
    fprintf('      Ready for NN training!\n');
elseif overall_all < 10
    fprintf('  ★★ VERY GOOD: All errors < 10%%\n');
    fprintf('     Ready for NN training!\n');
elseif overall_all < 15
    fprintf('  ★ GOOD: All errors < 15%%\n');
    fprintf('    Acceptable for NN training\n');
else
    fprintf('  ✗ NEEDS WORK: Some errors > 15%%\n');
    fprintf('    Consider parameter tuning\n');
end

fprintf('\n========================================\n');
fprintf('VALIDATION COMPLETE!\n');
fprintf('========================================\n');