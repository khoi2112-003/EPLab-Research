%% inverse_extended_bounds.m - Extended Rcoil bounds
% Previous run hit Rcoil = 1.0 (lower bound), so extend it

clear; clc; close all;
warning('off', 'all');

fprintf('========================================\n');
fprintf('EXTENDED BOUNDS OPTIMIZATION\n');
fprintf('========================================\n\n');

%% Load reference data
data_Ar = readmatrix('Ar_new.csv');
data_Ar_plus = readmatrix('Ar+_new.csv');
data_Arm = readmatrix('Arm_new.csv');
data_Arr = readmatrix('Arr_new.csv');
data_Arp = readmatrix('Arp_new.csv');
data_Te = readmatrix('Te Case 1.csv');
data_Tg = readmatrix('Tg Case 1.csv');

ref.power_species = data_Ar(:,1);
ref.n_Ar = data_Ar(:,2);
ref.n_Ar_plus = data_Ar_plus(:,2);
ref.n_Arm = data_Arm(:,2);
ref.n_Arr = data_Arr(:,2);
ref.n_Arp = data_Arp(:,2);
ref.power_Te = data_Te(:,1);
ref.Te = data_Te(:,2);
ref.power_Tg = data_Tg(:,1);
ref.Tg = data_Tg(:,2);

%% Setup base parameters
geom.R = 0.06; geom.L = 0.10; geom.Tg0 = 300;

base_params.geom = geom;
base_params.p0 = 0.266645;
base_params.Vgrid = 1000;
base_params.beta_i = 0.7;
base_params.beta_g = 0.3;
base_params.Q0 = 4.5e18;
base_params.s_sheath = 1.5e-3;
base_params.ps = struct();
base_params.f_RF = 13.56e6;
base_params.omega = 2*pi*base_params.f_RF;
base_params.Nturns = 5;
base_params.Rc = 0.07;
base_params.lc = geom.L;
base_params.f_sheath = 0.55;
base_params.power_model = 'circuit';

%% Parameter definitions - EXTENDED Rcoil bounds
param_names = {'Rcoil', 'sigma_h', 'sigma_nn', 'sigma_i', 'gamma_g', 'gamma_m', 'gamma_r', 'gamma_p'};

% Start from previous best
x0 = [1.0, ...                 % Rcoil - was at bound
      log10(5.3e-18), ...      % sigma_heating
      log10(5e-19), ...        % sigma_nn
      log10(1.2e-18), ...      % sigma_i
      log10(7.6e-4), ...       % gamma_g
      log10(1.05e-3), ...      % gamma_m
      log10(1.03e-3), ...      % gamma_r
      log10(5.02e-3)];         % gamma_p

% EXTENDED bounds - Rcoil can go down to 0.3 Ohm
%           Rcoil  sig_h   sig_nn  sig_i   gam_g   gam_m   gam_r   gam_p
lb =       [0.3,   -18,    -20,    -18.5,  -4,     -4,     -4,     -3];
ub =       [3,     -16.5,  -18,    -17,    -2.5,   -2,     -2,     -1.5];

%% Weights
weights.Ar = 1.0;
weights.Ar_plus = 1.5;
weights.Arm = 1.0;
weights.Arr = 1.0;
weights.Arp = 1.5;
weights.Te = 5.0;
weights.Tg = 5.0;

%% Phase 1: Test lower Rcoil values
fprintf('Phase 1: Testing lower Rcoil values...\n');

Rcoil_vals = [0.3, 0.5, 0.7, 0.9, 1.0, 1.2];
sigma_h_vals = [3e-18, 5e-18, 8e-18, 1e-17];
sigma_i_vals = [8e-19, 1e-18, 1.5e-18, 2e-18];

best_error = inf;
best_x = x0;

for Rc = Rcoil_vals
    for sh = sigma_h_vals
        for si = sigma_i_vals
            x_test = x0;
            x_test(1) = Rc;
            x_test(2) = log10(sh);
            x_test(4) = log10(si);
            
            err = run_objective(x_test, base_params, ref, weights);
            
            if err < best_error
                best_error = err;
                best_x = x_test;
                fprintf('  Rcoil=%.2f, sigma_h=%.1e, sigma_i=%.1e -> Error=%.1f\n', ...
                        Rc, sh, si, err);
            end
        end
    end
end

fprintf('\nGrid search best: Rcoil=%.2f, Error=%.1f\n', best_x(1), best_error);

%% Phase 2: Fine-tune
fprintf('\nPhase 2: Fine-tuning...\n');

options = optimset('Display', 'iter', 'TolX', 1e-3, 'TolFun', 0.3, ...
                   'MaxIter', 250, 'MaxFunEvals', 1000);

obj_fn = @(x) run_objective_bounded(x, base_params, ref, weights, lb, ub);
[x_opt, error_opt] = fminsearch(obj_fn, best_x, options);
x_opt = max(min(x_opt, ub), lb);

%% Results
Rcoil_opt = x_opt(1);
sigma_h_opt = 10^x_opt(2);
sigma_nn_opt = 10^x_opt(3);
sigma_i_opt = 10^x_opt(4);
gamma_g_opt = 10^x_opt(5);
gamma_m_opt = 10^x_opt(6);
gamma_r_opt = 10^x_opt(7);
gamma_p_opt = 10^x_opt(8);

fprintf('\n========================================\n');
fprintf('EXTENDED BOUNDS RESULT\n');
fprintf('========================================\n');
fprintf('Optimal parameters:\n');
fprintf('  Rcoil         = %.4f Ohm\n', Rcoil_opt);
fprintf('  sigma_heating = %.4e m^2\n', sigma_h_opt);
fprintf('  sigma_nn      = %.4e m^2\n', sigma_nn_opt);
fprintf('  sigma_i       = %.4e m^2\n', sigma_i_opt);
fprintf('  gamma_g       = %.6e\n', gamma_g_opt);
fprintf('  gamma_m       = %.6f\n', gamma_m_opt);
fprintf('  gamma_r       = %.6f\n', gamma_r_opt);
fprintf('  gamma_p       = %.6f\n', gamma_p_opt);

%% Validation
warning('on', 'all');

params = base_params;
params.Rcoil = Rcoil_opt;
params.sigma_heating = sigma_h_opt;
params.sigma_nn = sigma_nn_opt;
params.sigma_i = sigma_i_opt;
params.gamma_g = gamma_g_opt;
params.gamma_m = gamma_m_opt;
params.gamma_r = gamma_r_opt;
params.gamma_p = gamma_p_opt;

[Te_model, Tg_model, n_Ar_model, n_Ar_plus_model, n_Arm_model, n_Arr_model, n_Arp_model] = ...
    run_full_model(params, ref);

error_Te = abs(Te_model - ref.Te) ./ ref.Te * 100;
error_Tg = abs(Tg_model - ref.Tg) ./ ref.Tg * 100;
error_Ar = abs(n_Ar_model - ref.n_Ar) ./ ref.n_Ar * 100;
error_Ar_plus = abs(n_Ar_plus_model - ref.n_Ar_plus) ./ ref.n_Ar_plus * 100;
error_Arm = abs(n_Arm_model - ref.n_Arm) ./ ref.n_Arm * 100;
error_Arr = abs(n_Arr_model - ref.n_Arr) ./ ref.n_Arr * 100;
error_Arp = abs(n_Arp_model - ref.n_Arp) ./ ref.n_Arp * 100;

fprintf('\n========================================\n');
fprintf('FINAL VALIDATION ERRORS\n');
fprintf('========================================\n');
fprintf('%-8s | Mean   | Max    | RMS\n', 'Output');
fprintf('---------|--------|--------|--------\n');
fprintf('%-8s | %5.1f%% | %5.1f%% | %5.1f%%\n', 'Ar', mean(error_Ar), max(error_Ar), rms(error_Ar));
fprintf('%-8s | %5.1f%% | %5.1f%% | %5.1f%%\n', 'Ar+', mean(error_Ar_plus), max(error_Ar_plus), rms(error_Ar_plus));
fprintf('%-8s | %5.1f%% | %5.1f%% | %5.1f%%\n', 'Arm', mean(error_Arm), max(error_Arm), rms(error_Arm));
fprintf('%-8s | %5.1f%% | %5.1f%% | %5.1f%%\n', 'Arr', mean(error_Arr), max(error_Arr), rms(error_Arr));
fprintf('%-8s | %5.1f%% | %5.1f%% | %5.1f%%\n', 'Arp', mean(error_Arp), max(error_Arp), rms(error_Arp));
fprintf('%-8s | %5.1f%% | %5.1f%% | %5.1f%%\n', 'Te', mean(error_Te), max(error_Te), rms(error_Te));
fprintf('%-8s | %5.1f%% | %5.1f%% | %5.1f%%\n', 'Tg', mean(error_Tg), max(error_Tg), rms(error_Tg));

%% Plots
figure('Color', 'w', 'Position', [50 50 1400 500]);

subplot(1,3,1);
plot(ref.power_Te, ref.Te, 'ko-', 'LineWidth', 2, 'MarkerSize', 6);
hold on;
plot(ref.power_Te, Te_model, 's--', 'LineWidth', 2, 'MarkerSize', 5, 'Color', [0.85 0.33 0.10]);
xlabel('RF Power (W)'); ylabel('T_e (eV)');
title('Electron Temperature'); legend('Reference', 'Model', 'Location', 'best'); grid on;

subplot(1,3,2);
plot(ref.power_Tg, ref.Tg, 'ko-', 'LineWidth', 2, 'MarkerSize', 6);
hold on;
plot(ref.power_Tg, Tg_model, 's--', 'LineWidth', 2, 'MarkerSize', 5, 'Color', [0.85 0.33 0.10]);
xlabel('RF Power (W)'); ylabel('T_g (K)');
title('Gas Temperature'); legend('Reference', 'Model', 'Location', 'best'); grid on;

subplot(1,3,3);
plot(ref.power_Te, error_Te, 'b-', 'LineWidth', 2); hold on;
plot(ref.power_Tg, error_Tg, 'r-', 'LineWidth', 2);
xlabel('RF Power (W)'); ylabel('Error (%)');
title('Temperature Errors vs Power');
legend('T_e', 'T_g'); 
yline(5, 'k--', '5%'); grid on;

sgtitle(sprintf('Extended Bounds: Rcoil=%.2f, \\sigma_h=%.1e, \\sigma_i=%.1e', ...
        Rcoil_opt, sigma_h_opt, sigma_i_opt));
saveas(gcf, 'inverse_extended_temps.png');

% Species plot with densities AND errors
figure('Color', 'w', 'Position', [100 100 1400 700]);
species_names = {'Ar', 'Ar^+', 'Ar^m', 'Ar^r', 'Ar^p'};
ref_data = {ref.n_Ar, ref.n_Ar_plus, ref.n_Arm, ref.n_Arr, ref.n_Arp};
model_data = {n_Ar_model, n_Ar_plus_model, n_Arm_model, n_Arr_model, n_Arp_model};
errors = {error_Ar, error_Ar_plus, error_Arm, error_Arr, error_Arp};

% Top row: Densities
for i = 1:5
    subplot(2,5,i);
    semilogy(ref.power_species, ref_data{i}, 'ko-', 'LineWidth', 1.5, 'MarkerSize', 4);
    hold on;
    semilogy(ref.power_species, model_data{i}, 's--', 'LineWidth', 1.5, 'MarkerSize', 4, 'Color', [0.85 0.33 0.10]);
    xlabel('Power (W)'); ylabel('Density (m^{-3})');
    title(species_names{i}); grid on;
    if i == 1
        legend('Ref', 'Model', 'Location', 'best');
    end
end

% Bottom row: Errors
for i = 1:5
    subplot(2,5,i+5);
    plot(ref.power_species, errors{i}, 'r-', 'LineWidth', 2);
    xlabel('Power (W)'); ylabel('Error (%)');
    title(sprintf('%s Error (mean=%.1f%%)', species_names{i}, mean(errors{i})));
    yline(10, 'k--'); yline(20, 'k:');
    ylim([0 max(max(errors{i})*1.2, 25)]);
    grid on;
end

sgtitle(sprintf('Species Validation - Rcoil=%.2f Ω', Rcoil_opt), 'FontSize', 14);
saveas(gcf, 'inverse_extended_species.png');

% Summary bar chart
figure('Color', 'w', 'Position', [200 200 800 400]);
all_errors = [mean(error_Ar), mean(error_Ar_plus), mean(error_Arm), mean(error_Arr), ...
              mean(error_Arp), mean(error_Te), mean(error_Tg)];
b = bar(all_errors);
b.FaceColor = 'flat';
colors = parula(7);
for k = 1:7
    b.CData(k,:) = colors(k,:);
end
set(gca, 'XTickLabel', {'Ar', 'Ar^+', 'Ar^m', 'Ar^r', 'Ar^p', 'T_e', 'T_g'});
ylabel('Mean Relative Error (%)');
title(sprintf('Error Summary - Rcoil=%.2f, \\sigma_h=%.1e, \\sigma_i=%.1e', ...
      Rcoil_opt, sigma_h_opt, sigma_i_opt));
yline(5, 'g--', '5%', 'LineWidth', 2);
yline(10, 'r--', '10%', 'LineWidth', 2);
grid on;
saveas(gcf, 'inverse_extended_summary.png');

%% Output
fprintf('\n========================================\n');
fprintf('COPY THESE TO YOUR Ar_script.m:\n');
fprintf('========================================\n');
fprintf('params.Rcoil = %.6f;\n', Rcoil_opt);
fprintf('params.sigma_heating = %.6e;\n', sigma_h_opt);
fprintf('params.sigma_nn = %.6e;\n', sigma_nn_opt);
fprintf('params.sigma_i = %.6e;\n', sigma_i_opt);
fprintf('params.gamma_g = %.6e;\n', gamma_g_opt);
fprintf('params.gamma_m = %.6f;\n', gamma_m_opt);
fprintf('params.gamma_r = %.6f;\n', gamma_r_opt);
fprintf('params.gamma_p = %.6f;\n', gamma_p_opt);

%% Helper functions
function err = run_objective(x, base_params, ref, weights)
    Rcoil = x(1);
    sigma_heating = 10^x(2);
    sigma_nn = 10^x(3);
    sigma_i = 10^x(4);
    gamma_g = 10^x(5);
    gamma_m = 10^x(6);
    gamma_r = 10^x(7);
    gamma_p = 10^x(8);
    
    if Rcoil < 0.1, err = 1e6; return; end
    
    params = base_params;
    params.Rcoil = Rcoil;
    params.sigma_heating = sigma_heating;
    params.sigma_nn = sigma_nn;
    params.sigma_i = sigma_i;
    params.gamma_g = gamma_g;
    params.gamma_m = gamma_m;
    params.gamma_r = gamma_r;
    params.gamma_p = gamma_p;
    
    all_powers = unique([ref.power_species; ref.power_Te; ref.power_Tg]);
    all_powers = sort(all_powers);
    step = max(1, floor(length(all_powers)/12));
    all_powers = all_powers(1:step:end);
    N = length(all_powers);
    
    y0 = [6.4e19; 1e17; 7e16; 1e15; 3e17; 3.6; 310];
    opts = odeset('RelTol', 1e-4, 'AbsTol', [1e15 1e14 1e14 1e13 1e14 1e-3 1], ...
                  'NonNegative', 1:7, 'MaxStep', 0.01);
    
    Y_model = nan(N, 7);
    
    for k = 1:N
        params.PRF = all_powers(k);
        if k == 1
            tspan = [0, 0.2];
        else
            tspan = [0, 0.03];
            y0_new = Y_model(k-1, :)';
            if any(isnan(y0_new)) || any(y0_new <= 0)
                err = 1e6; return;
            end
            y0 = y0_new;
        end
        
        try
            sol = ode15s(@(t,z) rhs_global(t,z,params,[]), tspan, y0, opts);
            if isempty(sol.y) || sol.x(end) < tspan(2)*0.5
                err = 1e6; return;
            end
            Y_model(k,:) = sol.y(:,end)';
        catch
            err = 1e6; return;
        end
    end
    
    if any(isnan(Y_model(:))) || any(isinf(Y_model(:)))
        err = 1e6; return;
    end
    
    n_Ar_model = interp1(all_powers, Y_model(:,1), ref.power_species, 'linear', 'extrap');
    n_Arm_model = interp1(all_powers, Y_model(:,2), ref.power_species, 'linear', 'extrap');
    n_Arr_model = interp1(all_powers, Y_model(:,3), ref.power_species, 'linear', 'extrap');
    n_Arp_model = interp1(all_powers, Y_model(:,4), ref.power_species, 'linear', 'extrap');
    n_Ar_plus_model = interp1(all_powers, Y_model(:,5), ref.power_species, 'linear', 'extrap');
    Te_model = interp1(all_powers, Y_model(:,6), ref.power_Te, 'linear', 'extrap');
    Tg_model = interp1(all_powers, Y_model(:,7), ref.power_Tg, 'linear', 'extrap');
    
    if any(Te_model <= 0) || any(Tg_model <= 0)
        err = 1e6; return;
    end
    
    err_Ar = rms((n_Ar_model - ref.n_Ar) ./ ref.n_Ar) * 100;
    err_Ar_plus = rms((n_Ar_plus_model - ref.n_Ar_plus) ./ ref.n_Ar_plus) * 100;
    err_Arm = rms((n_Arm_model - ref.n_Arm) ./ ref.n_Arm) * 100;
    err_Arr = rms((n_Arr_model - ref.n_Arr) ./ ref.n_Arr) * 100;
    err_Arp = rms((n_Arp_model - ref.n_Arp) ./ ref.n_Arp) * 100;
    err_Te = rms((Te_model - ref.Te) ./ ref.Te) * 100;
    err_Tg = rms((Tg_model - ref.Tg) ./ ref.Tg) * 100;
    
    err = weights.Ar*err_Ar + weights.Ar_plus*err_Ar_plus + ...
          weights.Arm*err_Arm + weights.Arr*err_Arr + weights.Arp*err_Arp + ...
          weights.Te*err_Te + weights.Tg*err_Tg;
    
    if err > 1e4, err = 1e6; end
end

function err = run_objective_bounded(x, base_params, ref, weights, lb, ub)
    penalty = 0;
    for i = 1:length(x)
        if x(i) < lb(i)
            penalty = penalty + 1e4*(lb(i)-x(i))^2;
        elseif x(i) > ub(i)
            penalty = penalty + 1e4*(x(i)-ub(i))^2;
        end
    end
    err = run_objective(x, base_params, ref, weights) + penalty;
end

function [Te_model, Tg_model, n_Ar, n_Ar_plus, n_Arm, n_Arr, n_Arp] = run_full_model(params, ref)
    all_powers = unique([ref.power_species; ref.power_Te; ref.power_Tg]);
    all_powers = sort(all_powers);
    N = length(all_powers);
    
    y0 = [6.4e19; 1e17; 7e16; 1e15; 3e17; 3.6; 310];
    opts = odeset('RelTol', 1e-6, 'AbsTol', [1e14 1e13 1e13 1e12 1e13 1e-4 0.1], 'NonNegative', 1:7);
    
    Y_model = nan(N, 7);
    
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
    end
    
    n_Ar = interp1(all_powers, Y_model(:,1), ref.power_species, 'linear', 'extrap');
    n_Arm = interp1(all_powers, Y_model(:,2), ref.power_species, 'linear', 'extrap');
    n_Arr = interp1(all_powers, Y_model(:,3), ref.power_species, 'linear', 'extrap');
    n_Arp = interp1(all_powers, Y_model(:,4), ref.power_species, 'linear', 'extrap');
    n_Ar_plus = interp1(all_powers, Y_model(:,5), ref.power_species, 'linear', 'extrap');
    Te_model = interp1(all_powers, Y_model(:,6), ref.power_Te, 'linear', 'extrap');
    Tg_model = interp1(all_powers, Y_model(:,7), ref.power_Tg, 'linear', 'extrap');
end