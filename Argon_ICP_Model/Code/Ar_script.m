%% generate_nn_training_data.m
% Generate training data for NN surrogate model
% Sweep: Power x Pressure
% Outputs: Densities (nAr, nArm, nArr, nArp, ne) + Tg
% Te excluded due to known 0D model limitations (~3-5% error)

clear; clc; close all;

fprintf('========================================\n');
fprintf('NN SURROGATE TRAINING DATA GENERATION\n');
fprintf('Power x Pressure Sweep\n');
fprintf('========================================\n\n');

%% Step 1: Define Parameter Ranges
% Grid: 100 x 65 = 6500 points (expect ~5000+ successful)
% Restricted to stable pressure range (1-6 mTorr)
P_min = 200;   P_max = 1600;  N_power = 100;
p_min = 1;     p_max = 6;     N_pressure = 65;  % mTorr - restricted for stability

power_vec = linspace(P_min, P_max, N_power);
pressure_vec = linspace(p_min, p_max, N_pressure);  % mTorr

% Convert pressure to Pa (1 mTorr = 0.133322 Pa)
pressure_Pa_vec = pressure_vec * 0.133322;

N_total = N_power * N_pressure;
fprintf('Power: %d points (%.0f - %.0f W)\n', N_power, P_min, P_max);
fprintf('Pressure: %d points (%.1f - %.1f mTorr)\n', N_pressure, p_min, p_max);
fprintf('Total training points: %d\n\n', N_total);

%% Step 2: Setup Model Parameters (OPTIMIZED VALUES)
% These are from the n_is-fixed model optimization
% Validation errors: All outputs < 5% mean error at 2 mTorr

geom.R = 0.06; 
geom.L = 0.10; 
geom.Tg0 = 300;

base_params.geom = geom;
base_params.Vgrid = 1000;
base_params.beta_i = 0.7;
base_params.beta_g = 0.3;
base_params.Q0 = 4.5e18;
base_params.s_sheath = 1.5e-3;
base_params.ps = struct();

% RF and circuit parameters
base_params.f_RF = 13.56e6;
base_params.omega = 2*pi*base_params.f_RF;
base_params.Nturns = 5;
base_params.Rc = 0.07;
base_params.lc = geom.L;
base_params.f_sheath = 0.55;

% Power model: use circuit model
base_params.power_model = 'circuit';

% ========================================
% OPTIMIZED PARAMETERS (from inverse optimization)
% Validated at 2 mTorr with < 5% error on all outputs
% ========================================
base_params.Rcoil = 0.543183;
base_params.sigma_heating = 4.852125e-18;
base_params.sigma_nn = 8.236909e-19;
base_params.sigma_i = 7.446966e-19;
base_params.gamma_g = 1.000000e-04;
base_params.gamma_m = 0.001110;
base_params.gamma_r = 0.001945;
base_params.gamma_p = 0.005982;

% Reference pressure for scaling
base_params.p_ref = 0.266645;  % 2 mTorr in Pa

%% Step 3: Preallocate Results
results = struct();
results.Power = zeros(N_total, 1);
results.Pressure_mTorr = zeros(N_total, 1);
results.Pressure_Pa = zeros(N_total, 1);
results.nAr = zeros(N_total, 1);
results.nArm = zeros(N_total, 1);
results.nArr = zeros(N_total, 1);
results.nArp = zeros(N_total, 1);
results.ne = zeros(N_total, 1);
results.Te_model = zeros(N_total, 1);  % Keep for reference (not for NN training)
results.Tg = zeros(N_total, 1);
results.converged = zeros(N_total, 1);

%% Step 4: Run Parameter Sweep
opts = odeset('RelTol', 1e-6, ...
              'AbsTol', [1e14 1e13 1e13 1e12 1e13 1e-4 0.1], ...
              'NonNegative', 1:7);

fprintf('Running simulations...\n');
fprintf('Progress: ');

idx = 0;
tic;

for i = 1:N_pressure
    % Set pressure for this sweep
    params = base_params;
    params.p0 = pressure_Pa_vec(i);
    
    % Initial condition (reset for each pressure)
    % Scale initial neutral density with pressure using ideal gas law
    kB = 1.380649e-23;
    n_g0 = params.p0 / (kB * geom.Tg0);
    
    % Initial state: [nAr, nArm, nArr, nArp, ne, Te, Tg]
    y0 = [n_g0;           % nAr: most neutrals are ground state
          n_g0 * 1e-3;    % nArm: small metastable fraction
          n_g0 * 1e-3;    % nArr: small resonant fraction  
          n_g0 * 1e-5;    % nArp: very small 4p fraction
          n_g0 * 5e-3;    % ne: ionization fraction ~0.5%
          3.5;            % Te: typical electron temperature (eV)
          310];           % Tg: slightly above wall temperature (K)
    
    for j = 1:N_power
        idx = idx + 1;
        params.PRF = power_vec(j);
        
        % Time span (longer for first point, shorter for continuation)
        if j == 1
            tspan = [0, 1.0];  % Longer to reach steady state at new pressure
        else
            tspan = [0, 0.2];  % Shorter for continuation from previous point
        end
        
        try
            sol = ode15s(@(t,z) rhs_global(t,z,params,[]), tspan, y0, opts);
            y_final = sol.y(:, end);
            
            % Check for valid results
            if any(isnan(y_final)) || any(y_final < 0) || y_final(6) > 20 || y_final(7) > 1000
                error('Invalid solution values');
            end
            
            % Store results
            results.Power(idx) = power_vec(j);
            results.Pressure_mTorr(idx) = pressure_vec(i);
            results.Pressure_Pa(idx) = pressure_Pa_vec(i);
            results.nAr(idx) = y_final(1);
            results.nArm(idx) = y_final(2);
            results.nArr(idx) = y_final(3);
            results.nArp(idx) = y_final(4);
            results.ne(idx) = y_final(5);
            results.Te_model(idx) = y_final(6);
            results.Tg(idx) = y_final(7);
            results.converged(idx) = 1;
            
            % Use final state as IC for next power point
            y0 = y_final;
            
        catch ME
            % Mark as failed
            results.Power(idx) = power_vec(j);
            results.Pressure_mTorr(idx) = pressure_vec(i);
            results.Pressure_Pa(idx) = pressure_Pa_vec(i);
            results.converged(idx) = 0;
            fprintf('\nFailed at P=%.0fW, p=%.1f mTorr: %s\n', ...
                    power_vec(j), pressure_vec(i), ME.message);
            
            % Try to continue with reset IC
            y0 = [n_g0; n_g0*1e-3; n_g0*1e-3; n_g0*1e-5; n_g0*5e-3; 3.5; 310];
        end
    end
    
    % Progress indicator
    fprintf('%.0f%% ', i/N_pressure*100);
end

elapsed = toc;
fprintf('\nDone! Elapsed time: %.1f seconds (%.3f s per point)\n\n', elapsed, elapsed/N_total);

%% Step 5: Check Convergence
n_converged = sum(results.converged);
n_failed = N_total - n_converged;
fprintf('Convergence: %d/%d (%.1f%%) succeeded\n', n_converged, N_total, n_converged/N_total*100);

if n_failed > 0
    fprintf('WARNING: %d points failed to converge\n', n_failed);
    
    % Show which pressure/power combinations failed
    failed_idx = find(results.converged == 0);
    if length(failed_idx) <= 20
        fprintf('Failed points:\n');
        for fi = 1:length(failed_idx)
            fprintf('  P=%.0fW, p=%.1f mTorr\n', ...
                    results.Power(failed_idx(fi)), results.Pressure_mTorr(failed_idx(fi)));
        end
    end
end
fprintf('\n');

%% Step 6: Filter to Converged Points Only
mask = results.converged == 1;
clean_results = struct();
clean_results.Power = results.Power(mask);
clean_results.Pressure_mTorr = results.Pressure_mTorr(mask);
clean_results.Pressure_Pa = results.Pressure_Pa(mask);
clean_results.nAr = results.nAr(mask);
clean_results.nArm = results.nArm(mask);
clean_results.nArr = results.nArr(mask);
clean_results.nArp = results.nArp(mask);
clean_results.ne = results.ne(mask);
clean_results.Te_model = results.Te_model(mask);
clean_results.Tg = results.Tg(mask);

%% Step 7: Save to CSV (for Python/sklearn)
% NN Inputs: Power, Pressure
% NN Outputs: nAr, nArm, nArr, nArp, ne, Tg (6 outputs)
% Note: Te_model included for reference but NOT recommended for training

T = table(clean_results.Power, clean_results.Pressure_mTorr, clean_results.Pressure_Pa, ...
          clean_results.nAr, clean_results.nArm, clean_results.nArr, clean_results.nArp, ...
          clean_results.ne, clean_results.Tg, clean_results.Te_model, ...
          'VariableNames', {'Power_W', 'Pressure_mTorr', 'Pressure_Pa', ...
                           'nAr', 'nArm', 'nArr', 'nArp', 'ne', 'Tg', 'Te_model'});

writetable(T, 'nn_training_data.csv');
fprintf('Training data saved to nn_training_data.csv\n');

%% Step 8: Save to MAT file (for MATLAB use)
save('nn_training_data.mat', 'clean_results', 'power_vec', 'pressure_vec', 'base_params');
fprintf('Training data saved to nn_training_data.mat\n\n');

%% Step 9: Generate Statistics
fprintf('========================================\n');
fprintf('DATA STATISTICS\n');
fprintf('========================================\n');
fprintf('%-12s | %12s | %12s | %12s\n', 'Variable', 'Min', 'Max', 'Mean');
fprintf('-------------|--------------|--------------|-------------\n');
fprintf('%-12s | %12.0f | %12.0f | %12.0f\n', 'Power (W)', min(clean_results.Power), max(clean_results.Power), mean(clean_results.Power));
fprintf('%-12s | %12.2f | %12.2f | %12.2f\n', 'Pressure', min(clean_results.Pressure_mTorr), max(clean_results.Pressure_mTorr), mean(clean_results.Pressure_mTorr));
fprintf('%-12s | %12.2e | %12.2e | %12.2e\n', 'nAr', min(clean_results.nAr), max(clean_results.nAr), mean(clean_results.nAr));
fprintf('%-12s | %12.2e | %12.2e | %12.2e\n', 'nArm', min(clean_results.nArm), max(clean_results.nArm), mean(clean_results.nArm));
fprintf('%-12s | %12.2e | %12.2e | %12.2e\n', 'nArr', min(clean_results.nArr), max(clean_results.nArr), mean(clean_results.nArr));
fprintf('%-12s | %12.2e | %12.2e | %12.2e\n', 'nArp', min(clean_results.nArp), max(clean_results.nArp), mean(clean_results.nArp));
fprintf('%-12s | %12.2e | %12.2e | %12.2e\n', 'ne', min(clean_results.ne), max(clean_results.ne), mean(clean_results.ne));
fprintf('%-12s | %12.2f | %12.2f | %12.2f\n', 'Tg (K)', min(clean_results.Tg), max(clean_results.Tg), mean(clean_results.Tg));
fprintf('%-12s | %12.3f | %12.3f | %12.3f\n', 'Te (eV)', min(clean_results.Te_model), max(clean_results.Te_model), mean(clean_results.Te_model));

%% Step 10: Visualizations
figure('Color', 'w', 'Position', [100 100 1400 900]);

% Define pressures to plot
pressures_to_plot = [1, 2, 5, 10];  % mTorr
colors = lines(length(pressures_to_plot));

% Plot 1: ne vs Power at different pressures
subplot(2,3,1);
hold on;
for k = 1:length(pressures_to_plot)
    p_target = pressures_to_plot(k);
    mask_p = abs(clean_results.Pressure_mTorr - p_target) < 0.5;
    if any(mask_p)
        semilogy(clean_results.Power(mask_p), clean_results.ne(mask_p), '-', ...
                 'Color', colors(k,:), 'LineWidth', 1.5, ...
                 'DisplayName', sprintf('%.0f mTorr', p_target));
    end
end
xlabel('Power (W)');
ylabel('n_e (m^{-3})');
title('Electron Density');
legend('Location', 'best');
grid on;

% Plot 2: Tg vs Power at different pressures
subplot(2,3,2);
hold on;
for k = 1:length(pressures_to_plot)
    p_target = pressures_to_plot(k);
    mask_p = abs(clean_results.Pressure_mTorr - p_target) < 0.5;
    if any(mask_p)
        plot(clean_results.Power(mask_p), clean_results.Tg(mask_p), '-', ...
             'Color', colors(k,:), 'LineWidth', 1.5, ...
             'DisplayName', sprintf('%.0f mTorr', p_target));
    end
end
xlabel('Power (W)');
ylabel('T_g (K)');
title('Gas Temperature');
legend('Location', 'best');
grid on;

% Plot 3: Te vs Power (for reference - not NN target)
subplot(2,3,3);
hold on;
for k = 1:length(pressures_to_plot)
    p_target = pressures_to_plot(k);
    mask_p = abs(clean_results.Pressure_mTorr - p_target) < 0.5;
    if any(mask_p)
        plot(clean_results.Power(mask_p), clean_results.Te_model(mask_p), '-', ...
             'Color', colors(k,:), 'LineWidth', 1.5, ...
             'DisplayName', sprintf('%.0f mTorr', p_target));
    end
end
xlabel('Power (W)');
ylabel('T_e (eV)');
title('Electron Temperature (Reference Only)');
legend('Location', 'best');
grid on;

% Plot 4: 2D surface - ne
subplot(2,3,4);
[P_grid, p_grid] = meshgrid(power_vec, pressure_vec);
ne_grid = griddata(clean_results.Power, clean_results.Pressure_mTorr, ...
                   log10(clean_results.ne), P_grid, p_grid);
surf(P_grid, p_grid, ne_grid, 'EdgeColor', 'none');
xlabel('Power (W)');
ylabel('Pressure (mTorr)');
zlabel('log_{10}(n_e)');
title('Electron Density Surface');
colorbar;
view(45, 30);

% Plot 5: 2D surface - Tg
subplot(2,3,5);
Tg_grid = griddata(clean_results.Power, clean_results.Pressure_mTorr, ...
                   clean_results.Tg, P_grid, p_grid);
surf(P_grid, p_grid, Tg_grid, 'EdgeColor', 'none');
xlabel('Power (W)');
ylabel('Pressure (mTorr)');
zlabel('T_g (K)');
title('Gas Temperature Surface');
colorbar;
view(45, 30);

% Plot 6: Data coverage with convergence
subplot(2,3,6);
scatter(clean_results.Power, clean_results.Pressure_mTorr, 30, clean_results.ne, 'filled');
xlabel('Power (W)');
ylabel('Pressure (mTorr)');
title(sprintf('Training Data Coverage (%d points)', n_converged));
colorbar;
colormap(gca, 'parula');
c = colorbar;
c.Label.String = 'n_e (m^{-3})';

sgtitle(sprintf('NN Training Data: %d Points (%.0f-%.0f W, %.0f-%.0f mTorr)', ...
        n_converged, P_min, P_max, p_min, p_max), 'FontSize', 14, 'FontWeight', 'bold');

saveas(gcf, 'nn_training_data_visualization.png');
fprintf('\nVisualization saved to nn_training_data_visualization.png\n');

%% Step 11: Species plots
figure('Color', 'w', 'Position', [150 150 1200 500]);

species = {'nAr', 'nArm', 'nArr', 'nArp', 'ne'};
species_labels = {'Ar', 'Ar^m', 'Ar^r', 'Ar^p', 'n_e'};
species_data = {clean_results.nAr, clean_results.nArm, clean_results.nArr, ...
                clean_results.nArp, clean_results.ne};

for s = 1:5
    subplot(1,5,s);
    hold on;
    for k = 1:length(pressures_to_plot)
        p_target = pressures_to_plot(k);
        mask_p = abs(clean_results.Pressure_mTorr - p_target) < 0.5;
        if any(mask_p)
            semilogy(clean_results.Power(mask_p), species_data{s}(mask_p), '-', ...
                     'Color', colors(k,:), 'LineWidth', 1.5);
        end
    end
    xlabel('Power (W)');
    ylabel('Density (m^{-3})');
    title(species_labels{s});
    grid on;
end
legend(arrayfun(@(x) sprintf('%d mTorr', x), pressures_to_plot, 'UniformOutput', false), ...
       'Location', 'best');
sgtitle('Species Densities vs Power at Different Pressures', 'FontSize', 12);
saveas(gcf, 'nn_training_species.png');

%% Summary
fprintf('\n========================================\n');
fprintf('SUMMARY\n');
fprintf('========================================\n');
fprintf('Training points generated: %d / %d (%.1f%%)\n', n_converged, N_total, n_converged/N_total*100);
fprintf('\nInput features (2):\n');
fprintf('  - Power: %.0f - %.0f W\n', P_min, P_max);
fprintf('  - Pressure: %.1f - %.1f mTorr\n', p_min, p_max);
fprintf('\nOutput targets (6) - for NN training:\n');
fprintf('  - nAr   (ground state Ar density)\n');
fprintf('  - nArm  (metastable Ar density)\n');
fprintf('  - nArr  (resonant Ar density)\n');
fprintf('  - nArp  (4p state Ar density)\n');
fprintf('  - ne    (electron density)\n');
fprintf('  - Tg    (gas temperature)\n');
fprintf('\nNOTE: Te_model is included for reference but\n');
fprintf('should NOT be used as NN target due to ~3-5%%\n');
fprintf('overestimation in the 0D global model.\n');
fprintf('========================================\n');

fprintf('\nFiles created:\n');
fprintf('  - nn_training_data.csv (for Python)\n');
fprintf('  - nn_training_data.mat (for MATLAB)\n');
fprintf('  - nn_training_data_visualization.png\n');
fprintf('  - nn_training_species.png\n');