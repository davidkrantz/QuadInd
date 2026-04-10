% paper_timing_comparison.m - Timing comparison: tabulated vs direct evaluation paper

% Setup
clear; close all;
set(groot,'defaultAxesTickLabelInterpreter','latex');
set(groot,'defaulttextinterpreter','latex');
set(groot,'defaultLegendInterpreter','latex');
scriptDir = fileparts(mfilename('fullpath'));
projectDir = fileparts(scriptDir);
addpath(projectDir);

% Parameters
nth = 40;
nph = 40;
upsamp_fac = 1;
savefig = 0;

% Target counts for the sweep
nTargets = logspace(1,5,20);
nRepeats = 3;  % Timing repeats per measurement (take median)

M_pool = 500;  % Grid resolution for generating the target pool
FS = 16;

% Geometry & kernel
fprintf('=== TIMING COMPARISON: TABULATED vs DIRECT ===\n\n');

geom = quadest.geometry.Capsule('R', 1, 'L', 6, 'kappa', 3);
%geom = quadest.geometry.Spheroid('a', 0.05, 'c', 0.1);
kernel = quadest.kernel.StokesStresslet();

fprintf('Geometry: Capsule (R = 1, L = 6, kappa = 3)\n');
fprintf('Kernel:   %s (p = %.1f)\n', kernel.kernelName(), kernel.singularityOrder());
fprintf('Grid:     %d x %d = %d nodes\n\n', nth, nph, nth*nph);

% Precompute error estimates (timed)
fprintf('Precomputing error estimates...\n');
tic;
estimator = quadest.errorest.ErrorEstimator(geom, kernel, ...
    'nth', nth, 'nph', nph, ...
    'upsampFactors', upsamp_fac, ...
    'interpolateRoots', true);
t_precomp = toc;
fprintf('  Precomputation time: %.2f s\n', t_precomp);

grid_est = estimator.getGrid();
Nb = grid_est.numPoints();

% Define oscillatory density
[theta_mat, phi_mat] = grid_est.meshgrid();

sigma1 = 2.1 + sin(8*theta_mat) + sin(12*phi_mat);
sigma2 = 2 + sin(6*theta_mat) .* cos(6*phi_mat);
sigma3 = sin(5*theta_mat) .* exp(-cos(phi_mat).^2) + 1.03;
density = [sigma1(:), sigma2(:), sigma3(:)];

% Create large pool of exterior target points
fprintf('\nCreating target pool (%d x %d grid_est)...\n', M_pool, M_pool);
xv = linspace(-1.5 * geom.maxRadius(), 1.5 * geom.maxRadius(), M_pool);
zv = linspace(-1.5 * geom.maxHeight(), 1.5 * geom.maxHeight(), M_pool);
[X, Z] = meshgrid(xv, zv);
targets_pool = [X(:), zeros(M_pool^2, 1), Z(:)];

mask_ext = geom.isExterior(targets_pool);
targets_pool = targets_pool(mask_ext, :);
nPool = size(targets_pool, 1);

fprintf('  Pool size: %d exterior targets\n', nPool);

% Timing sweep
allCounts = nTargets;
allCounts = allCounts(allCounts <= nPool);

t_tabulated = NaN(length(allCounts), 1);
t_direct    = NaN(length(allCounts), 1);

fprintf('\n--- Timing sweep ---\n');
fprintf('  %-8s  %12s  %12s  %8s\n', 'Targets', 'Tabulated', 'Direct', 'Speedup');
fprintf('  %s\n', repmat('-', 1, 48));

for i = 1:length(allCounts)
    n = allCounts(i);
    idx = round(linspace(1, nPool, n));
    sub_targets = targets_pool(idx, :);

    % Tabulated evaluation (always measured)
    times_tab = zeros(nRepeats, 1);
    for r = 1:nRepeats
        tic;
        est_tab = estimator.evaluate(sub_targets, density);
        times_tab(r) = toc;
    end
    t_tabulated(i) = median(times_tab);

    % Direct evaluation
    times_dir = zeros(nRepeats, 1);
    for r = 1:nRepeats
        tic;
        est_dir = estimator.evaluateDirect(sub_targets, density);
        times_dir(r) = toc;
    end
    t_direct(i) = median(times_dir);

    speedup = t_direct(i) / t_tabulated(i);
    fprintf('  %-8d  %10.4f s  %10.4f s  %7.0fx\n', ...
        n, t_tabulated(i), t_direct(i), speedup);
end

% Compute per-target costs and summarize
mask_valid_both = ~isnan(t_direct);
if any(mask_valid_both)
    cost_tab_per = median(t_tabulated(mask_valid_both) ./ allCounts(mask_valid_both).');
    cost_dir_per = median(t_direct(mask_valid_both) ./ allCounts(mask_valid_both).');
    fprintf('\n  Per-target cost (median):\n');
    fprintf('    Tabulated: %.2f us/target  (%.0f targets/s)\n', 1e6 * cost_tab_per, 1 / cost_tab_per);
    fprintf('    Direct:    %.2f us/target  (%.0f targets/s)\n', 1e6 * cost_dir_per, 1 / cost_dir_per);
    fprintf('    Speedup:   %.0fx\n', cost_dir_per / cost_tab_per);
end

fprintf('\n=== TIMING COMPARISON COMPLETE ===\n');
fprintf('  Precomputation: %.2f s\n', t_precomp);
if any(mask_valid_both)
    fprintf('  Speedup (per-target): ~%.0fx\n', cost_dir_per / cost_tab_per);
end
fprintf('\n');

% Plot: log-log timing comparison
close all;
figure('DefaultAxesFontSize',FS);

% Tabulated (blue circles + line)
loglog(allCounts, t_tabulated, 'o-', 'Color', [0.0, 0.45, 0.74], ...
    'MarkerFaceColor', [0.0, 0.45, 0.74], 'LineWidth', 1.5, 'MarkerSize', 7);
hold on;

% Direct (red squares + line)
mask_dir = ~isnan(t_direct);
loglog(allCounts(mask_dir), t_direct(mask_dir), 's-', 'Color', [0.85, 0.33, 0.1], ...
    'MarkerFaceColor', [0.85, 0.33, 0.1], 'LineWidth', 1.5, 'MarkerSize', 7);

% Precomputation time (horizontal dashed black line)
xlims = [min(allCounts), max(allCounts)];
plot(xlims, [t_precomp, t_precomp], 'k--', 'LineWidth', 1.5);

% Reference slope (dotted black line)
loglog(allCounts,1e-4*allCounts,'k:','LineWidth',1.5);

xlim(xlims);
xticks([1e1 1e2 1e3 1e4 1e5]);
yticks([1e-4 1e-3 1e-2 1e-1 1e0 1e1 1e2 1e3 1e4]);
xlabel('Number of targets, $N$', 'FontSize', FS,'Interpreter', 'latex');
ylabel('Time (s)', 'FontSize', FS,'Interpreter', 'latex');
legend({'Tabulated', 'Direct', 'Precomputation'}, ...
    'Interpreter', 'latex', 'Location', 'northwest', 'FontSize', 14);
grid on;
annotation('textarrow',[0.744 0.65],[0.55 0.615],'String','$\mathcal{O}(N)$','fontsize',FS,'interpreter','latex')
set(gca, 'FontSize', FS);

if savefig
    if ~exist('../figs', 'dir'); mkdir('../figs'); end
    disp('saving figures...');
    exportgraphics(figure(1),'../figs/capsule_timing_comparison.pdf','Resolution',400);
    disp('sucessfully saved figures');
end