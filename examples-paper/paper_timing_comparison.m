% paper_timing_comparison.m - Timing comparison: tabulated vs direct evaluation

%% Setup
clear; close all;
set(groot,'defaultAxesTickLabelInterpreter','latex');
set(groot,'defaulttextinterpreter','latex');
set(groot,'defaultLegendInterpreter','latex');
scriptDir = fileparts(mfilename('fullpath'));
projectDir = fileparts(scriptDir);
addpath(projectDir);

% Parameters
nth = 50;
nph = 50;
upsamp_fac = 1;

% Target counts for the sweep
nTargets = logspace(2,5,20);
nRepeats = 3;  % Timing repeats per measurement (take median)

M_pool = 500;  % Grid resolution for generating the target pool
M_contour = 200;  % Grid resolution for contour plot
ref_upsamp = 8;  % Upsampling factor for reference solution

FS = 16;
savefig = 0;

%% Geometry & kernel
fprintf('=== TIMING COMPARISON: TABULATED vs DIRECT ===\n\n');

geom = quadind.geometry.Capsule('R', 1, 'L', 6, 'kappa', 3);
kernel = quadind.kernel.StokesStresslet();

fprintf('Geometry: Capsule (R = 1, L = 6, kappa = 3)\n');
fprintf('Kernel:   %s (p = %.1f)\n', kernel.kernelName(), kernel.singularityOrder());
fprintf('Grid:     %d x %d = %d nodes\n\n', nth, nph, nth*nph);

%% Precompute error indicators (timed)
fprintf('Precomputing error indicators...\n');
tic;
evaluator = quadind.IndicatorEvaluator(geom, kernel, ...
    'nth', nth, 'nph', nph, ...
    'upsampFactors', upsamp_fac, ...
    'interpolateRoots', true);
t_precomp = toc;
fprintf('  Precomputation time: %.2f s\n', t_precomp);

grid_surf = evaluator.getGrid();
Nb = grid_surf.numPoints();

% Define near-singular density (smooth Poisson-kernel-like near-source)
[theta_mat, phi_mat] = grid_surf.meshgrid();
peakCenter = [0.68, 0.00, 0.0];
sourceRadius = 0.47;
peakScale = 1e4;
sinTheta = sin(theta_mat);
cosTheta = cos(theta_mat);
xhat = sinTheta .* cos(phi_mat);
yhat = sinTheta .* sin(phi_mat);
zhat = cosTheta;
peakCenter = peakCenter / norm(peakCenter);
peakDot = peakCenter(1)*xhat + peakCenter(2)*yhat + peakCenter(3)*zhat;
peakRaw = (1 - sourceRadius^2) ./ ...
    (1 - 2*sourceRadius*peakDot + sourceRadius^2).^(3/2);
peak = peakScale * ...
    (peakRaw - min(peakRaw(:))) / (max(peakRaw(:)) - min(peakRaw(:)));
sigma1 = 1.0 + 14.0*peak .* xhat;
sigma2 = 1.1 + 10.0*peak .* zhat;
sigma3 = 0.9 + 12.0*peak .* yhat;
density = [sigma1(:), sigma2(:), sigma3(:)];

% Create large pool of exterior target points
fprintf('\nCreating target pool (%d x %d grid_surf)...\n', M_pool, M_pool);
xv = linspace(-4.0, 4.5, M_pool);
zv = linspace(-4.5, 4.5, M_pool);
[X, Z] = meshgrid(xv, zv);
targets_pool = [X(:), zeros(M_pool^2, 1), Z(:)];

mask_ext = geom.isExterior(targets_pool);
targets_pool = targets_pool(mask_ext, :);
nPool = size(targets_pool, 1);

% Create points for contour plot
fprintf('\nCreating contour points (%d x %d grid_est)...\n', M_contour, M_contour);
xv = linspace(-3.0, 3.0, M_contour);
zv = linspace(-4.5, 4.5, M_contour);
[X, Z] = meshgrid(xv, zv);
contour_points = [X(:), zeros(M_contour^2, 1), Z(:)];
mask_ext_contour = geom.isExterior(contour_points);

fprintf('  Pool size: %d exterior targets\n', nPool);

%% Timing sweep
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
        est_tab = evaluator.evaluate(sub_targets, density);
        times_tab(r) = toc;
    end
    t_tabulated(i) = median(times_tab);

    % Direct evaluation
    times_dir = zeros(nRepeats, 1);
    for r = 1:nRepeats
        tic;
        est_dir = evaluator.evaluateDirect(sub_targets, density);
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

%% Compute error at all targets for final tabulated and direct indicators (for plotting)
fprintf('Computing errors for final tabulated and direct indicators...\n');
est_tab = evaluator.evaluate(contour_points(mask_ext_contour, :), density);
est_dir = evaluator.evaluateDirect(contour_points(mask_ext_contour, :), density);
u_dir = kernel.evaluateOnGrid(contour_points(mask_ext_contour, :), grid_surf, density);
u_ref = kernel.evaluateUpsampled(contour_points(mask_ext_contour, :), grid_surf, density, ref_upsamp);
u_err = vecnorm(u_dir - u_ref, 2, 2);

%% Plot: log-log timing comparison
close all;

figure;
figure('DefaultAxesFontSize',FS);

% Tabulated (blue circles + line)
loglog(allCounts, t_tabulated, 'o-', 'Color', [0.0, 0.45, 0.74], ...
    'MarkerFaceColor', [0.0, 0.45, 0.74], 'LineWidth', 1.5, 'MarkerSize', 7);
hold on;

% Precomputation time (horizontal dashed black line)
xlims = [min(allCounts), max(allCounts)];

% Reference slope (dotted black line)
loglog(allCounts,1e-6*allCounts,'k:','LineWidth',1.5);

xlim(xlims);
xticks([1e1 1e2 1e3 1e4 1e5]);
xlabel('Number of targets, $M$', 'FontSize', FS,'Interpreter', 'latex');
ylabel('Time (s)', 'FontSize', FS,'Interpreter', 'latex');
grid on;
annotation('textarrow',[0.44 0.5],[0.6 0.53],'String','$\mathcal{O}(M)$','fontsize',FS,'interpreter','latex')
set(gca, 'FontSize', FS);

% Plot contour of error for tabulated and direct indicators
figure('DefaultAxesFontSize',FS);
quadind.util.Plotting.plotErrorContour(xv, zv, u_err, est_tab, mask_ext_contour, geom, 'showLabels', false);
figure('DefaultAxesFontSize',FS);
quadind.util.Plotting.plotErrorContour(xv, zv, u_err, est_dir, mask_ext_contour, geom , 'showLabels', false);

close(1);
alignfigs;

if savefig
    if ~exist('../figs', 'dir'); mkdir('../figs'); end
    disp('saving figures...');
    exportgraphics(figure(2),'../figs/capsule_timing_comparison.pdf','Resolution',400);
    exportgraphics(figure(3),'../figs/capsule_contour_tabulated.pdf','Resolution',1000);
    disp('sucessfully saved figures');
end
