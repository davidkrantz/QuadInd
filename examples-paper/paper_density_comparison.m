% paper_density_comparison.m - Compare error estimates with and without density

%% Setup path
clear; close all;
set(groot, 'defaultAxesTickLabelInterpreter', 'latex');
set(groot, 'defaulttextinterpreter', 'latex');
set(groot, 'defaultLegendInterpreter', 'latex');
scriptDir = fileparts(mfilename('fullpath'));
projectDir = fileparts(scriptDir);
addpath(projectDir);

%% Parameters
nth = 80;
nph = 40;
N = nth * nph;  % 3200 points per body

% Geometry
a = 0.1;
c = 0.5;
center2 = [0, 0, 1.0025];  % Center of body 2

% Reference and evaluation
refFactor = 5;

% Target grid
Ngrid = 300;
xv = linspace(-0.06, 0.06, Ngrid);
zv = linspace(0.45, 0.55, Ngrid);

% Plot settings
FS = 16;
levels = 0:-2:-10;
savefig = 0;

fprintf('=== DENSITY COMPARISON FIGURE ===\n\n');

%% Step 1: Load density from data.mat
fprintf('Step 1: Load density\n');

data = load(fullfile(scriptDir, 'density_comparison_data.mat'), 'Q');
Q = data.Q;
assert(numel(Q) == 6 * N, ...
    'Expected Q to have %d elements (6 x %d), got %d', 6*N, N, numel(Q));

% Unpack: Q = [body1_c1; body1_c2; body1_c3; body2_c1; body2_c2; body2_c3]
density1 = [Q(1:N), Q(N+1:2*N), Q(2*N+1:3*N)];
density2 = [Q(3*N+1:4*N), Q(4*N+1:5*N), Q(5*N+1:6*N)];

fprintf('  density1: size = [%d, %d], norm = %.4e\n', size(density1, 1), size(density1, 2), norm(density1(:)));
fprintf('  density2: size = [%d, %d], norm = %.4e\n', size(density2, 1), size(density2, 2), norm(density2(:)));

%% Step 2: Setup geometry, kernel, estimator
fprintf('\nStep 2: Setup geometry and precompute estimator\n');

geom = quadest.geometry.Spheroid('a', a, 'c', c);
kernel = quadest.kernel.StokesStresslet();

tic;
estimator = quadest.errorest.ErrorEstimator(geom, kernel, ...
    'nth', nth, 'nph', nph, ...
    'upsampFactors', 1, ...
    'interpolateRoots', true);
fprintf('  Precomputation time: %.2f s\n', toc);

grid = estimator.getGrid();
assert(grid.numPoints() == N, 'Grid size mismatch: expected %d, got %d', N, grid.numPoints());

%% Step 3: Build target grid
fprintf('\nStep 3: Build target grid\n');

[X, Z] = meshgrid(xv, zv);
targets_all = [X(:), zeros(numel(X), 1), Z(:)];

% Exterior to both bodies
mask_ext1 = geom.isExterior(targets_all);
mask_ext2 = geom.isExterior(targets_all - center2);
mask_ext = mask_ext1 & mask_ext2;

targets_ext = targets_all(mask_ext, :);
targets_ext_shifted = targets_ext - center2;  % In body 2's frame

fprintf('  Grid: %d x %d = %d total points\n', Ngrid, Ngrid, Ngrid^2);
fprintf('  Exterior to both bodies: %d points\n', sum(mask_ext));

%% Step 4: Compute measured quadrature error
fprintf('\nStep 4: Compute measured error (direct vs reference, factor %d)\n', refFactor);

tic;
% Body 1 (grid at origin, targets as-is)
fprintf('  Body 1: direct quadrature...\n');
u1_direct = kernel.evaluateOnGrid(targets_ext, grid, density1);
fprintf('  Body 1: reference (factor %d)...\n', refFactor);
u1_ref = kernel.evaluateUpsampled(targets_ext, grid, density1, refFactor);

% Body 2 (grid at origin, shifted targets)
fprintf('  Body 2: direct quadrature...\n');
u2_direct = kernel.evaluateOnGrid(targets_ext_shifted, grid, density2);
fprintf('  Body 2: reference (factor %d)...\n', refFactor);
u2_ref = kernel.evaluateUpsampled(targets_ext_shifted, grid, density2, refFactor);
fprintf('  Kernel evaluation time: %.2f s\n', toc);
%%
% Per-body error (max over 3 components for scalar per target)
%err1 = max(abs(u1_direct - u1_ref), [], 2);
%err2 = max(abs(u2_direct - u2_ref), [], 2);

err1 = sum((u1_direct-u1_ref).^2,2).^(1/2);
err2 = sum((u2_direct-u2_ref).^2,2).^(1/2);

% Combined: max over both bodies
measured_error = max(err1, err2);

fprintf('  Body 1 error: max = %.2e, median = %.2e\n', max(err1), median(err1));
fprintf('  Body 2 error: max = %.2e, median = %.2e\n', max(err2), median(err2));
fprintf('  Combined error: max = %.2e, median = %.2e\n', max(measured_error), median(measured_error));

%% Step 5: Compute error estimates
fprintf('\nStep 5: Compute error estimates\n');

tic;
% With density
est1_wd = estimator.evaluate(targets_ext, density1);
est2_wd = estimator.evaluate(targets_ext_shifted, density2);
est_with_density = max(est1_wd, est2_wd);

% Without density (unit density -> uniform estimate)
unit_density = ones(N, 3);
est1_nd = estimator.evaluate(targets_ext, unit_density);
est2_nd = estimator.evaluate(targets_ext_shifted, unit_density);
est_without_density = max(est1_nd, est2_nd);
fprintf('  Estimation time: %.2f s\n', toc);

fprintf('  With density:    max = %.2e, median = %.2e\n', max(est_with_density), median(est_with_density));
fprintf('  Without density: max = %.2e, median = %.2e\n', max(est_without_density), median(est_without_density));

%% Step 6: Plot 2-panel figure
close all;
fprintf('\nStep 6: Generate figure\n');

panelData = {est_with_density, est_without_density};

for ip = 1:2
    figure('DefaultAxesFontSize',FS);

    estimates = panelData{ip};

    % Reshape to grid (log10 scale)
    err_full = NaN(numel(mask_ext), 1);
    err_full(mask_ext) = log10(abs(measured_error) + eps);
    Err = reshape(err_full, Ngrid, Ngrid);

    est_full = NaN(numel(mask_ext), 1);
    est_full(mask_ext) = log10(abs(estimates) + eps);
    Est = reshape(est_full, Ngrid, Ngrid);

    % Clamp and clean
    minLevel = levels(end);
    Err(Err < minLevel) = minLevel;
    Est(Est < minLevel) = minLevel;
    Err(isinf(Err)) = -16;
    Est(isinf(Est)) = -16;

    % Plot
    contourf(X, Z, Err, levels, 'LineColor', 'none');
    hold on;
    [C, h] = contour(X, Z, Est, levels(2:end), '-k', 'LineWidth', 2);
    clabel(C, h, 'LabelSpacing', 50000, 'FontSize', 14, 'Color', 'k', 'Interpreter', 'latex');

    if ip == 2
        figure(1);
        [C, h] = contour(fliplr(X),Z,fliplr(Est),levels(2:end), '--', 'Color', [0.2 0.8 0.2], 'LineWidth', 2);
        clabel(C, h, 'LabelSpacing', 50000, 'FontSize', 14, 'Color', [0.2 0.8 0.2], 'Interpreter', 'latex');
        h1 = plot(nan, nan, 'k-',  'LineWidth', 2); hold on;
        h2 = plot(nan, nan, '--', 'Color', [0.2 0.8 0.2], 'LineWidth', 2);
        legend([h1, h2], {'With density modifier', 'Without density modifier'}, ...
            'Interpreter', 'latex', 'Location', 'southwest', 'FontSize', FS);
        figure(2);
    end

    % Colorbar
    cbar = colorbar('Ticks', fliplr(levels), 'TickLabelInterpreter', 'latex', 'FontSize', FS);
    xlabel(cbar, '$\log_{10}(\textrm{Absolute error})$', 'FontSize', FS, 'Interpreter', 'latex');

    % Colormap
    colormap(gca, quadest.util.Plotting.divergingColormap(length(levels) - 1));
    
    % Plot first spheroid
    quadest.util.Plotting.plotSurfaceMesh(grid,'FlipYZ',true);

    % Plot second sheroid (shifted)
    Xtmp = reshape(grid.x(:,1), nth, nph);
    Ytmp = reshape(grid.x(:,2), nth, nph);
    Ztmp = reshape(grid.x(:,3), nth, nph);
    Ztmp = Ztmp + center2(3);
    Xtmp = [Xtmp, Xtmp(:,1)];
    Ytmp = [Ytmp, Ytmp(:,1)];
    Ztmp = [Ztmp, Ztmp(:,1)];
    surf(Xtmp, Ztmp, Ytmp, ...
        'FaceColor', 0.8*[1 1 1], ...
        'EdgeColor', [0.3 0.3 0.3], ...
        'FaceAlpha', 1);
    axis equal;

    view(0, 90);
    clim([levels(end), levels(1)]);
    xlim([min(xv), max(xv)]);
    ylim([min(zv), max(zv)]);

    % Labels
    xlabel('$x$', 'Interpreter', 'latex', 'FontSize', FS);
    ylabel('$z$', 'Interpreter', 'latex', 'FontSize', FS);
    set(gca, 'FontSize', FS);
end

if savefig
    if ~exist('../figs', 'dir'); mkdir('../figs'); end
    disp('saving figures...');
    exportgraphics(figure(1),'../figs/capsule_density_comparison.pdf','Resolution',800);
    disp('sucessfully saved figures');
end

alignfigs;

fprintf('\n=== DONE ===\n');
