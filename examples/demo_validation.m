% demo_validation.m - Validation example with error comparison plots
%
% This script validates error estimates by:
%   1. Computing stresslet potential with adaptive quadrature
%   2. Comparing to reference solution (stresslet identity or high upsampling)
%   3. Visualizing actual errors alongside estimates
%
% Similar to legacy/Ex_Quick_ErrEstimates.m but using the new quadest API.
%
% Run from the QuadEst directory:
%   >> run('examples/demo_validation.m')
%
% Authors: Pritpal 'Pip' Matharu, David Krantz
% Max Planck Institute for Mathematics in the Sciences / KTH

%% Setup path
scriptDir = fileparts(mfilename('fullpath'));
projectDir = fileparts(scriptDir);
addpath(projectDir);

%% User Parameters
% -------------------------------------------------------------------------
% Geometry selection: 'spheroid', 'peanut', or 'capsule'
geometryShape = 'capsule';

% Density flag:
%   0 = Stresslet identity (constant density) → reference is zero
%   1 = Analytic density → reference via high upsampling
dflag = 1;

% Tolerance for classification
TOL = 1e-6;

% Upsampling factors to precompute (direct quad factor 1 added automatically)
upsamp_fac = 2:4;

% Discretization (number of quadrature points)
nth = 40;   % Theta (Gauss-Legendre) points
nph = 40;   % Phi (trapezoidal) points

% Reference upsampling factor (for dflag = 1)
refFactor = 10;

% Special quadrature placeholder factor (high upsampling)
sqFactor = 5;

% Plot settings
FS = 16;  % Font size
% -------------------------------------------------------------------------

%% Step 1: Setup geometry
fprintf('=== ERROR ESTIMATE VALIDATION ===\n\n');
fprintf('Step 1: Setup geometry\n');

switch lower(geometryShape)
    case 'spheroid'
        geom = quadest.geometry.Spheroid('a', 0.05, 'c', 0.1);
        % Grid bounds for spheroid
        M = 100; L = 0.2;
        xv = linspace(-L/2, L/2, M);
        zv = linspace(-L, L, M);
    case 'peanut'
        geom = quadest.geometry.Peanut();
        % Larger grid for peanut
        M = 100; L = 5;
        xv = linspace(-L/1.5, L/1.5, M);
        zv = linspace(-L, L, M);
    case 'capsule'
        geom = quadest.geometry.Capsule('R', 1, 'L', 6, 'kappa', 3);
        % Grid bounds for capsule: slightly beyond R and L/2
        M = 100;
        xv = linspace(-4*geom.maxRadius(), 4*geom.maxRadius(), M);
        zv = linspace(-1.5*geom.maxHeight(), 1.5*geom.maxHeight(), M);
    otherwise
        error('Unknown geometry: %s', geometryShape);
end
fprintf('  Geometry: %s\n', geometryShape);

%% Step 2: Select kernel
fprintf('\nStep 2: Select kernel\n');

kernel = quadest.kernel.StokesStresslet();
fprintf('  Kernel: %s (p = %.1f)\n', kernel.kernelName(), kernel.singularityOrder());

%% Step 3: Precompute error estimates
fprintf('\nStep 3: Precompute error estimates\n');

tic;
estimator = quadest.errorest.ErrorEstimator(geom, kernel, ...
    'nth', nth, 'nph', nph, ...
    'upsampFactors', upsamp_fac, ...
    'interpolateRoots', true);
fprintf('  Precomputation time: %.2f seconds\n', toc);

grid_surf = estimator.getGrid();

% Get the actual upsampling factors used (includes factor 1 for direct quad)
allUpsampFactors = estimator.upsampFactors;

%% Step 4: Create target grid (xz-plane at y=0)
fprintf('\nStep 4: Create target grid\n');

[X, Z] = meshgrid(xv, zv);
yp = 0.0;
Y = yp * ones(size(X));
targets_all = [X(:), Y(:), Z(:)];

% Filter to exterior points
mask_ext = geom.isExterior(targets_all);
targets = targets_all(mask_ext, :);

fprintf('  Grid: %d × %d = %d total points\n', M, M, M*M);
fprintf('  Exterior points: %d\n', sum(mask_ext));

%% Step 5: Define density
fprintf('\nStep 5: Define density\n');

[theta_mat, phi_mat] = grid_surf.meshgrid();

switch dflag
    case 1  % Analytic density
        sigma1 = 2.1 + sin(8*theta_mat) + sin(12*phi_mat);
        sigma2 = 2 + sin(6*theta_mat) .* cos(6*phi_mat);
        sigma3 = sin(5*theta_mat) .* exp(-cos(phi_mat).^2) + 1.03;
        density = [sigma1(:), sigma2(:), sigma3(:)];
        densityType = 'Analytic';
    otherwise  % Stresslet identity
        Nb = grid_surf.numPoints();
        density = [300*ones(Nb,1), 100*ones(Nb,1), 200*ones(Nb,1)];
        densityType = 'Stresslet Identity (constant)';
end

fprintf('  Density type: %s\n', densityType);

% Check density resolution
quadest.util.Diagnostics.checkDensityResolution(density, nph);

%% Step 6: Classify targets using error estimates
fprintf('\nStep 6: Classify targets\n');

tic;
[estimates, classification] = estimator.evaluate(targets, density, 'tol', TOL);
fprintf('  Classification time: %.2f seconds\n', toc);

% Count targets in each region
fprintf('\n  Classification (tol = %.0e):\n', TOL);
for k = 1:length(allUpsampFactors)
    nPts = sum(classification.masks{k});
    if allUpsampFactors(k) == 1
        fprintf('    Direct quadrature:  %d points\n', nPts);
    else
        fprintf('    Upsampling %d×:      %d points\n', allUpsampFactors(k), nPts);
    end
end
nSQ = sum(classification.isSQ);
fprintf('    Special quadrature: %d points\n', nSQ);

%% Step 7: Compute solution using adaptive quadrature
fprintf('\nStep 7: Compute with adaptive quadrature\n');

u = zeros(size(targets));
tic;

% Direct and upsampling regions (iterate over all factors including 1)
for k = 1:length(allUpsampFactors)
    mask_k = classification.masks{k};
    if any(mask_k)
        u(mask_k, :) = kernel.evaluateUpsampled(targets(mask_k, :), grid_surf, density, allUpsampFactors(k));
    end
end

% Special quadrature region (placeholder: use high upsampling)
if any(classification.isSQ)
    u(classification.isSQ, :) = kernel.evaluateUpsampled(targets(classification.isSQ, :), grid_surf, density, sqFactor);
end

fprintf('  Computation time: %.2f seconds\n', toc);

%% Step 8: Compute reference solution
fprintf('\nStep 8: Compute reference solution\n');

tic;
switch dflag
    case 1  % Analytic density → use high upsampling
        uref = kernel.evaluateUpsampled(targets, grid_surf, density, refFactor);
        fprintf('  Reference: High upsampling (factor = %d)\n', refFactor);
    otherwise  % Stresslet identity → zero solution
        uref = zeros(size(u));
        fprintf('  Reference: Stresslet identity (exact = 0)\n');
end
fprintf('  Reference computation time: %.2f seconds\n', toc);

%% Step 9: Compute errors
fprintf('\nStep 9: Compute errors\n');

% L2 error at each target (adaptive quadrature result)
error_vec = sqrt(sum((u - uref).^2, 2));

% Per-factor errors: compute the potential at ALL exterior targets using
% each upsampling factor separately, then compare to reference
fprintf('  Computing per-factor errors...\n');
errors_cell = cell(1, length(allUpsampFactors));
for k = 1:length(allUpsampFactors)
    u_k = kernel.evaluateUpsampled(targets, grid_surf, density, allUpsampFactors(k));
    errors_cell{k} = sqrt(sum((u_k - uref).^2, 2));
end

% Statistics per region
fprintf('\n  Error statistics:\n');
for k = 1:length(allUpsampFactors)
    if any(classification.masks{k})
        err_k = error_vec(classification.masks{k});
        if allUpsampFactors(k) == 1
            fprintf('    Direct:    max = %.2e, mean = %.2e\n', max(err_k), mean(err_k));
        else
            fprintf('    Upsamp %d×: max = %.2e, mean = %.2e\n', allUpsampFactors(k), max(err_k), mean(err_k));
        end
    end
end
if any(classification.isSQ)
    err_sq = error_vec(classification.isSQ);
    fprintf('    SQ:        max = %.2e, mean = %.2e\n', max(err_sq), mean(err_sq));
end

%% Step 10: Visualization
fprintf('\nStep 10: Generate plots\n');

% Embed solution/error back into full grid for plotting
u_full = zeros(size(targets_all));
u_full(mask_ext, :) = u;
error_full = NaN(size(targets_all, 1), 1);
error_full(mask_ext) = error_vec;

%% Figure 1: Quadrature regions
figure('Name', 'Quadrature Regions');
quadest.util.Plotting.plotQuadRegions(xv, zv, classification, mask_ext, geom, 'FontSize', FS);
title(sprintf('Quadrature Regions (tol = %.0e)', TOL), 'FontSize', FS);

%% Figure 2: Error vs Distance
figure('Name', 'Error vs Distance');
quadest.util.Plotting.plotErrorVsDistance(u, uref, targets, grid_surf, classification, TOL, 'FontSize', FS);

%% Figure 3: Error contour per upsampling factor
% Each panel shows the actual error using that factor + its estimated error
quadest.util.Plotting.plotErrorContour(xv, zv, errors_cell, classification.estimates, ...
    mask_ext, geom, 'upsampFactors', allUpsampFactors, 'levels', 0:-2:-10, 'FontSize', FS);
quadest.util.Plotting.alignfigs;

%% Summary
fprintf('\n=== VALIDATION COMPLETE ===\n');
fprintf('Key results:\n');
fprintf('  - All direct quadrature errors should be < %.0e (tolerance)\n', TOL);
fprintf('  - Error estimate contours should align with actual error level sets\n');
fprintf('  - Error should decay with distance from surface\n');

% Check validation criterion
if any(classification.masks{1})
    maxDirectErr = max(error_vec(classification.masks{1}));
    if maxDirectErr < TOL
        fprintf('\n  ✓ PASS: Max direct quadrature error (%.2e) < tolerance (%.0e)\n', maxDirectErr, TOL);
    else
        fprintf('\n  ✗ FAIL: Max direct quadrature error (%.2e) >= tolerance (%.0e)\n', maxDirectErr, TOL);
    end
end
