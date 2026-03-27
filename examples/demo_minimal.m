% demo_minimal.m - Minimal example of quadest package usage
%
% This script demonstrates the basic workflow:
%   1. Define geometry
%   2. Select kernel
%   3. Precompute error estimates
%   4. Define density and targets
%   5. Evaluate error estimates
%
% Run from the QuadEst directory:
%   >> run('examples/demo_minimal.m')

%% Setup path
% Add the parent directory to access the +quadest package
scriptDir = fileparts(mfilename('fullpath'));
projectDir = fileparts(scriptDir);
addpath(projectDir);

%% Step 1: Define geometry
% Create a prolate spheroid with semi-axes a (xy-plane) and c (z-axis)
geom = quadest.geometry.Spheroid('a', 0.05, 'c', 0.1);
fprintf('Geometry: ');
disp(geom);

%% Step 2: Select kernel
% Use the Stokes stresslet (double-layer) kernel
kernel = quadest.kernel.StokesStresslet();
fprintf('Kernel: %s (singularity order p = %.1f)\n\n', ...
    kernel.kernelName(), kernel.singularityOrder());

%% Step 3: Precompute error estimates
% This builds gridded interpolants for various upsampling factors
fprintf('Precomputing error estimates...\n');
tic;
estimator = quadest.errorest.ErrorEstimator(geom, kernel, ...
    'nth', 40, 'nph', 60, ...           % quadrature grid resolution
    'upsampFactors', 1:2, ...           % upsampling levels to precompute
    'interpolateRoots', true);          % cache theta roots
precompTime = toc;
fprintf('Precomputation time: %.2f seconds\n\n', precompTime);

%% Step 4: Define density and target points
% Use stresslet identity density (constant)
grid = estimator.getGrid();
density = ones(grid.numPoints(), 3);  % [N×3] array

% Single target point outside the spheroid
target = [0.1, 0, 0];
fprintf('Target point: [%.2f, %.2f, %.2f]\n', target);

%% Step 5: Evaluate error estimate
tic;
estimate = estimator.evaluate(target, density);
evalTime = toc;

fprintf('Error estimate: %.2e\n', estimate);
fprintf('Evaluation time: %.4f seconds\n\n', evalTime);

%% With classification
% If we provide a tolerance, we get classification information
tol = 1e-6;
[estimate, classification] = estimator.evaluate(target, density, 'tol', tol);

fprintf('With tolerance %.0e:\n', tol);
fprintf('  Error estimate: %.2e\n', estimate);
fprintf('  Recommended upsampling: %d\n', classification.upsampfac);
fprintf('  Needs special quadrature: %s\n', mat2str(classification.isSQ));

%% Multiple targets
fprintf('\n--- Multiple targets ---\n');
targets = [
    0.1, 0, 0;       % Far from surface
    0.06, 0, 0;      % Closer
    0.052, 0, 0;     % Very close
    0, 0, 0.15;      % Along z-axis
];

estimates = estimator.evaluate(targets, density);
[~, class] = estimator.evaluate(targets, density, 'tol', tol);

fprintf('Target                 Estimate    Upsampling   SQ?\n');
fprintf('------------------------------------------------------\n');
for i = 1:size(targets, 1)
    fprintf('[%5.3f, %5.3f, %5.3f]   %.2e      %d          %s\n', ...
        targets(i,:), estimates(i), class.upsampfac(i), mat2str(class.isSQ(i)));
end

fprintf('\nDone!\n');
