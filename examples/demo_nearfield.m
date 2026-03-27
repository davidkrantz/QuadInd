% demo_nearfield.m - Near-field error estimation behavior
%
% This script demonstrates error estimation for targets at varying
% distances from the surface, focusing on the near-singular regime.
%
% Run from the QuadEst directory:
%   >> run('examples/demo_nearfield.m')

%% Setup path
scriptDir = fileparts(mfilename('fullpath'));
projectDir = fileparts(scriptDir);
addpath(projectDir);

%% Setup
fprintf('=== NEAR-FIELD DEMONSTRATION ===\n\n');

% Geometry: spheroid
geom = quadest.geometry.Spheroid('a', 0.05, 'c', 0.1);
kernel = quadest.kernel.StokesStresslet();

fprintf('Geometry: ');
disp(geom);

% Precompute
estimator = quadest.errorest.ErrorEstimator(geom, kernel, ...
    'nth', 40, 'nph', 60, ...
    'upsampFactors', 1:6, ...
    'interpolateRoots', true);

grid_est = estimator.getGrid();
density = ones(grid_est.numPoints(), 3);  % Stresslet identity

%% Create targets at varying distances
% Start from surface and move outward in x-direction
a = 0.05;  % Semi-axis
distances = logspace(-4, 0, 50);  % 0.0001 to 1
targets = [(a + distances)', zeros(length(distances), 2)];

fprintf('Testing %d targets at distances from %.4f to %.2f\n\n', ...
    length(distances), min(distances), max(distances));

%% Evaluate estimates for each upsampling factor
tol = 1e-6;
[estimates, classification] = estimator.evaluate(targets, density, 'tol', tol);

%% Plot results
figure('Position', [100, 100, 1000, 400]);

% Plot 1: Error estimate vs distance
subplot(1, 2, 1);
loglog(distances, estimates, 'b-', 'LineWidth', 2);
hold on;
yline(tol, 'r--', 'LineWidth', 1.5, 'Label', sprintf('tol = %.0e', tol));
xlabel('Distance from surface');
ylabel('Error estimate');
title('Error Estimate vs Distance');
grid on;
legend('Estimate', 'Tolerance', 'Location', 'northeast');

% Plot 2: Recommended upsampling vs distance
subplot(1, 2, 2);
upsamp = classification.upsampfac;
upsamp(classification.isSQ) = 7;  % Mark SQ as 7 for plotting

semilogx(distances, upsamp, 'ko-', 'LineWidth', 1.5, 'MarkerSize', 6);
xlabel('Distance from surface');
ylabel('Recommended upsampling factor');
title(sprintf('Quadrature Method (tol = %.0e)', tol));
yticks(1:7);
yticklabels({'1 (direct)', '2', '3', '4', '5', '6', 'SQ'});
ylim([0.5, 7.5]);
grid on;

sgtitle('Near-Field Error Estimation');

%% Print summary
fprintf('Distance ranges by quadrature method:\n');
for k = 1:6
    mask = (upsamp == k);
    if any(mask)
        fprintf('  Upsampling %d: distance %.4f to %.4f\n', ...
            k, min(distances(mask)), max(distances(mask)));
    end
end
if any(classification.isSQ)
    fprintf('  Special Quad: distance %.4f to %.4f\n', ...
        min(distances(classification.isSQ)), max(distances(classification.isSQ)));
end

fprintf('\nDone!\n');
