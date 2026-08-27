% demo_nearfield.m - Near-field error-indicator behavior
%
% This script demonstrates error-indicator evaluation for targets at varying
% distances from the surface, focusing on the near-singular regime.
%
% Run from the QuadInd directory:
%   >> run('examples/demo_nearfield.m')

%% Setup path
scriptDir = fileparts(mfilename('fullpath'));
projectDir = fileparts(scriptDir);
addpath(projectDir);

%% Setup
fprintf('=== NEAR-FIELD DEMONSTRATION ===\n\n');

% Geometry: spheroid
geom = quadind.geometry.Spheroid('a', 0.05, 'c', 0.1);
kernel = quadind.kernel.StokesStresslet();

fprintf('Geometry: ');
disp(geom);

% Precompute
evaluator = quadind.IndicatorEvaluator(geom, kernel, ...
    'nth', 40, 'nph', 60, ...
    'upsampFactors', 1:6, ...
    'interpolateRoots', true);

grid_est = evaluator.getGrid();
density = ones(grid_est.numPoints(), 3);  % Stresslet identity

%% Create targets at varying distances
% Start from surface and move outward in x-direction
a = 0.05;  % Semi-axis
distances = logspace(-4, 0, 50);  % 0.0001 to 1
targets = [(a + distances)', zeros(length(distances), 2)];

fprintf('Testing %d targets at distances from %.4f to %.2f\n\n', ...
    length(distances), min(distances), max(distances));

%% Evaluate indicators for each upsampling factor
tol = 1e-6;
[indicators, classification] = evaluator.evaluate(targets, density, 'tol', tol);

%% Plot results
figure('Position', [100, 100, 1000, 400]);

% Plot 1: Error indicator vs distance
subplot(1, 2, 1);
loglog(distances, indicators, 'b-', 'LineWidth', 2);
hold on;
yline(tol, 'r--', 'LineWidth', 1.5, 'Label', sprintf('tol = %.0e', tol));
xlabel('Distance from surface');
ylabel('Error indicator');
title('Error Indicator vs Distance');
grid on;
legend('Indicator', 'Tolerance', 'Location', 'northeast');

% Plot 2: Recommended upsampling vs distance
subplot(1, 2, 2);
upsamp = classification.upsamplingFactor;
sqPlotValue = max(evaluator.upsampFactors) + 1;
upsamp(classification.requiresSpecialQuadrature) = sqPlotValue;

semilogx(distances, upsamp, 'ko-', 'LineWidth', 1.5, 'MarkerSize', 6);
xlabel('Distance from surface');
ylabel('Recommended upsampling factor');
title(sprintf('Quadrature Method (tol = %.0e)', tol));
yticks(1:sqPlotValue);
yticklabels([compose('%d', evaluator.upsampFactors), "SQ"]);
ylim([0.5, sqPlotValue + 0.5]);
grid on;

sgtitle('Near-Field Error Indicators');

%% Print summary
fprintf('Distance ranges by quadrature method:\n');
for k = evaluator.upsampFactors
    mask = (upsamp == k);
    if any(mask)
        fprintf('  Upsampling %d: distance %.4f to %.4f\n', ...
            k, min(distances(mask)), max(distances(mask)));
    end
end
if any(classification.requiresSpecialQuadrature)
    fprintf('  Special Quad: distance %.4f to %.4f\n', ...
        min(distances(classification.requiresSpecialQuadrature)), max(distances(classification.requiresSpecialQuadrature)));
end

fprintf('\nDone!\n');
