% demo_validation.m - Practical indicator validation against converged references.

scriptDir = fileparts(mfilename('fullpath'));
addpath(fileparts(scriptDir));

geometry = quadind.geometry.Spheroid('a', 0.05, 'c', 0.1);
kernel = quadind.kernel.StokesStresslet();
nth = 24;
nph = 32;
tol = 1e-6;
factors = 1:3;
config = quadind.util.Config('ntab', 35, 'nztab', 35);

fprintf('Building QuadInd tables...\n');
evaluator = quadind.IndicatorEvaluator(geometry, kernel, ...
    'nth', nth, 'nph', nph, 'upsampFactors', factors, 'config', config);
baseGrid = evaluator.getGrid();
density = analyticDensity(baseGrid);

radial = linspace(0.06, 0.14, 28);
axial = linspace(-0.13, 0.13, 28);
[R, Z] = meshgrid(radial, axial);
targets = [R(:), zeros(numel(R),1), Z(:)];
targets = targets(geometry.isExterior(targets), :);

[indicator, classification] = evaluator.evaluate(targets, density, 'tol', tol);

% Evaluate the manufactured density directly on each refined grid.
referenceFactors = [4, 7];
reference = cell(size(referenceFactors));
for k = 1:numel(referenceFactors)
    fineGrid = baseGrid.upsample(referenceFactors(k));
    reference{k} = kernel.evaluateOnGrid(targets, fineGrid, analyticDensity(fineGrid));
end
referenceUncertainty = vecnorm(reference{2} - reference{1}, 2, 2);
referenceValue = reference{2};

directValue = kernel.evaluateOnGrid(targets, baseGrid, density);
measuredError = vecnorm(directValue - referenceValue, 2, 2);
reliable = referenceUncertainty <= max(1e-1 * measuredError, 100*eps);
if ~all(reliable)
    warning('quadind:demo_validation:referenceConvergence', ...
        '%d of %d targets did not meet the reference-convergence criterion.', ...
        nnz(~reliable), numel(reliable));
end

falseSafe = reliable & indicator < tol & measuredError >= tol;
fprintf('Targets: %d; converged references: %d; false-safe classifications: %d\n', ...
    numel(reliable), nnz(reliable), nnz(falseSafe));
if any(falseSafe)
    falseSafeIndex = find(falseSafe);
    fprintf('\nFalse-safe target details (indicator < tol <= measured error):\n');
    fprintf('  Index          x          y          z       Indicator    Measured error\n');
    for k = 1:numel(falseSafeIndex)
        i = falseSafeIndex(k);
        fprintf('  %5d  %9.3e  %9.3e  %9.3e  %12.5e  %14.5e\n', ...
            i, targets(i,1), targets(i,2), targets(i,3), ...
            indicator(i), measuredError(i));
    end
end
fprintf('Targets requiring an external special quadrature method: %d\n', ...
    nnz(classification.requiresSpecialQuadrature));

figure('Name', 'QuadInd validation');
loglog(measuredError(reliable), indicator(reliable), '.');
hold on;
limits = [1e-14, 1e1];
loglog(limits, limits, 'k--');
xlim(limits); ylim(limits); grid on; axis square;
xlabel('Measured Euclidean-norm quadrature error');
ylabel('QuadInd error indicator');
title('Direct quadrature: measured error versus indicator');

function density = analyticDensity(surfaceGrid)
[theta, phi] = surfaceGrid.meshgrid();
density = [2.1 + sin(4*theta(:)) + sin(4*phi(:)), ...
    2 + sin(3*theta(:)).*cos(3*phi(:)), ...
    1.03 + cos(2*theta(:)).*cos(2*phi(:))];
end
