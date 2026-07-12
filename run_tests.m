% RUN_TESTS Run all QuadEst tests, including legacy comparison tests.
%
% Usage from the project root:
%   run_tests
%

verbose = false;
stopOnFailure = true;
includeLegacy = false; % should not pass legacy tests

init;

testClasses = {
    'quadest.test.TestUtil'
    'quadest.test.TestGeometry'
    'quadest.test.TestGrid'
    'quadest.test.TestKernel'
    'quadest.test.TestErrorEstimation'
};

if includeLegacy
    testClasses{end+1} = 'quadest.test.TestLegacyComparison';
end

results = struct();
results.tests = {};
results.passed = 0;
results.failed = 0;
results.errors = 0;
results.skipped = 0;
results.startTime = datetime('now');

if verbose
    fprintf('\n========================================\n');
    fprintf('  QuadEst Full Test Suite\n');
    fprintf('  %s\n', char(results.startTime));
    fprintf('========================================\n\n');
end

for i = 1:numel(testClasses)
    className = testClasses{i};

    if verbose
        fprintf('Running %s...\n', className);
    end

    try
        testObj = feval(className);
        classResults = testObj.run('verbose', verbose);

        results.tests = [results.tests, classResults.tests];
        results.passed = results.passed + classResults.passed;
        results.failed = results.failed + classResults.failed;
        results.errors = results.errors + classResults.errors;
        results.skipped = results.skipped + classResults.skipped;

        if stopOnFailure && (classResults.failed > 0 || classResults.errors > 0)
            if verbose
                fprintf('\nStopping on failure.\n');
            end
            break;
        end
    catch ME
        results.errors = results.errors + 1;
        results.tests = [results.tests, {struct( ...
            'name', className, ...
            'status', 'ERROR', ...
            'message', ME.message, ...
            'time', 0)}];

        if verbose
            fprintf('  ERROR loading test class: %s\n', ME.message);
        end

        if stopOnFailure
            break;
        end
    end

    if verbose
        fprintf('\n');
    end
end

results.endTime = datetime('now');
results.totalTime = seconds(results.endTime - results.startTime);
results.total = results.passed + results.failed + results.errors + results.skipped;

if results.failed == 0 && results.errors == 0
    statusStr = 'PASSED';
else
    statusStr = 'FAILED';
end

results.summary = sprintf('%s: %d/%d tests passed (%.2f s)', ...
    statusStr, results.passed, results.total, results.totalTime);

fprintf('========================================\n');
fprintf('  %s\n', results.summary);
if results.failed > 0
    fprintf('  Failed: %d\n', results.failed);
end
if results.errors > 0
    fprintf('  Errors: %d\n', results.errors);
end
if results.skipped > 0
    fprintf('  Skipped: %d\n', results.skipped);
end
fprintf('========================================\n\n');

if results.failed > 0 || results.errors > 0
    error('quadest:run_tests:failed', ...
        'Test suite failed: %d failed, %d errors.', ...
        results.failed, results.errors);
end
