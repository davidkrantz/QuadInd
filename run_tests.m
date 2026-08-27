function results = run_tests()
%RUN_TESTS Run the maintained QuadInd matlab.unittest suite.

root = fileparts(mfilename('fullpath'));
addpath(root);
addpath(fullfile(root, 'tests'));
warningState = warning('query', 'quadind:Diagnostics');
warningCleanup = onCleanup(@() warning(warningState.state, 'quadind:Diagnostics')); %#ok<NASGU>
warning('off', 'quadind:Diagnostics');
suite = matlab.unittest.TestSuite.fromFolder(fullfile(root, 'tests'));
runner = matlab.unittest.TestRunner.withTextOutput();
results = runner.run(suite);

if any([results.Failed])
    error('quadind:tests:failed', '%d of %d tests failed.', ...
        nnz([results.Failed]), numel(results));
else
    disp('All tests passed!')
    thumbsup()
end
end

function thumbsup()
    disp('⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⢀⣤⣤⣀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀');
    disp('⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠸⣿⣿⣿⣷⡀⠀⠀⠀⠀⠀⠀⠀⠀⠀');
    disp('⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⣿⣿⣿⣿⣧⠀⠀⠀⠀⠀⠀⠀⠀⠀');
    disp('⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⢀⣿⣿⣿⣿⣿⠀⠀⠀⠀⠀⠀⠀⠀⠀');
    disp('⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⢀⣾⣿⣿⣿⣿⡟⠀⠀⠀⠀⠀⠀⠀⠀⠀');
    disp('⠀⠀⠀⠀⠀⠀⠀⠀⠀⣰⣿⣿⣿⣿⣿⣿⣷⣶⣶⣶⣤⣤⣀⡀⠀⠀');
    disp('⠀⣀⣀⣀⣀⡀⠀⢀⣾⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⡆⠀');
    disp('⢰⣿⣿⣿⣿⡇⠀⢸⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣅⠀');
    disp('⢸⣿⣿⣿⣿⡇⠀⢸⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⠄');
    disp('⣿⣿⣿⣿⣿⡇⠀⢸⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⠋⠀');
    disp('⣿⣿⣿⣿⣿⡇⠀⢸⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⠀⠀');
    disp('⢸⣿⣿⣿⣿⡇⠀⢸⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⡟⠋⠀⠀');
    disp('⢸⣿⣿⣿⣿⡇⠀⠀⠙⠿⠿⣿⣿⣿⣿⣿⣿⣿⣿⣿⡿⠿⠃⠀⠀⠀');
    disp('⠀⠛⠛⠛⠛⠃⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀');
end