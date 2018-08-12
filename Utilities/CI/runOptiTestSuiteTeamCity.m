function runOptiTestSuiteTeamCity(doExit)
% runTestSuiteTeamCity
%
%   Entry point for running OPTI unit tests in CI

if (nargin < 1 || isempty(doExit)), doExit = true; end

% Clear MATLAB to default state
try
    restoredefaultpath;
%     matlabrc;
catch ME
    printf(2, 'Error initializing MATLAB! %s', ME.message);
    ciexit(3, doExit);
end

% Install the toolbox
try
    curDir = cd;
    cd('../../');
    opti_Install(false, false, false, false);
    cd(curDir);
catch ME
    fprintf(2, 'Error installing OPTI Toolbox! %s', ME.message);
    ciexit(2, doExit);
end
% 
% % Run the unit tests
% try
%     import matlab.unittest.TestRunner
%     import matlab.unittest.TestSuite
%     import matlab.unittest.plugins.TAPPlugin
%     import matlab.unittest.plugins.ToFile
% 
%     suite = TestSuite.fromClass(?barvec_tests);
%     runner = TestRunner.withTextOutput;
%     tapFile = 'BARONTestOutput.tap';
%     plugin = TAPPlugin.producingVersion13(ToFile(tapFile));
%     runner.addPlugin(plugin)
%     runner.run(suite);
% catch ME
%     fprintf(2, 'Error running unit tests! %s', ME.message);
%     ciexit(1, doExit);
% end
% 
% % Display unit test output
% try
%     display(fileread(tapFile));
% catch ME
%     fprintf(2, 'Error reading unit test results! %s', ME.message);
%     ciexit(1, doExit);
% end

% Success
ciexit(0, doExit);



function ciexit(code, doExit)
if (doExit)
    exit(code);
elseif (code == 0)
    fprintf('Completed successfully\n');
else
    error('Failed with code %d', code);
end
    