function opti_printUpdateInfo
% Print Download / Update Information for OPTI

fprintf('\n------------------------------------------------\n')
fprintf('OPTI Toolbox Update Instructions:\n\n');

fprintf(' Perform ONE of the following options (whichever you prefer)\n\n');

fprintf('Option 1 - SourceTree /  Git GUI\n');
fprintf('1: Ensure you have the OPTI master branched checked out\n');
fprintf('2: Click "Pull"\n');

fprintf('\nOption 2 - Git Command Line\n');
f = which('opti_Install.m');
if (~isempty(f))
    fprintf('1: $ cd "%s"\n', fileparts(f));
else
    fprintf('1: $ cd optiToolboxInstallLoc     (as per your installation)\n');
end
fprintf('2: $ git checkout master\n');
fprintf('3: $ git pull\n');

fprintf('\nOption 3 - Manual Update (not recommended - use Git!)\n');
fprintf('1: In your favourite browser, navigate to: https://github.com/jonathancurrie/OPTI/releases/latest\n');
fprintf('2: Download the Source Code ("Source code (zip)")\n');
fprintf('3: Delete your existing OPTI installation\n');
fprintf('4: Unzip the new installation to the same place you originally installed OPTI\n');

fprintf('\n Once you have completed ONE of the above methods, in MATLAB run:\n');
fprintf('>> opti_Install.m\n');

fprintf('------------------------------------------------\n')