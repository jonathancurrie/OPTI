%% Steps for a new Release (my use only..)

% Local
% 1) Update optiver version, release notes and date
% 2) Rebuild MEX files with the latest optiver version
% 3) Run opti_Dist_Test (2018b + latest MATLAB) - will automatically package the MEX files and update contents

% Git
% 1) Commit the latest changes with the new version
% 2) Tag the release using the format "OPTI_Toolbox_vX.XX_Released"
% 3) Push the changes to Github
% 4) Under OPTI/Releases, draft a new release. Select the tag, title "OPTI Toolbox vX.XX" and attach the zipped mex files.

% Website
% 1) Update OPTI Wiki page with new features
% 2) Update homepage with new version