function results = run_smoke_checks()
%RUN_SMOKE_CHECKS Run a small automated smoke-test subset for the repo.

repo_root = fileparts(fileparts(mfilename('fullpath')));
addpath(genpath(fullfile(repo_root, 'src')));
original_dir = pwd;
cleanup_obj = onCleanup(@() cd(original_dir)); %#ok<NASGU>
cd(fullfile(repo_root, 'tests'));

results = runtests({fullfile(repo_root, 'tests', 'TestEulerMethod.m'), ...
                    fullfile(repo_root, 'tests', 'TestNewtonMethod.m'), ...
                    fullfile(repo_root, 'tests', 'TestWeekTab.m'), ...
                    fullfile(repo_root, 'tests', 'TestSirModel.m')});
disp(table(results))
end
