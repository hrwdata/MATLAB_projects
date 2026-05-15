classdef TestEulerMethod < matlab.unittest.TestCase
    methods (Test)
        function computesExpectedGrid(testCase)
            repo_root = fileparts(fileparts(mfilename('fullpath')));
            addpath(genpath(fullfile(repo_root, 'src')));
            original_dir = pwd;
            cleanup_obj = onCleanup(@() cd(original_dir)); %#ok<NASGU>
            cd(fullfile(repo_root, 'tests'));
            original_format = evalc('format');
            [x, y, f, table_out] = eulers(@(x,y) x + y, 0, 0.2, 1, 0.1);
            testCase.verifyEqual(x, [0 0.1 0.2], 'AbsTol', 1e-12);
            testCase.verifyEqual(y(1), 1, 'AbsTol', 1e-12);
            testCase.verifyEqual(numel(f), numel(x));
            testCase.verifySize(table_out, [3 3]);
            testCase.verifyEqual(evalc('format'), original_format);
        end
    end
end
