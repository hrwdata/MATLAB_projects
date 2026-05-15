classdef TestLineSearch < matlab.unittest.TestCase
    methods (Test)
        function returnsPositiveStepSizes(testCase)
            addpath(genpath('src'));
            testCase.assumeFalse(isempty(ver('symbolic')), 'Symbolic Math Toolbox is required for line-search tests.');
            syms x1 x2
            f = x1^2 + x2^2;
            x0 = [1; 1];
            d = [-1; -1];
            testCase.verifyGreaterThan(armijo(f, x0, d), 0);
            testCase.verifyGreaterThan(goldstein(f, x0, d), 0);
            testCase.verifyGreaterThan(wolfe(f, x0, d), 0);
        end

        function acceptsExplicitEmptyOptionalArgs(testCase)
            addpath(genpath('src'));
            testCase.assumeFalse(isempty(ver('symbolic')), 'Symbolic Math Toolbox is required for line-search tests.');
            syms x1 x2
            f = x1^2 + x2^2;
            x0 = [1; 1];
            d = [-1; -1];
            testCase.verifyGreaterThan(armijo(f, x0, d, [], []), 0);
            testCase.verifyGreaterThan(goldstein(f, x0, d, [], []), 0);
            testCase.verifyGreaterThan(wolfe(f, x0, d, [], []), 0);
        end
    end
end
