classdef TestPCGrad < matlab.unittest.TestCase
    methods (Test)
        function returnsFiniteOutputs(testCase)
            addpath(genpath('src'));
            testCase.assumeFalse(isempty(ver('symbolic')), 'Symbolic Math Toolbox is required for pcgrad.');
            syms x1 x2
            f = (x1 - 1)^2 + (x2 + 2)^2;
            [x_star, z_star, x_values, iterations] = pcgrad(f, [0; 0], 2, 'A', 'F', 1e-10, 10);
            testCase.verifySize(x_star, [2 1]);
            testCase.verifyGreaterThanOrEqual(iterations, 1);
            testCase.verifyTrue(all(isfinite(x_values(:))));
            testCase.verifyTrue(isfinite(z_star));
        end

        function appliesDefaultsWhenOptionalArgsOmitted(testCase)
            addpath(genpath('src'));
            testCase.assumeFalse(isempty(ver('symbolic')), 'Symbolic Math Toolbox is required for pcgrad.');
            syms x1 x2
            f = (x1 - 1)^2 + (x2 + 2)^2;
            [x_star, z_star, x_values, iterations] = pcgrad(f, [0; 0], 2);
            testCase.verifySize(x_star, [2 1]);
            testCase.verifyGreaterThanOrEqual(iterations, 1);
            testCase.verifyTrue(all(isfinite(x_values(:))));
            testCase.verifyTrue(isfinite(z_star));
        end
    end
end
