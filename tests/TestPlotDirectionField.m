classdef TestPlotDirectionField < matlab.unittest.TestCase
    methods (Test)
        function usesDefaultStepWhenOmitted(testCase)
            addpath(genpath('src'));
            fig = figure('Visible', 'off');
            testCase.onCleanup(@() close(fig));
            handle = plotd(@(t,y) t + y, 0, 1, -1, 1);
            testCase.verifyClass(handle, 'matlab.graphics.chart.primitive.Quiver');
        end
    end
end
