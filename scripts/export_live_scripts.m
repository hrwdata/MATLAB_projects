function export_live_scripts()
%EXPORT_LIVE_SCRIPTS Export selected live scripts to HTML.
% This helper relies on an internal MATLAB API and may require updates
% across MATLAB releases.

if exist('matlab.internal.liveeditor.openAndConvert', 'file') ~= 2
    error('export_live_scripts:UnsupportedMATLAB', ...
          ['Live Script export requires matlab.internal.liveeditor.openAndConvert. ', ...
           'Use a MATLAB release that provides it or export manually from the editor.']);
end

targets = {
    'examples/tutorial/matlab_tutorial.mlx'
    'examples/tutorial/matlab_principles.mlx'
    'examples/bifurcation/bifurcation.mlx'
    'examples/epidemiology/SIR_fit_AL.mlx'
};

for idx = 1:numel(targets)
    input_file = targets{idx};
    [folder, name] = fileparts(input_file);
    output_file = fullfile(folder, [name '.html']);
    try
        matlab.internal.liveeditor.openAndConvert(input_file, output_file);
    catch ME
        error('export_live_scripts:ConversionFailed', ...
              'Failed to export %s: %s', input_file, ME.message);
    end
end
end
