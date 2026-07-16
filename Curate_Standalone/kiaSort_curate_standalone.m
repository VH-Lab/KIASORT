function fig = kiaSort_curate_standalone(outputFolder, cfg)
% KIASORT_CURATE_STANDALONE - open KIASORT curation in its own window
%
%   FIG = KIASORT_CURATE_STANDALONE(OUTPUTFOLDER)
%   FIG = KIASORT_CURATE_STANDALONE(OUTPUTFOLDER, CFG)
%
%   Opens KIASORT's interactive curation interface (kiaSort_curate_results) in a new
%   standalone figure for a completed sort located in OUTPUTFOLDER (the folder that
%   contains RES_Sorted and Sorted_Samples). This is a convenience entry point for
%   launching curation outside the main KIASORT GUI - for example from another tool -
%   without having to build the figure, panel and config yourself.
%
%   CFG is an optional KIASORT config struct. Any fields it provides override the
%   defaults from kiaSort_main_configs / kiaSort_extended_configs /
%   kiaSort_hidden_configs. At minimum you should pass the samplingFrequency (and
%   numChannels) used for the sort, because the curation UI needs them to render the
%   time axis and waveforms correctly:
%
%       kiaSort_curate_standalone(outputFolder, struct('samplingFrequency',30000, ...
%           'numChannels',384));
%
%   OUTPUTFOLDER is set on the returned config automatically, so the curation UI
%   reads its results from there.
%
%   Returns FIG, the uifigure hosting the curation UI.
%
%   See also: kiaSort_curate_results, run_kiasort_nogui, kiaSort_main_configs

    if nargin < 1 || isempty(outputFolder)
        error('kiaSort_curate_standalone:noOutputFolder', ...
            'An outputFolder (containing RES_Sorted) is required.');
    end
    if ~isfolder(fullfile(outputFolder, 'RES_Sorted'))
        error('kiaSort_curate_standalone:noResults', ...
            'No RES_Sorted folder found in %s. Run KIASORT first.', outputFolder);
    end

    % Start from the standard KIASORT defaults, then overlay any provided cfg fields.
    baseCfg = kiaSort_main_configs();
    baseCfg = kiaSort_extended_configs(baseCfg);
    baseCfg = kiaSort_hidden_configs(baseCfg);
    if nargin >= 2 && isstruct(cfg)
        f = fieldnames(cfg);
        for i = 1:numel(f)
            baseCfg.(f{i}) = cfg.(f{i});
        end
    end
    cfg = baseCfg;

    cfg.outputFolder = outputFolder;
    if ~isfield(cfg, 'inputFolder') || isempty(cfg.inputFolder)
        cfg.inputFolder = outputFolder;
    end
    if ~isfield(cfg, 'altResFolder')
        cfg.altResFolder = '';
    end

    figColor = [0.94 0.94 0.94];
    fig = uifigure('Name', 'KIASORT Curation', ...
        'Position', [100 100 1150 720], 'Color', figColor);

    % kiaSort_curate_results places its widgets directly into this grid using
    % fixed Layout.Row/Column (leftPanel at cols [1 5] rows [2 3], Save button at
    % row 1 col 4, the Visualization/Plots panels at col 6). It does NOT set the
    % grid's dimensions itself, so the parent must be a [3 6] grid with the same
    % proportions as the main GUI's tab3 grid - otherwise the layout collapses
    % (e.g. an enormous Save button filling the top).
    panel = uigridlayout(fig, [3 6], ...
        'RowHeight',   {'fit','fit','1x'}, ...
        'ColumnWidth', {'fit','fit','fit','fit','1x','4x'}, ...
        'Padding',     10, ...
        'RowSpacing',  10);

    kiaSort_curate_results(cfg, panel, figColor, fig);
end
