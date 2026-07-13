function run_kiasort_nogui(dataFilePath, outputFolder, channelMapFile, cfg_overrides, varargin)
% Run kiaSort end-to-end without the GUI.
%
%   dataFilePath   - .dat / .bin file
%   outputFolder   - directory for results
%   channelMapFile - optional .mat with channel map fields (empty to skip)
%   cfg_overrides  - optional struct of cfg fields to override after defaults
%
%   Optional name-value pairs:
%     'progressfcn' - a function handle f(pct, msg) called during sorting. PCT is a
%                     fraction in [0,1] over the WHOLE run (the three stages are
%                     mapped onto consecutive slices of that range) and MSG is a
%                     short status string. Use it to drive a progress bar.
%     'verbose'     - logical (default true). When false, the summary/"Done" prints
%                     are suppressed (KIASORT still writes its log file).

if nargin < 2
    error('dataFilePath and outputFolder are required');
end
if nargin < 3, channelMapFile = []; end
if nargin < 4 || isempty(cfg_overrides), cfg_overrides = struct(); end

% Parse optional name-value pairs.
progressFcn = [];
verbose = true;
for i = 1:2:numel(varargin)
    key = lower(varargin{i});
    val = varargin{i+1};
    switch key
        case 'progressfcn'
            progressFcn = val;
        case 'verbose'
            verbose = logical(val);
        otherwise
            error('run_kiasort_nogui: unknown option "%s".', varargin{i});
    end
end

if ~exist(dataFilePath, 'file')
    error('Input data file not found: %s', dataFilePath);
end
if ~exist(outputFolder, 'dir'), mkdir(outputFolder); end

cfg = kiaSort_main_configs();
cfg = kiaSort_extended_configs(cfg);
cfg = kiaSort_hidden_configs(cfg);

ovrFields = fieldnames(cfg_overrides);
for i = 1:numel(ovrFields)
    cfg.(ovrFields{i}) = cfg_overrides.(ovrFields{i});
end

[cfg.inputFolder, ~, ~] = fileparts(dataFilePath);
cfg.fullFilePath = dataFilePath;
cfg.outputFolder = outputFolder;

hp = sorting_hyperparameters_in();

channel_mapping = [];
channel_inclusion = [];
channel_locations = [];
if ~isempty(channelMapFile) && exist(channelMapFile, 'file')
    cfg.channel_info = channelMapFile;
    [channel_mapping, channel_locations, channel_inclusion] = ...
        load_channel_map(channelMapFile, cfg);
elseif isfield(cfg, 'numChannels')
    channel_mapping   = 1:cfg.numChannels;
    channel_inclusion = true(cfg.numChannels, 1);
    channel_locations = [];
end

cfg.num_channel_extract = derive_num_channel_extract(channel_locations, ...
    cfg.waveform_radius, cfg.num_channel_extract);

if verbose
    fprintf('Input: %s\nOutput: %s\nChannels: %d, fs: %d Hz, dtype: %s, BP: [%d %d]\n', ...
        dataFilePath, outputFolder, cfg.numChannels, cfg.samplingFrequency, ...
        cfg.dataType, cfg.bandpass(1), cfg.bandpass(2));
end

% Wrap the caller's progress callback so each stage reports into its own slice of an
% overall [0,1] range (extract, then optionally sort_samples, then sortData).
doSort  = ~cfg.sort_only;
nStages = 2 + double(doSort);
extractExtra = {};
sortExtra    = {};
dataExtra    = {};
if ~isempty(progressFcn)
    if i_hasVararg('kiaSort_main_extract_sample_data')
        extractExtra = {'progressfcn', @(p, m) progressFcn((0 + p) / nStages, m)};
    end
    if doSort && i_hasVararg('kiaSort_main_sort_samples')
        sortExtra = {'progressfcn', @(p, m) progressFcn((1 + p) / nStages, m)};
    end
    if i_hasVararg('kiaSort_main_sortData')
        base = 2; if ~doSort, base = 1; end
        dataExtra = {'progressfcn', @(p, m) progressFcn((base + p) / nStages, m)};
    end
end

kiaSort_main_extract_sample_data(cfg.fullFilePath, cfg.outputFolder, cfg, ...
    'channel_mapping',   channel_mapping, ...
    'channel_inclusion', channel_inclusion, ...
    'channel_locations', channel_locations, ...
    extractExtra{:});

if doSort
    kiaSort_main_sort_samples(cfg.outputFolder, cfg, hp, sortExtra{:});
end

kiaSort_main_sortData(cfg.fullFilePath, cfg.outputFolder, cfg, dataExtra{:});

if verbose
    fprintf('Done. Results in %s\n', outputFolder);
end

end

function tf = i_hasVararg(fnname)
% True if the function FNNAME accepts varargin (can take extra name-value pairs).
tf = false;
try
    tf = nargin(fnname) < 0;
catch
end
end
