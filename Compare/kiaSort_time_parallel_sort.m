function t = kiaSort_time_parallel_sort(inputPath, outputPath, numChannels, fs, varargin)
% KIASORT_TIME_PARALLEL_SORT - time the parallel stage-3 orchestrator (tuning only)
%
%   t = kiaSort_time_parallel_sort(inputPath, outputPath, numChannels, fs, ...)
%
% Runs kiaSort_sortData_parallel once (post-hoc skipped) and returns/prints the
% elapsed seconds. No serial baseline, no comparison - use this to sweep worker
% count and chunk size quickly once correctness is already validated (via
% kiaSort_validate_parallel_sort). outputPath must contain RES_Samples and
% Sorted_Samples (stages 1-2).
%
% Options:
%   'numWorkers'    (4)
%   'dataType'      ('int16')
%   'cfg_overrides' (struct())  - e.g. struct('denoising',true,'sortingChunkDuration',120)
%
% Example sweep at fixed chunk size:
%   for w = [4 6 8]
%       kiaSort_time_parallel_sort(inputPath, outputPath, 384, 30000, ...
%           'numWorkers', w, 'cfg_overrides', struct('denoising',true,'sortingChunkDuration',120));
%   end
%
% See also: KIASORT_VALIDATE_PARALLEL_SORT, KIASORT_SORTDATA_PARALLEL

    p.numWorkers    = 4;
    p.dataType      = 'int16';
    p.cfg_overrides = struct();
    for i = 1:2:numel(varargin)
        key = varargin{i};
        if ~isfield(p, key)
            error('kiaSort_time_parallel_sort: unknown option "%s".', key);
        end
        p.(key) = varargin{i+1};
    end

    cfg = kiaSort_main_configs();
    cfg = kiaSort_extended_configs(cfg);
    cfg = kiaSort_hidden_configs(cfg);
    cfg.numChannels       = numChannels;
    cfg.samplingFrequency = fs;
    cfg.dataType          = p.dataType;
    fn = fieldnames(p.cfg_overrides);
    for i = 1:numel(fn), cfg.(fn{i}) = p.cfg_overrides.(fn{i}); end
    cfg.fullFilePath = inputPath;
    cfg.inputFolder  = fileparts(inputPath);
    cfg.outputFolder = outputPath;

    tt = tic;
    kiaSort_sortData_parallel(inputPath, outputPath, cfg, ...
        'numWorkers', p.numWorkers, 'skipPostHoc', true, 'verbose', false);
    t = toc(tt);

    fprintf('[time] %2d workers, %g s chunks -> %.1f s\n', ...
        p.numWorkers, cfg.sortingChunkDuration, t);
end
