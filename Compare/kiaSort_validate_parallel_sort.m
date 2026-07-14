function result = kiaSort_validate_parallel_sort(inputPath, outputPath, numChannels, fs, varargin)
% KIASORT_VALIDATE_PARALLEL_SORT - verify chunk-parallel stage 3 matches serial
%
%   result = kiaSort_validate_parallel_sort(inputPath, outputPath, numChannels, fs, ...)
%
% Runs stage 3 (sortData) both SERIALLY and in PARALLEL on the SAME frozen
% stage-1/stage-2 output already sitting in outputPath (RES_Samples +
% Sorted_Samples from a prior run), then compares the two raw results. Stage 3
% is deterministic given the frozen templates, so an exact match confirms the
% parallelization (kiaSort_sortData_parallel) is correct. Post-hoc is skipped in
% both runs so the comparison is of the raw stage-3 output.
%
% This is the clean way to validate the parallelization: it isolates stage 3
% from the (stochastic) stage-1 chunk draw and stage-2 UMAP by reusing their
% frozen output for both runs.
%
%   inputPath   - the .bin/.dat used for the sort
%   outputPath  - folder with RES_Samples/channel_info.mat and
%                 Sorted_Samples/sorted_samples.mat from a prior run
%   numChannels - channel count
%   fs          - sampling rate (Hz)
%
% Options (name/value):
%   'numWorkers'    (4)          parallel workers to use
%   'dataType'      ('int16')    binary data type
%   'cfg_overrides' (struct())   extra cfg fields applied after defaults; pass
%                                the same options your real run used (e.g.
%                                denoising / extremeNoise / sortingChunkDuration)
%                                so the sort is representative. Not required for
%                                correctness of the comparison (both runs use
%                                the identical cfg), only for realism.
%   'verbose'       (true)
%
% Returns the kiaSort_compare_sortings result struct and prints PASS/FAIL.
% Writes RES_Sorted_serial (serial) and RES_Sorted (parallel) under outputPath.
%
% See also: KIASORT_SORTDATA_PARALLEL, KIASORT_COMPARE_SORTINGS, KIASORT_MAIN_SORTDATA

    p.numWorkers    = 4;
    p.dataType      = 'int16';
    p.cfg_overrides = struct();
    p.verbose       = true;
    for i = 1:2:numel(varargin)
        key = varargin{i};
        if ~isfield(p, key)
            error('kiaSort_validate_parallel_sort: unknown option "%s".', key);
        end
        p.(key) = varargin{i+1};
    end

    if ~exist(fullfile(outputPath,'RES_Samples','channel_info.mat'),'file') || ...
       ~exist(fullfile(outputPath,'Sorted_Samples','sorted_samples.mat'),'file')
        error(['outputPath must already contain RES_Samples/channel_info.mat and ' ...
               'Sorted_Samples/sorted_samples.mat (from a prior stages 1-2 run).']);
    end

    % Build cfg the way the pipeline does; the SAME cfg is used for both runs,
    % so the serial-vs-parallel comparison is valid regardless of how closely it
    % matches the original run.
    cfg = kiaSort_main_configs();
    cfg = kiaSort_extended_configs(cfg);
    cfg = kiaSort_hidden_configs(cfg);
    cfg.numChannels       = numChannels;
    cfg.samplingFrequency = fs;
    cfg.dataType          = p.dataType;
    fn = fieldnames(p.cfg_overrides);
    for i = 1:numel(fn), cfg.(fn{i}) = p.cfg_overrides.(fn{i}); end

    % Path fields the stage functions expect (run_kiasort_nogui normally sets
    % these; map_input_file uses cfg.outputFolder for its log).
    cfg.fullFilePath = inputPath;
    cfg.inputFolder  = fileparts(inputPath);
    cfg.outputFolder = outputPath;

    % Clean any prior validation outputs.
    for sd = {'RES_Sorted', 'RES_Sorted_serial'}
        d = fullfile(outputPath, sd{1});
        if exist(d, 'dir'), rmdir(d, 's'); end
    end

    if p.verbose, fprintf('[validate] SERIAL stage 3 -> RES_Sorted_serial ...\n'); end
    tSerial = tic;
    kiaSort_main_sortData(inputPath, outputPath, cfg, ...
        'resultSubdir', 'RES_Sorted_serial', 'skip_posthoc', true);
    tSerial = toc(tSerial);

    if p.verbose, fprintf('[validate] PARALLEL stage 3 (%d workers) -> RES_Sorted ...\n', p.numWorkers); end
    tParallel = tic;
    kiaSort_sortData_parallel(inputPath, outputPath, cfg, ...
        'numWorkers', p.numWorkers, 'skipPostHoc', true, 'verbose', p.verbose);
    tParallel = toc(tParallel);

    fprintf('[validate] timing: serial %.1f s | parallel %.1f s | speedup %.2fx (%d workers)\n', ...
        tSerial, tParallel, tSerial / max(tParallel, eps), p.numWorkers);

    result = kiaSort_compare_sortings(outputPath, outputPath, fs, ...
        'subdirA', 'RES_Sorted_serial', 'subdirB', 'RES_Sorted', ...
        'nameA', 'serial', 'nameB', 'parallel', 'verbose', p.verbose);
    result.timing = struct('serial_s', tSerial, 'parallel_s', tParallel, ...
        'speedup', tSerial / max(tParallel, eps), 'numWorkers', p.numWorkers);

    if result.exact.same_spike_times && isequal(result.exact.same_labels, true)
        fprintf('\n[validate] PASS: parallel stage 3 is identical to serial.\n');
    else
        fprintf(['\n[validate] FAIL: parallel and serial differ - see the report above. ' ...
                 '(same_spike_times=%d)\n'], result.exact.same_spike_times);
    end
end
