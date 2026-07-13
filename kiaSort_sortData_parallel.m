function kiaSort_sortData_parallel(inputPath, outputPath, cfg, varargin)
% KIASORT_SORTDATA_PARALLEL - run stage-3 sorting (kiaSort_main_sortData) in parallel over chunks
%
%   kiaSort_sortData_parallel(inputPath, outputPath, cfg, ...)
%
% Chunk-level parallel wrapper around kiaSort_main_sortData. The recording's
% sorting chunks are partitioned into contiguous ranges, each range is sorted
% on its own parallel worker (running the UNCHANGED per-chunk logic via
% kiaSort_main_sortData with 'chunk_range'/'resultSubdir'/'skip_posthoc'), and
% the per-worker HDF5 outputs are then merged, in chunk order, into RES_Sorted.
% The post-hoc drift-merge/curate is run once at the end.
%
% Chunks are independent (KIASORT reconciles drift post-hoc, not incrementally),
% so the merged result is expected to match a serial run's spike indices and
% labels exactly - validate with kiaSort_compare_sortings.
%
%   inputPath   - the .bin/.dat data file
%   outputPath  - the KIASORT output folder (must already contain RES_Samples
%                 and Sorted_Samples from stages 1-2)
%   cfg         - the fully-built KIASORT config (as stage 3 receives it)
%
% Name/value options:
%   'numWorkers' (see below) - number of parallel workers / chunk partitions.
%       Default: min(4, physical cores - 1), or cfg.numParallelWorker if it is a
%       positive number. Each worker needs roughly one serial per-chunk memory
%       footprint, so MEMORY - not cores - is the limit: start low (2-4) and
%       raise it while watching RAM, and keep the sorting chunk small.
%   'progressfcn' - f(frac,msg), called on the CLIENT as each part completes.
%   'keepParts' (false) - keep the RES_Sorted_partK folders after merging.
%   'verbose' (true)
%
% If the Parallel Computing Toolbox is unavailable or numWorkers<=1, this falls
% back to a single serial kiaSort_main_sortData call (full range).
%
% See also: KIASORT_MAIN_SORTDATA, KIASORT_COMPARE_SORTINGS

    progressFcn = [];
    keepParts   = false;
    verbose     = true;
    numWorkers  = [];
    for i = 1:2:numel(varargin)
        switch lower(varargin{i})
            case 'progressfcn', progressFcn = varargin{i+1};
            case 'numworkers',  numWorkers  = varargin{i+1};
            case 'keepparts',   keepParts   = logical(varargin{i+1});
            case 'verbose',     verbose     = logical(varargin{i+1});
            otherwise, warning('kiaSort_sortData_parallel: unknown option "%s".', varargin{i});
        end
    end

    % Resolve worker count (memory-conservative default).
    if isempty(numWorkers)
        if isfield(cfg,'numParallelWorker') && isnumeric(cfg.numParallelWorker) ...
                && isscalar(cfg.numParallelWorker) && cfg.numParallelWorker >= 1
            numWorkers = round(cfg.numParallelWorker);
        else
            numWorkers = max(1, min(4, feature('numcores') - 1));
        end
    end

    % Number of chunks - mirror kiaSort_main_sortData exactly.
    m = map_input_file(inputPath, cfg);
    num_samples   = size(m.Data.data, 2);
    clear m;
    num_chunk_pts = cfg.sortingChunkDuration * cfg.samplingFrequency;
    num_chunks    = ceil(num_samples / num_chunk_pts);

    havePCT = ~isempty(ver('parallel')) && license('test','Distrib_Computing_Toolbox');

    % Fall back to plain serial stage 3 when we can't or shouldn't parallelise.
    if ~havePCT || numWorkers <= 1 || num_chunks <= 1
        if verbose
            fprintf('kiaSort_sortData_parallel: running serially (workers=%d, chunks=%d, PCT=%d).\n', ...
                numWorkers, num_chunks, havePCT);
        end
        extra = {};
        if ~isempty(progressFcn), extra = {'progressfcn', progressFcn}; end
        kiaSort_main_sortData(inputPath, outputPath, cfg, extra{:});
        return;
    end

    % Partition 1:num_chunks into contiguous ranges (chunk order preserved).
    numParts = min(numWorkers, num_chunks);
    bounds   = round(linspace(0, num_chunks, numParts + 1));
    ranges   = {};
    partSubdirs = {};
    for k = 1:numParts
        c0 = bounds(k) + 1;
        c1 = bounds(k + 1);
        if c1 >= c0
            ranges{end+1}      = [c0 c1];                              %#ok<AGROW>
            partSubdirs{end+1} = sprintf('RES_Sorted_part%d', numel(ranges)); %#ok<AGROW>
        end
    end
    nP = numel(ranges);

    if verbose
        fprintf('kiaSort_sortData_parallel: %d chunks across %d workers (%d parts).\n', ...
            num_chunks, numWorkers, nP);
    end

    % Clean any stale part folders from a previous run.
    for k = 1:nP
        pdir = fullfile(outputPath, partSubdirs{k});
        if exist(pdir, 'dir'), rmdir(pdir, 's'); end
    end

    pool = gcp('nocreate');
    if isempty(pool)
        pool = parpool(numWorkers);
    end

    % Launch one future per part. i_runRange pins each worker to a single
    % computational thread so N workers don't oversubscribe the CPU with N x
    % BLAS threads each.
    futures = parallel.FevalFuture.empty(0, nP);
    for k = 1:nP
        futures(k) = parfeval(pool, @i_runRange, 0, ...
            inputPath, outputPath, cfg, ranges{k}, partSubdirs{k}); %#ok<AGROW>
    end

    % Wait, reporting progress as each part finishes (fetchNext rethrows worker errors).
    for done = 1:nP
        fetchNext(futures);
        if ~isempty(progressFcn)
            progressFcn(done / nP, sprintf('Sorted part %d/%d', done, nP));
        elseif verbose
            fprintf('  part %d/%d complete.\n', done, nP);
        end
    end

    % Merge part outputs (in chunk order) into RES_Sorted.
    i_mergeParts(outputPath, partSubdirs, verbose);

    % Remove part folders unless asked to keep them.
    if ~keepParts
        for k = 1:nP
            pdir = fullfile(outputPath, partSubdirs{k});
            if exist(pdir, 'dir'), rmdir(pdir, 's'); end
        end
    end

    % Post-hoc drift merge + curate, once, on the merged result (mirror sortData).
    if isfield(cfg,'postHocProcessing') && cfg.postHocProcessing && ...
            ~(isfield(cfg,'sort_only') && cfg.sort_only)
        kiaSort_drift_merge_posthoc_iterative(outputPath, ...
            'overwrite', true, 'verbose', false, 'mainArgs', {'debugFigs', false});
        try
            kiaSort_post_sort_curate(outputPath, ...
                'ccg_cleaning', true, 'merging', true, ...
                'xcorrThreshold', 0.9, 'verbose', false);
        catch
        end
    end

    if verbose
        fprintf('kiaSort_sortData_parallel: done. Merged results in %s\n', ...
            fullfile(outputPath, 'RES_Sorted'));
    end
end

% ------------------------------------------------------------------------
function i_runRange(inputPath, outputPath, cfg, range, subdir)
% Worker task: sort one contiguous chunk range into its own RES_Sorted subfolder,
% skipping the post-hoc step (the client runs that once after merging).
    try, maxNumCompThreads(1); catch, end
    kiaSort_main_sortData(inputPath, outputPath, cfg, ...
        'chunk_range', range, 'resultSubdir', subdir, 'skip_posthoc', true);
end

% ------------------------------------------------------------------------
function i_mergeParts(outputPath, partSubdirs, verbose)
% Concatenate each per-part HDF5 field, in part (=chunk) order, into RES_Sorted.
% spike indices are already absolute, so this is a straight ordered concat and
% reproduces the byte layout a serial run would have written.
    finalDir = fullfile(outputPath, 'RES_Sorted');
    if ~exist(finalDir, 'dir')
        mkdir(finalDir);
    else
        old = dir(fullfile(finalDir, '*.h5'));
        for i = 1:numel(old), delete(fullfile(finalDir, old(i).name)); end
    end

    % Union of field files present across all parts (a part with no spikes for
    % a field simply won't have produced that file).
    fieldFiles = {};
    for k = 1:numel(partSubdirs)
        d = dir(fullfile(outputPath, partSubdirs{k}, '*.h5'));
        for i = 1:numel(d), fieldFiles{end+1} = d(i).name; end %#ok<AGROW>
    end
    fieldFiles = unique(fieldFiles);

    for i = 1:numel(fieldFiles)
        fname = fieldFiles{i};
        [~, base] = fileparts(fname);
        dset = ['/' base];
        data = [];
        for k = 1:numel(partSubdirs)
            pf = fullfile(outputPath, partSubdirs{k}, fname);
            if exist(pf, 'file')
                data = cat(1, data, h5read(pf, dset));   % concat along rows = chunk order
            end
        end
        if isempty(data), continue; end
        outFile = fullfile(finalDir, fname);
        h5create(outFile, dset, size(data));
        h5write(outFile, dset, data);
    end

    if verbose
        fprintf('  merged %d field(s) from %d parts into RES_Sorted.\n', ...
            numel(fieldFiles), numel(partSubdirs));
    end
end
