function result = kiaSort_compare_sortings(A, B, fs, varargin)
% KIASORT_COMPARE_SORTINGS - quantitatively compare two spike sortings
%
%   result = kiaSort_compare_sortings(A, B, fs, ...)
%
% Compares two sortings by matching their units via spike-train coincidence
% and reporting per-unit agreement, so you can quantify how similar two runs
% are instead of eyeballing the curation GUI. Works for two independent sorts
% (e.g. validating a code change such as parallelisation) and for GROUND-TRUTH
% comparison (pass the ground truth as A) - it is the same algorithm.
%
% A and B may each be:
%   * a results folder path (char/string) - loaded with kiaSort_load_results
%     (the folder that contains the RES_Sorted subfolder), or
%   * a struct with fields .spike_idx and .unifiedLabels (as returned by
%     kiaSort_load_results), or
%   * a struct with fields .spikeIdx and .labels (handy for ground truth:
%     struct('spikeIdx', gt_samples, 'labels', gt_unit_ids)).
%
% Spike indices are in SAMPLES; FS is the sampling rate (Hz), needed to turn
% the coincidence window into samples.
%
% Name/value options:
%   'window_ms'       (0.4)  Two spikes count as coincident if within this
%                            many milliseconds.
%   'match_threshold' (0.5)  Agreement score at/above which a unit pair is
%                            called a match (SpikeInterface uses 0.5).
%   'ignore_labels'   ([])   Labels to drop from BOTH sortings before
%                            comparing (e.g. 0 for a noise/unsorted label).
%   'nameA'/'nameB'   ('A'/'B')  Names used in the printed report.
%   'verbose'         (true) Print a human-readable report.
%
% The agreement score for a unit pair is the SpikeInterface convention
%       agreement = n_match / (n_A + n_B - n_match)         (= TP/(TP+FN+FP))
% where n_match is the number of coincident spikes (each spike matched at most
% once within a pair), n_A / n_B are the two units' spike counts.
%
% RESULT is a struct with:
%   .unitsA, .unitsB          unit ids in each sorting (after ignore_labels)
%   .countsA, .countsB        spike counts per unit
%   .agreement                nA-by-nB agreement-score matrix
%   .coincidence              nA-by-nB coincident-spike-count matrix
%   .matches                  table of 1-to-1 matched pairs (greedy on
%                             agreement): unitA, unitB, agreement, nMatch,
%                             nOnlyA (misses), nOnlyB (false positives)
%   .unmatchedA, .unmatchedB  units with no match >= threshold, each with the
%                             best partial overlap and a split/merge hint
%   .summary                  overall counts and mean matched agreement
%   .exact                    when the two sortings have identical spike times:
%                             .same_spike_times, .same_labels, .label_agreement
%                             (this is the strong test for a deterministic
%                             change like stage-3 parallelisation - reuse the
%                             SAME stage-1/2 output for both runs and expect
%                             same_labels == true)
%
% Example - validate a deterministic change (expect an exact match):
%   r = kiaSort_compare_sortings(serialDir, parallelDir, 30000);
%   assert(r.exact.same_spike_times && r.exact.same_labels)
%
% Example - ground truth on simulated data:
%   gt = struct('spikeIdx', gtSamples, 'labels', gtUnitIds);
%   r  = kiaSort_compare_sortings(gt, sortedDir, 30000, 'nameA','truth','nameB','kiaSort');
%
% See also: kiaSort_load_results, estimateContamination, refractoryMatrix

    opts.window_ms       = 0.4;
    opts.match_threshold = 0.5;
    opts.ignore_labels   = [];
    opts.nameA           = 'A';
    opts.nameB           = 'B';
    opts.verbose         = true;
    opts = i_parseOpts(opts, varargin);

    if nargin < 3 || isempty(fs)
        error('kiaSort_compare_sortings:fs', 'FS (sampling rate in Hz) is required.');
    end

    [spkA, lblA] = i_getSorting(A);
    [spkB, lblB] = i_getSorting(B);

    if ~isempty(opts.ignore_labels)
        keepA = ~ismember(lblA, opts.ignore_labels);
        keepB = ~ismember(lblB, opts.ignore_labels);
        spkA = spkA(keepA); lblA = lblA(keepA);
        spkB = spkB(keepB); lblB = lblB(keepB);
    end

    win = max(1, round(opts.window_ms * 1e-3 * fs));

    unitsA = unique(lblA(:))';
    unitsB = unique(lblB(:))';
    nA = numel(unitsA);
    nB = numel(unitsB);

    % Presort each unit's spike times once.
    spikesA = cell(nA,1); countsA = zeros(nA,1);
    for i = 1:nA
        spikesA{i} = sort(spkA(lblA == unitsA(i)));
        countsA(i) = numel(spikesA{i});
    end
    spikesB = cell(nB,1); countsB = zeros(nB,1);
    for j = 1:nB
        spikesB{j} = sort(spkB(lblB == unitsB(j)));
        countsB(j) = numel(spikesB{j});
    end

    % Coincidence + agreement matrices.
    coincidence = zeros(nA, nB);
    agreement   = zeros(nA, nB);
    for i = 1:nA
        for j = 1:nB
            m = i_countCoincident(spikesA{i}, spikesB{j}, win);
            coincidence(i,j) = m;
            denom = countsA(i) + countsB(j) - m;
            if denom > 0
                agreement(i,j) = m / denom;
            end
        end
    end

    % Greedy 1-to-1 assignment on the agreement matrix (highest first).
    matchB     = zeros(nA,1);   % index into unitsB (0 = unmatched)
    matchScore = zeros(nA,1);
    if nA > 0 && nB > 0
        [vals, order] = sort(agreement(:), 'descend');
        [ii, jj] = ind2sub([nA nB], order);
        usedB = false(nB,1);
        for k = 1:numel(order)
            if vals(k) <= 0, break; end
            a = ii(k); b = jj(k);
            if matchB(a) == 0 && ~usedB(b)
                matchB(a) = b; matchScore(a) = vals(k); usedB(b) = true;
            end
        end
    end
    isMatched = matchScore >= opts.match_threshold & matchB > 0;

    % Matched-pair table.
    mi = find(isMatched);
    matches = struct('unitA',{},'unitB',{},'agreement',{},'nMatch',{}, ...
        'nOnlyA',{},'nOnlyB',{});
    for t = 1:numel(mi)
        a = mi(t); b = matchB(a); m = coincidence(a,b);
        matches(t).unitA     = unitsA(a);
        matches(t).unitB     = unitsB(b);
        matches(t).agreement = agreement(a,b);
        matches(t).nMatch    = m;
        matches(t).nOnlyA    = countsA(a) - m;   % misses (if A is ground truth)
        matches(t).nOnlyB    = countsB(b) - m;   % false positives
    end

    % Unmatched A units - report best partial overlap + split hint.
    ua = find(~isMatched);
    unmatchedA = i_unmatchedReport(ua, unitsA, countsA, agreement, coincidence, unitsB, 'A');

    % Unmatched B units (never claimed as a match target).
    claimedB = false(nB,1); claimedB(matchB(isMatched)) = true;
    ub = find(~claimedB);
    unmatchedB = i_unmatchedReport(ub, unitsB, countsB, agreement', coincidence', unitsA, 'B');

    % Summary.
    summary = struct();
    summary.nUnitsA        = nA;
    summary.nUnitsB        = nB;
    summary.nMatched       = numel(mi);
    summary.meanAgreement  = i_safeMean([matches.agreement]);
    summary.window_ms      = opts.window_ms;
    summary.matchThreshold = opts.match_threshold;
    % Spike-level coverage of A by B (fraction of A spikes with a coincident B).
    summary.spikeCoverageA = i_totalCoverage(spkA, spkB, win);
    summary.spikeCoverageB = i_totalCoverage(spkB, spkA, win);

    % Exact-match check (deterministic-change validation).
    exact = struct('same_spike_times', false, 'same_labels', [], 'label_agreement', []);
    if numel(spkA) == numel(spkB)
        [sa, ia] = sort(spkA(:)); [sb, ib] = sort(spkB(:));
        exact.same_spike_times = isequal(sa, sb);
        if exact.same_spike_times
            la = lblA(ia); lb = lblB(ib);
            exact.label_agreement = mean(double(la(:)) == double(lb(:)));
            exact.same_labels     = (exact.label_agreement == 1);
        end
    end

    result = struct();
    result.unitsA = unitsA; result.unitsB = unitsB;
    result.countsA = countsA; result.countsB = countsB;
    result.agreement = agreement;
    result.coincidence = coincidence;
    result.matches = matches;
    result.unmatchedA = unmatchedA;
    result.unmatchedB = unmatchedB;
    result.summary = summary;
    result.exact = exact;
    result.nameA = opts.nameA; result.nameB = opts.nameB;

    if opts.verbose
        i_printReport(result);
    end
end

% ------------------------------------------------------------------------
function [spk, lbl] = i_getSorting(X)
    if ischar(X) || isstring(X)
        s = kiaSort_load_results(char(X));
        spk = s.spike_idx; lbl = s.unifiedLabels;
    elseif isstruct(X)
        if isfield(X, 'spike_idx') && isfield(X, 'unifiedLabels')
            spk = X.spike_idx; lbl = X.unifiedLabels;
        elseif isfield(X, 'spikeIdx') && isfield(X, 'labels')
            spk = X.spikeIdx; lbl = X.labels;
        else
            error('kiaSort_compare_sortings:input', ...
                ['A/B struct must have fields (spike_idx & unifiedLabels) ' ...
                 'or (spikeIdx & labels).']);
        end
    else
        error('kiaSort_compare_sortings:input', ...
            'A/B must be a results folder path or a sorting struct.');
    end
    spk = double(spk(:));
    lbl = double(lbl(:));
    if numel(spk) ~= numel(lbl)
        error('kiaSort_compare_sortings:input', ...
            'spike_idx and labels must have the same length.');
    end
end

% ------------------------------------------------------------------------
function m = i_countCoincident(a, b, win)
% Number of spikes in A that have a spike in B within +/- win samples, with
% each B spike matched at most once (two-pointer over two sorted vectors).
    m = 0;
    ia = 1; ib = 1;
    na = numel(a); nb = numel(b);
    while ia <= na && ib <= nb
        d = a(ia) - b(ib);
        if abs(d) <= win
            m = m + 1; ia = ia + 1; ib = ib + 1;   % consume both
        elseif d > win
            ib = ib + 1;                            % b too early
        else
            ia = ia + 1;                            % a too early
        end
    end
end

% ------------------------------------------------------------------------
function cov = i_totalCoverage(a, b, win)
% Fraction of spikes in A (any label) that have a coincident spike in B.
    if isempty(a), cov = NaN; return; end
    a = sort(a(:)); b = sort(b(:));
    m = i_countCoincident(a, b, win);
    cov = m / numel(a);
end

% ------------------------------------------------------------------------
function rep = i_unmatchedReport(idx, units, counts, agreement, coincidence, otherUnits, side) %#ok<INUSD>
% Best partial overlap for each unmatched unit, plus a split/merge hint:
% how many of the OTHER sorting's units capture >10% of this unit's spikes.
    rep = struct('unit',{},'nSpikes',{},'bestOther',{},'bestAgreement',{}, ...
        'nContributors',{},'coverage',{});
    for t = 1:numel(idx)
        u = idx(t);
        [bestAg, bestj] = max(agreement(u,:));
        if isempty(bestAg), bestAg = 0; bestj = []; end
        contribFrac = coincidence(u,:) / max(1, counts(u));
        nContrib = sum(contribFrac > 0.10);
        rep(t).unit          = units(u);
        rep(t).nSpikes       = counts(u);
        if ~isempty(bestj)
            rep(t).bestOther = otherUnits(bestj);
        else
            rep(t).bestOther = NaN;
        end
        rep(t).bestAgreement = bestAg;
        rep(t).nContributors = nContrib;               % >1 suggests a split
        rep(t).coverage      = sum(coincidence(u,:)) / max(1, counts(u));
    end
end

% ------------------------------------------------------------------------
function m = i_safeMean(v)
    if isempty(v), m = NaN; else, m = mean(v); end
end

% ------------------------------------------------------------------------
function opts = i_parseOpts(opts, args)
    if isempty(args), return; end
    if mod(numel(args),2) ~= 0
        error('kiaSort_compare_sortings:opts', 'Options must be name/value pairs.');
    end
    for i = 1:2:numel(args)
        name = args{i};
        if ~isfield(opts, name)
            error('kiaSort_compare_sortings:opts', 'Unknown option "%s".', name);
        end
        opts.(name) = args{i+1};
    end
end

% ------------------------------------------------------------------------
function i_printReport(r)
    fprintf('\n=== kiaSort sorting comparison: %s vs %s ===\n', r.nameA, r.nameB);
    fprintf('  window = %.2f ms, match threshold = %.2f\n', ...
        r.summary.window_ms, r.summary.matchThreshold);
    fprintf('  %s: %d units, %d spikes | %s: %d units, %d spikes\n', ...
        r.nameA, r.summary.nUnitsA, sum(r.countsA), ...
        r.nameB, r.summary.nUnitsB, sum(r.countsB));

    if r.exact.same_spike_times
        if isequal(r.exact.same_labels, true)
            fprintf('  EXACT MATCH: identical spike times and identical labels.\n');
        else
            fprintf('  Identical spike times; label agreement = %.4f%s\n', ...
                r.exact.label_agreement, i_relabelNote(r.exact.label_agreement));
        end
    end

    fprintf('  matched units (agreement >= %.2f): %d / %d\n', ...
        r.summary.matchThreshold, r.summary.nMatched, r.summary.nUnitsA);
    fprintf('  mean matched agreement: %.4f\n', r.summary.meanAgreement);
    fprintf('  spike coverage: %s->%s %.4f, %s->%s %.4f\n', ...
        r.nameA, r.nameB, r.summary.spikeCoverageA, ...
        r.nameB, r.nameA, r.summary.spikeCoverageB);

    if ~isempty(r.matches)
        fprintf('\n  %-10s %-10s %-10s %-9s %-9s %-9s\n', ...
            [r.nameA '#'], [r.nameB '#'], 'agree', 'nMatch', 'onlyA', 'onlyB');
        for t = 1:numel(r.matches)
            m = r.matches(t);
            fprintf('  %-10g %-10g %-10.4f %-9d %-9d %-9d\n', ...
                m.unitA, m.unitB, m.agreement, m.nMatch, m.nOnlyA, m.nOnlyB);
        end
    end

    if ~isempty(r.unmatchedA)
        fprintf('\n  unmatched %s units (missed, if %s is truth):\n', r.nameA, r.nameA);
        i_printUnmatched(r.unmatchedA, r.nameB);
    end
    if ~isempty(r.unmatchedB)
        fprintf('\n  unmatched %s units (false positives, if %s is truth):\n', r.nameB, r.nameA);
        i_printUnmatched(r.unmatchedB, r.nameA);
    end
    fprintf('\n');
end

function i_printUnmatched(rep, otherName)
    for t = 1:numel(rep)
        u = rep(t);
        hint = '';
        if u.nContributors > 1
            hint = sprintf(' [split across %d %s units]', u.nContributors, otherName);
        end
        fprintf('    unit %g: %d spikes, best %s unit %g (agree %.3f, coverage %.3f)%s\n', ...
            u.unit, u.nSpikes, otherName, u.bestOther, u.bestAgreement, u.coverage, hint);
    end
end

function s = i_relabelNote(agr)
    if agr < 1
        s = ' (spike partition may match under a unit-id relabelling)';
    else
        s = '';
    end
end
