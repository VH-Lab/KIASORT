# Compare

Quantitatively compare two KIASORT spike sortings — or a sorting against
ground truth — instead of eyeballing the curation GUI.

## `kiaSort_compare_sortings`

Matches the units of two sortings by spike-train coincidence and reports
per-unit agreement, matched/split/orphaned units, and overall coverage.
Comparing two sorts and comparing against ground truth are the same
operation: pass the ground truth as one side.

```matlab
% Two independent runs (folders each containing a RES_Sorted subfolder):
r = kiaSort_compare_sortings(runDirA, runDirB, 30000);   % fs = 30 kHz

% Ground truth (spike samples + true unit ids) vs a sort:
gt = struct('spikeIdx', gtSamples, 'labels', gtUnitIds);
r  = kiaSort_compare_sortings(gt, sortDir, 30000, ...
        'nameA','truth', 'nameB','kiaSort');
```

Inputs (`A`, `B`) may each be a results-folder path (loaded via
`kiaSort_load_results`), a struct with `.spike_idx`/`.unifiedLabels`, or a
struct with `.spikeIdx`/`.labels`. `fs` is the sampling rate in Hz.

The agreement score for a unit pair is the SpikeInterface convention
`n_match / (n_A + n_B - n_match)` (i.e. `TP/(TP+FN+FP)`), with spikes counted
as coincident within `window_ms` (default 0.4 ms) and each spike matched at
most once.

### Validating a deterministic code change

For a change that should NOT alter results (e.g. parallelising stage 3
`sortData`), reuse the **same** stage-1/stage-2 output for both runs and expect
an exact match — the randomness in KIASORT lives in stage 2's UMAP, not in
stage 3's template matching:

```matlab
r = kiaSort_compare_sortings(serialDir, parallelDir, 30000);
assert(r.exact.same_spike_times && r.exact.same_labels)
```

Options: `window_ms` (0.4), `match_threshold` (0.5), `ignore_labels` (e.g.
`0` to drop a noise label), `nameA`/`nameB`, `verbose`.
