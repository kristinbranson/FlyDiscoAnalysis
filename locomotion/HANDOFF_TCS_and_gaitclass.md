# Handoff: Wire TCS into pipeline + test gait_class

## Task

Two items before the full VNC re-run:
1. Wire existing TCS (Tripod Coordination Strength) computation into all pipeline outputs
2. Test that gait_class fractions are correct end-to-end

## Background

TCS and gait_class are both computed per walk in `compute_WalkFeatures.m` but TCS is silently dropped during aggregation. Gait_class is wired through but has never been validated on real data.

## Current state of TCS

**Computation exists and is unit-tested:**
- `locomotion/compute_TCS.m` — computes TCS per Wosnitza 2012 (t2/t1 for tripod swing overlap)
- Called from `locomotion/compute_WalkFeatures.m` L219-221
- Stored as `walkfeaturestruct(ct).TCS` — a struct with:
  - `.tripod_LM`, `.tripod_RM` — per-tripod-group results
  - `.both` — combined, with `.data` (1×n event TCS values in [0,1]), `.mean`, `.std`, `.n`, `.event_times`, `.event_speeds`, `.n_steps`, `.n_nontripod`
- Unit test: `locomotion/tests/test_compute_TCS.m` (8 tests, synthetic data, all pass)

**The gap — dropped at two aggregation stages:**

Stage 1: Per-fly aggregation in `compute_WalkFeatures.m` L249-273:
```matlab
for fld = 1:numel(flds)
    if ~isstruct(...)           % non-struct → concat
    elseif any(strcmp(pfflist...))  % perframe features → computePerFlywalkperframefeatures
    elseif startsWith(flds{fld}, 'phase')  % phase → computePerFlyphasefeatures
    elseif strcmp(flds{fld}, 'gait_class') % gait → computePerFlygaitclass
    end  % TCS is a struct, not in pfflist, not 'phase*', not 'gait_class' → FALLS THROUGH
end
```

Stage 2: Per-exp aggregation in `computeWalkMetrics.m` L51-73 — identical pattern, same gap.

**TCS is NOT in:**
- `walk_metrics.perfly` or `walk_metrics.perexp`
- `locostatsperexp_onfloor.mat`
- `locomotion_walkstruct.mat` walk_struct
- `locomotionmetricsswingstanceboutstats.mat`

## Implementation plan for TCS

### Step 1: Create `locomotion/computePerFlyTCS.m`
Pattern after `computePerFlygaitclass.m`. Input: `walkfeaturestruct` (struct array with `.TCS` per walk). Output: per-fly aggregation struct with:
- `.data` — concatenated `TCS.both.data` across all walks
- `.event_speeds` — concatenated
- `.mean`, `.std`, `.n` — aggregate
- `.n_steps`, `.n_nontripod` — summed

TCS is scalar [0,1], use linear (not circular) mean/std.

### Step 2: Create `locomotion/computePerExpTCS.m`
Input: `perwalk_metrics`, `perfly_metrics`. Output: per-exp aggregation. Concatenate `.data` across flies, compute exp-level mean/std/n.

### Step 3: Add handlers to aggregation loops

In `compute_WalkFeatures.m` after L269 (the gait_class handler):
```matlab
elseif strcmp(flds{fld}, 'TCS')
    perflywalkfeatures.TCS = computePerFlyTCS(walkfeaturestruct);
```

In `computeWalkMetrics.m` after L71:
```matlab
elseif strcmp(flds{fld}, 'TCS')
    perexp_metrics.TCS = computePerExpTCS(perwalk_metrics, perfly_metrics);
```

### Step 4: Add `combineTCSMetrics` to LimbBoutAnalyzer

In `LimbBoutAnalyzer.m`, new private method after `combineGaitMetrics` (L1512):
```matlab
function statsperexp = combineTCSMetrics(~, walk_metrics, led_label, statsperexp)
    funname = sprintf('TCS__walk__%s__all', led_label);
    currstruct = struct;
    if ~isempty(fields(walk_metrics.perexp)) && isfield(walk_metrics.perexp, 'TCS')
        tcs = walk_metrics.perexp.TCS;
        currstruct.mean = tcs.mean;
        currstruct.std = tcs.std;
        currstruct.Z = tcs.n;
    else
        currstruct.mean = NaN;
        currstruct.std = NaN;
        currstruct.Z = 0;
    end
    statsperexp.(funname) = currstruct;
end
```

Call it from `computeStatsPerExp` alongside `combineGaitMetrics`.

### Step 5: Add TCS to `buildWalkStructForCondition`

In the walk_struct field population (L832+), add `TCS` field per walk from `walkfeaturestruct(w).TCS.both.mean` (scalar per walk).

### Step 6: Test

**Unit test (synthetic):**
- Create `locomotion/tests/test_TCS_aggregation.m`
- Build synthetic `walkfeaturestruct` with known TCS values (e.g., walk1: TCS.both.data=[0.8, 0.9], walk2: TCS.both.data=[0.7])
- Pass through `computePerFlyTCS` → verify mean = mean([0.8, 0.9, 0.7]), std matches, n=3
- Pass through `computePerExpTCS` with 2 flies → verify concatenation correct
- Pass through `combineTCSMetrics` → verify field name = `TCS__walk__LEDon__all`, mean/std/Z correct
- Test empty walks (no TCS events) → should produce NaN mean, Z=0

**Integration test:**
- Run full pipeline on test experiment `VNC2_JRC_SS57983_RigD_20230913T120134`
- Check TCS appears in:
  - `locostatsperexp_onfloor.mat` fields `TCS__walk__LEDon__all`, `TCS__walk__LEDoff__all`
  - `locomotion_walkstruct.mat` → `walk_struct_ON.TCS`, `walk_struct_OFF.TCS`
- Manually compute TCS for one walk/fly and compare to pipeline value

## Gait_class testing

Gait_class IS wired through the full pipeline. It needs validation, not implementation.

**What to verify:**
1. Run on test experiment
2. Check `locostatsperexp_onfloor.mat` has fields: `gait_class__walk__LEDon__tripod_frac`, `...__tetrapod_frac`, `...__grounded_frac`, `...__airborne_frac`, `...__other_frac` (and LEDoff variants)
3. Fractions should sum to ~1.0 for each LED condition
4. Manually count tripod/tetrapod frames for one walk and compare to pipeline output
5. Check that `walk_struct` in `locomotion_walkstruct.mat` does NOT contain gait_class (it's currently walk-level only in locostatsperexp, not in walk_struct — verify this is intentional or if it should be added)

**How gait_class flows:**
- `compute_gait_class.m` → per-frame codes (1-5)
- `compute_WalkFeatures.m` L226-235 → per-walk counts stored in `walkfeaturestruct(ct).gait_class`
- `computePerFlygaitclass.m` → per-fly counts
- `computePerExpgaitclass.m` → per-exp counts
- `combineGaitMetrics()` (LimbBoutAnalyzer L1474-1512) → locostatsperexp fractions

## Key files

| File | Lines | Role |
|---|---|---|
| `locomotion/compute_TCS.m` | all | TCS computation (DO NOT MODIFY) |
| `locomotion/tests/test_compute_TCS.m` | all | TCS unit test (existing, passes) |
| `locomotion/compute_WalkFeatures.m` | 219-221, 249-273 | TCS call + per-fly aggregation gap |
| `locomotion/computeWalkMetrics.m` | 51-73 | Per-exp aggregation gap |
| `locomotion/LimbBoutAnalyzer.m` | 425-487, 1412-1512 | `computeStatsPerExp`, combine methods |
| `locomotion/LimbBoutAnalyzer.m` | 791-1267 | `buildWalkStructForCondition` |
| `locomotion/computePerFlygaitclass.m` | all | Pattern to follow for TCS per-fly |
| `locomotion/computePerExpgaitclass.m` | all | Pattern to follow for TCS per-exp |
| `locomotion/compute_gait_class.m` | all | Gait classification (DO NOT MODIFY) |

## Test experiment

`/groups/branson/home/robiea/Projects_data/FlyDisco/Bubble_data/20260112_testing_locomotioncomputeperframestats/VNC2_JRC_SS57983_RigD_20230913T120134`

Settings: `settingsdir = '/groups/branson/home/robiea/Code_versioned/BransonFlyDiscoSettings/settings'`, `analysis_protocol = '20260326_flybubble_LED_VNC2'`

## Order

1. Implement TCS wiring (steps 1-5)
2. Write TCS unit test + run
3. Integration test TCS on test experiment
4. Test gait_class on same experiment
5. Verify both in all output files
