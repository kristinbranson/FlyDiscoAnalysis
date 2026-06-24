# Handoff: Add metachronal lag to pipeline

## Task

Add metachronal lag as a new per-walk metric in the LimbBoutAnalyzer pipeline, wired into all outputs.

## Definition

Metachronal lag = time delay (ms) between consecutive ipsilateral swing onsets in the posterior-to-anterior direction (Strauss & Heisenberg 1990).

Four ipsilateral pairs, posterior leg is the reference:
- `RH_to_RM`: RH swing onset to next RM swing onset
- `RM_to_RF`: RM swing onset to next RF swing onset
- `LH_to_LM`: LH swing onset to next LM swing onset
- `LM_to_LF`: LM swing onset to next LF swing onset

Leg indices in `legorder = {'RF','RM','RH','LH','LM','LF'}` (1-6):
- RH→RM = leg 3 → leg 2
- RM→RF = leg 2 → leg 1
- LH→LM = leg 4 → leg 5
- LM→LF = leg 5 → leg 6

Per-walk value = mean lag (ms) across all paired events in that walk. Also compute a pooled `ipsi_all` = mean of all 4 pairs.

Linear scalar (milliseconds), standard mean/std aggregation. Positive = posterior leads (normal metachronal wave). Negative = reversed.

## Existing code to reuse

### `computePhaseLag.m` + `findClosestSteps.m`

Already implement cross-limb swing onset pairing. Called from `compute_WalkFeatures.m` L190-192.

**`findClosestSteps(refdata, stepdata, pre_pad_proportion)`** (in `locomotion/`):
- `refdata`: swing onset times for the reference (posterior) limb within a walk
- `stepdata`: cell(nlimbs) of swing onset times for all limbs
- For each reference swing onset, finds the closest swing onset on each other limb within a window of `[refstep - frame_window, refstep + period)` where `period = diff(refdata)` (the reference limb's step period)
- Returns `matched_stepdata(nlimbs, nevents-1)`: matched onset frame for each limb and each reference event
- Returns NaN if no swing onset found in the window

**Pairing rule**: the anterior leg's swing onset must occur within one reference-leg step cycle. This naturally enforces "leg 2 onset before leg 1's next onset." The `findClosestSteps` search window `[refstep, refstep + period)` achieves this (with `pre_pad_proportion = 0`).

**`computePhaseLag.m`** does the same pairing but converts frame differences to phase (radians) and uses circular stats. Metachronal lag keeps it in milliseconds, uses linear stats.

### `compute_TCS.m`

Different pairing strategy (lead-leg anchored, tripod-grouped). Not directly reusable, but the wiring pattern through `compute_WalkFeatures` → `computePerFlyTCS` → `computePerExpTCS` → `combineTCSMetrics` is a good template for the aggregation chain.

## Where to compute

### New function: `locomotion/computeMetachronalLag.m`

```matlab
function mcl = computeMetachronalLag(currflyboutdata, walk_t0, walk_t1, timestamps)
```

**Inputs:**
- `currflyboutdata`: from `obj.limbBoutData(fly)`, has `.perlimb(limb).swing.start_indices`
- `walk_t0`, `walk_t1`: walk bout boundaries (trajectory frame indices)
- `timestamps`: per-fly timestamps from `trx.movie_timestamps{1}(firstframe:endframe)`

**Logic** (for each of the 4 ipsilateral pairs):
1. Get swing onset times for the posterior (reference) limb within `[walk_t0, walk_t1]`
2. Get swing onset times for the anterior (target) limb (all of them — `findClosestSteps` does the windowing)
3. Call `findClosestSteps(ref_swing_t0s, stepdata, 0)` with `pre_pad_proportion = 0` (no backward search, only forward within one period)
4. Extract matched onset frames for the target limb
5. Compute lag in ms: `(timestamps(matched_target) - timestamps(ref_swing_t0s(1:end-1))) * 1000`
6. Store per-event lags, mean, std, n

**Output struct:**
```matlab
mcl.RH_to_RM.data   = [lag1, lag2, ...] % per-event lags in ms
mcl.RH_to_RM.mean   = mean(data, 'omitnan')
mcl.RH_to_RM.std    = std(data, 'omitnan')
mcl.RH_to_RM.n      = sum(~isnan(data))
% ... same for RM_to_RF, LH_to_LM, LM_to_LF
mcl.ipsi_all.mean    = mean of all 4 pair means
mcl.ipsi_all.n       = sum of all 4 pair n's
```

**Edge cases:**
- Fewer than 2 swing onsets for the reference limb in a walk → NaN for that pair
- No matching anterior onset within the period window → NaN for that event (from `findClosestSteps`)
- Very short walks → likely NaN for all pairs

### Note on timestamps

`compute_WalkFeatures.m` does not currently receive timestamps. Either:
- (a) Pass `timestamps` from `computeWalkMetrics` through to `compute_WalkFeatures` (adds a parameter), OR
- (b) Compute timestamps inside `computeMetachronalLag` from `obj.trx` (pass the LimbBoutAnalyzer object or trx), OR
- (c) Return lag in frames, convert to ms during aggregation

Option (a) is cleanest. `computeWalkMetrics.m` already has access to `obj.trx` and `trx.movie_timestamps`. Add timestamps as a parameter to `compute_WalkFeatures` and pass through. Existing callers would need updating.

Check how `computeStepFeatures.m` gets timestamps — it receives `currfly_timestamps` from `computeboutmetrics2.m` L58:
```matlab
timestamps = trx.movie_timestamps{1};
currfly_timestamps = timestamps(trx.firstframes(fly):trx.endframes(fly));
```

Same pattern would work in `compute_WalkFeatures`.

## Where to wire

### 1. `compute_WalkFeatures.m`

Call `computeMetachronalLag` per walk, store in `walkfeaturestruct(ct).metachronal_lag`. Add after the TCS block (~L221):

```matlab
mcl = computeMetachronalLag(currflyboutdata, walk_t0, walk_t1, currfly_timestamps);
walkfeaturestruct(ct).metachronal_lag = mcl;
```

In the per-fly aggregation loop (~L270), add a branch:
```matlab
elseif strcmp(flds{fld}, 'metachronal_lag')
    perflywalkfeatures.metachronal_lag = computePerFlyMCL(walkfeaturestruct);
```

### 2. New aggregation functions

**`computePerFlyMCL.m`**: Pool per-walk MCL across walks for one fly. For each pair, concatenate `.data` arrays across walks, compute fly-level mean/std/n. Also compute `frm_mean_fly` (mean of pooled events) and `walk_mean_fly` (mean of per-walk means). Same dual-aggregation pattern as phase features in `computePerFlyphasefeatures.m`.

**`computePerExpMCL.m`**: Pool per-fly MCL across flies for one experiment. Weight by n (number of events per fly). Same pattern as `computePerExpTCS.m`.

### 3. `computeWalkMetrics.m`

Add branch in per-exp aggregation loop (after TCS branch):
```matlab
elseif strcmp(flds{fld}, 'metachronal_lag')
    perexp_metrics.metachronal_lag = computePerExpMCL(perwalk_metrics, perfly_metrics);
```

### 4. `LimbBoutAnalyzer.m` — `computeStatsPerExp`

Add `combineMCLMetrics` call after `combineTCSMetrics`:
```matlab
statsperexp = obj.combineMCLMetrics(walk_metrics, cm.label, statsperexp);
```

**New private method `combineMCLMetrics`**: Flatten into locostatsperexp fields:
```
metachronal_lag__walk__LEDon__RH_to_RM   (.mean, .std, .Z)
metachronal_lag__walk__LEDon__RM_to_RF
metachronal_lag__walk__LEDon__LH_to_LM
metachronal_lag__walk__LEDon__LM_to_LF
metachronal_lag__walk__LEDon__ipsi_all
```
(and same for LEDoff)

### 5. `LimbBoutAnalyzer.m` — `buildWalkStructForCondition`

Add MCL to walk_struct, after the TCS block. For each walk, store the per-pair means:
```matlab
walk_struct.metachronal_lag_RH_to_RM(w) = pw(w).metachronal_lag.RH_to_RM.mean;
% ... etc for each pair
```

## Testing

### Unit test: `test_metachronal_lag.m`

1. **Synthetic data with known lag**: Create ground contact signals for 3 ipsilateral limbs with fixed offsets (e.g., hind leads middle by 10ms, middle leads front by 10ms). Verify computed lag matches.
2. **Zero lag (synchronous)**: All limbs swing at the same time → lag = 0.
3. **Reversed wave**: Front leads hind → negative lag.
4. **Single swing onset**: Only 1 swing in walk → NaN (need >= 2 for period).
5. **Missing match**: Anterior limb doesn't swing during reference period → NaN for that event.

### Integration test: add to `test_new_metrics_integration.m`

6. MCL fields exist in locostatsperexp (5 pairs x 2 LED conditions = 10 fields)
7. MCL in walk_struct (5 fields per walk)
8. MCL values are reasonable (positive, typically 5-30ms for Drosophila at walking speeds)
9. Spot-check: recompute MCL for one walk from raw swing onsets, compare to pipeline output

### Existing tests to verify no regression

- `test_TCS_endtoend.m` — exercises full pipeline including `compute_WalkFeatures`
- `test_duty_factor.m` — exercises bout metrics path
- `test_new_metrics_integration.m` — exercises all outputs

## Overlap with other files

- `compute_WalkFeatures.m`: adding MCL call + per-fly aggregation branch. TCS and gait_class already have similar blocks — MCL follows the same pattern.
- `computeWalkMetrics.m`: adding per-exp aggregation branch. Same pattern as TCS.
- `LimbBoutAnalyzer.m`: adding `combineMCLMetrics` method + call in `computeStatsPerExp` + walk_struct fields. Same regions as TCS edits.

No overlap with duty_factor code (different pipeline path: bout vs walk).

## Reference

- Strauss & Heisenberg 1990: metachronal lag definition, power-law relationship with period
- Wosnitza et al. 2012: step cycle definitions, TCS
- DeAngelis et al. 2019: ipsilateral coupling approaches tripod only at fast speeds
- `locomotion_metric_primer_v3.md`: literature review with code-vs-paper comparisons
