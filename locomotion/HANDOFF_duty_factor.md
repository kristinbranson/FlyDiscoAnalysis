# Handoff: Add duty factor to pipeline

## Task

Add duty factor as a new per-step metric in the LimbBoutAnalyzer pipeline, wired into all outputs.

## Definition

`duty_factor = stance_duration / (stance_duration + swing_duration)` per step

- Scalar in [0,1]. Value ~0.6-0.7 at slow speeds, decreasing toward ~0.5 at fast speeds (Mendes et al. 2013).
- Linear (not circular) aggregation: regular mean/std.

## Where to compute

**`locomotion/computeStepFeatures.m`** — this function already receives `step_t0s, step_t1s, stance_t0s, stance_t1s` and `currfly_timestamps`.

Currently computes step duration (L22-24):
```matlab
[durations_frames,durations_time] = computeBoutDurations(step_t0s,step_t1s,currfly_timestamps);
```

Stance duration is computable from `stance_t0s, stance_t1s`. Swing duration = step_duration - stance_duration.

Add after the existing duration computation (~L24):
```matlab
[stance_durations_frames, stance_durations_time] = computeBoutDurations(stance_t0s, stance_t1s, currfly_timestamps);
swing_durations_time = durations_time - stance_durations_time;
boutfeatures.duty_factor = stance_durations_time ./ durations_time;
```

**Edge cases to handle:**
- `durations_time == 0` → duty_factor = NaN
- Very short steps where swing ≈ 0 → duty_factor ≈ 1 (valid)
- NaN timestamps → NaN duration → NaN duty_factor

## Where to wire

### 1. `computeboutmetrics2.m`
Step features from `computeStepFeatures` flow through `computeboutmetrics2.m` into `bout_metrics.perfly(fly).perlimb(limb).step.stepfeatures.duty_factor`. This should happen automatically since `computeStepFeatures` returns all fields in `boutfeatures` and they get stored as `.stepfeatures`.

Verify: check how stepfeatures are stored in `computeboutmetrics2.m` and whether all fields propagate.

### 2. `LimbBoutAnalyzer.m` — `combineStepMetrics` (L1340-1410)
Add `duty_factor` to the list of step scalar features that get flattened into locostatsperexp. Look at how `amplitude_BL`, `distance_BL` etc. are handled — duty_factor follows the same pattern.

Output field: `duty_factor__step__LEDon__pair1` (and pair2, pair3, all, limb1-6)

### 3. `LimbBoutAnalyzer.m` — `buildWalkStructForCondition` (L981+)
Add `duty_factor` to step_struct fields alongside existing `step_duration`, `stance_duration`, `swing_duration`.

### 4. `locomotionmetricsswingstanceboutstats.mat`
This is saved by `LimbBoutAnalyzer.saveResults()`. duty_factor should appear in `bout_metrics.perfly(fly).perlimb(limb).step.stepfeatures.duty_factor` automatically if step 1 works.

## Testing

### Unit test: `locomotion/tests/test_duty_factor.m`
1. **Known values:** stance=60ms, swing=40ms, step=100ms → duty_factor = 0.6 exactly
2. **Edge case:** step=0ms → duty_factor = NaN
3. **Edge case:** swing=0ms (all stance) → duty_factor = 1.0
4. **Edge case:** stance=0ms (all swing) → duty_factor = 0.0
5. **Multiple steps:** verify vector operation produces correct per-step values
6. **Aggregation:** pass known duty_factor values through combineStepMetrics, verify mean/std/Z correct
7. **Per-pair aggregation:** verify pair1 = mean of limb1+limb6 duty factors (check pair definition)

### Integration test on test experiment
- Run on `VNC2_JRC_SS57983_RigD_20230913T120134`
- Verify duty_factor in:
  - `locostatsperexp_onfloor.mat` → fields like `duty_factor__step__LEDon__pair1`
  - `locomotion_walkstruct.mat` → `step_struct_ON.duty_factor`, `step_struct_OFF.duty_factor`
  - `locomotionmetricsswingstanceboutstats.mat` → `bout_metrics.perfly(1).perlimb(1).step.stepfeatures.duty_factor`
- Manually compute for one fly/limb: take stance_t0s/t1s and step_t0s/t1s, compute durations, verify ratio matches
- Verify duty_factor decreases with speed: scatter plot duty_factor vs velmag_ctr across steps

## Key files

| File | Lines | Role |
|---|---|---|
| `locomotion/computeStepFeatures.m` | 22-24 | Add duty_factor computation here |
| `locomotion/computeboutmetrics2.m` | ~109-132 | Verify stepfeatures propagation |
| `locomotion/LimbBoutAnalyzer.m` | 1340-1410 | `combineStepMetrics` — add duty_factor |
| `locomotion/LimbBoutAnalyzer.m` | 981-1030 | `buildWalkStructForCondition` step_struct — add field |
| `locomotion/computeBoutDurations.m` | all | Used to compute stance durations from timestamps |

## Test experiment

`/groups/branson/home/robiea/Projects_data/FlyDisco/Bubble_data/20260112_testing_locomotioncomputeperframestats/VNC2_JRC_SS57983_RigD_20230913T120134`

Settings: `settingsdir = '/groups/branson/home/robiea/Code_versioned/BransonFlyDiscoSettings/settings'`, `analysis_protocol = '20260326_flybubble_LED_VNC2'`

## Pair definitions

From `computeboutmetrics2.m`: pairs = [1,6; 2,5; 3,4] = (RF,LF), (RM,LM), (RH,LH)
