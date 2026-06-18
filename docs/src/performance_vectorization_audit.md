# Performance and vectorization audit (GPU-ready style)

This note focuses on hot paths that can be expressed with standard array operations (`reshape`, `sum`, `roll`, elementwise products, matrix products) so they can be migrated later to array backends such as CuPy or Torch with minimal rewrite.

## Scope

- Included: core engine paths (`beads.py`, `barostats.py`, `properties.py`).
- Excluded intentionally: other specialized workflows outside the core profiling set used here.

## Refactor principles for backend portability

1. Prefer whole-array expressions over Python loops.
2. Prefer pre-shaped `(N, 3)` views for Cartesian operations.
3. Avoid Python container construction (`list`, `set`) inside hot loops.
4. Keep logic in terms of backend-agnostic primitives (`roll`, `sum`, broadcasted `*`, `@`/dot-like operations).

## Updated hotspots

### 1) Ring-polymer spring estimators (`ipi/engine/beads.py`)

- `get_vpath` now computes bead differences with one `np.roll` and one reduction.
- `get_fpath` now uses the same `dq` and a shifted subtraction to accumulate pair contributions.

These forms map directly to CuPy/Torch tensor ops (`roll`, elementwise arithmetic, reductions).

### 2) Kinetic stress estimators (`ipi/engine/barostats.py`)

- `get_kstress` and `kstress_mts_sc` now use:
  - `(q - qc).reshape(-1, 3)`
  - force reshape to `(-1, 3)`
  - tensor contraction via existing `noddot`
  - diagonal kinetic term via vectorized component operations

This removes triple nested Python loops and aligns implementation style with `kstress_mts`.

### 3) Bead-mass consistency check (`ipi/engine/properties.py`)

- `get_kinmd` now performs one vectorized consistency check (`np.allclose(dm3, dm3[0])`) for the bead-selected branch.
- Removes repeated per-atom list/set construction.

## Profiling-guided optimization plan

### Baseline summary (from internal profiling runs)

Representative runs (`classical_md_direct`, `pimd_npt_noff`) indicate most cost in non-forcefield overhead is concentrated in:

1. MTS integrator recursion and orchestration (`dynamics.py`: `step`, `mtsprop`, `mtsprop_ba`, `mtsprop_ab`, `pstep`).
2. Barostat stepping (`barostats.py`: `pstep`, `qcstep`, `stress_mts`/`kstress_mts` chain).
3. Thermostat kernel (`thermostats.py`: Langevin `step`).
4. Free ring-polymer propagation (`normalmodes.py`: `free_qstep` path).

### Phase 1 — Low risk, high return (1–2 PRs)

1. **Reduce repeated `dstrip` and temporary creation in hot loops**
   - Cache local stripped arrays inside `barostats.pstep`, `barostats.qcstep`, and integrator-level stepping.
   - Avoid repeated reshape/reindex operations if unchanged inside an MTS sub-step.
2. **Vectorize remaining barostat stress helper paths**
   - Bring all SC/MTS stress paths to one shared vectorized contraction utility.
3. **Micro-optimize thermostat step**
   - Minimize temporary vectors in `ThermoLangevin.step` by reusing buffers where possible.

**Success criteria:** >=10% reduction in wall-time for `pimd_npt_noff` profile case with unchanged energies/trajectories within tolerance.

### Phase 2 — Structural cleanup of integrator overhead (2–3 PRs)

1. **Flatten MTS recursion overhead**
   - Evaluate iterative schedule construction for fixed `nmts` to reduce Python call overhead from `mtsprop_*` recursion.
2. **Consolidate repeated `pconstraints`/small-step dispatches**
   - Group no-op or identical operations in inner loops.
3. **Precompute level-dependent constants**
   - Store frequently reused `pdt`-scaled factors and masks once per step.

**Success criteria:** additional >=10–15% reduction in MD-step overhead for profiling inputs.

### Phase 3 — Backend portability layer (NumPy -> CuPy/Torch ready)

1. Introduce an internal array namespace shim (`xp`) for selected kernels (`roll`, `sum`, `exp`, `reshape`, contraction).
2. Provide backend-aware `noddot` / contraction wrappers.
3. Keep conversions centralized at explicit boundaries (`dstrip`/input-output marshalling).

**Success criteria:** ability to run selected kernels with NumPy and CuPy/Torch-compatible arrays behind one API.

### Validation protocol for each phase

- Re-run lightweight profiles:
  - `ipi_tests/profiling/classical_md_direct`
  - `ipi_tests/profiling/pimd_npt_noff`
- Run existing unit/regression subset relevant to motion/barostat/thermostat.
- Track:
  - wall time per step,
  - cumulative time of top 20 functions,
  - numerical drift of conserved quantities.

## Migration thoughts (NumPy -> CuPy/Torch)

- The rewritten sections are now mostly backend-friendly primitives.
- A practical next step is introducing a tiny array-backend adapter (e.g., `xp`) for `roll/sum/reshape` and a backend-aware `noddot` implementation.
- Keep `dstrip` boundaries explicit; these are likely spots where backend conversion decisions should be centralized.

## Notes

- This pass is focused on the core runtime paths covered by the profiling setups used here.
- Optimization should stay profiling-driven, with measurable before/after metrics at each phase.
