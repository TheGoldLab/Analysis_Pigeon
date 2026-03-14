"""
Simulation utilities for the Pigeon project.
Python translations of getPigeon_simulatedData.m and getPigeon_simulatedDataTable.m.

Generative model
----------------
Each trial is a 1-D random walk starting at 0.  At every step the position
increments by a sample from N(generative_mean, generative_std).  The trial
ends (boundary crossing) when the absolute position first exceeds a
threshold ``bound``.  An NDT (non-decision time) of additional steps is
added after the crossing to produce the final RT.

Key parameters
--------------
generative_mean : float
    Drift rate.  SNR = generative_mean / generative_std.
generative_std : float
    Per-step noise SD (shared across SNR levels within a simulation).
bound_mean : float
    Mean of the per-trial bound distribution.
bound_std : float
    Trial-to-trial SD of the bound (0 = fixed bound).
ndt_min, ndt_max : int
    Non-decision time drawn uniformly from [ndt_min, ndt_max].
    NDT formula matches MATLAB: randi(hi-lo+1) - lo + 1,
    i.e. Python: rng.integers(1, hi-lo+2) + (1-lo).

Constants
---------
MIN_BOUND : float
    Floor on the sampled bound (prevents degenerate zero-bound trials).
BLOCK_DEFAULTS : dict
    Per-block step/coin budget defaults, used by ``get_simulated_data_table``
    to truncate simulated blocks to the same step budget as real data.

Usage notes
-----------
- For bias-correction simulations (``BiasCorrection.ipynb``), call
  ``simulate_trials`` directly with ``num_trials=1_000_000``.  Do NOT use
  ``get_simulated_data_table`` for this — it truncates trials via
  ``BLOCK_DEFAULTS['steps_per_block'] = 600``.
- ``get_simulated_data_table`` is designed for fitting/plotting simulations
  that need the same DataFrame format as ``get_data_table``.
"""

import os
import numpy as np
import pandas as pd

from .data import WALL, get_bounds, apply_bias_corrections
from .pathDefaults import BIAS_CORRECTION_FILE


MIN_BOUND = 0.03  # minimum allowed bound (prevents degenerate zero-bound trials)

# Per-block reward/penalty defaults (1-indexed block numbers, matching MATLAB blockDefaults)
BLOCK_DEFAULTS = {
    1: dict(steps_per_block=600, p2p_coins=0, coins_gained=1, coins_lost=0,  steps_lost=0),
    2: dict(steps_per_block=600, p2p_coins=0, coins_gained=1, coins_lost=4,  steps_lost=0),
    3: dict(steps_per_block=600, p2p_coins=0, coins_gained=1, coins_lost=0,  steps_lost=30),
    4: dict(steps_per_block=600, p2p_coins=0, coins_gained=1, coins_lost=0,  steps_lost=0),
    5: dict(steps_per_block=600, p2p_coins=0, coins_gained=1, coins_lost=4,  steps_lost=0),
    6: dict(steps_per_block=600, p2p_coins=0, coins_gained=1, coins_lost=0,  steps_lost=30),
    7: dict(steps_per_block=600, p2p_coins=0, coins_gained=1, coins_lost=1,  steps_lost=0),
}


def simulate_trials(
    generative_mean=0.05,
    generative_std=0.15,
    num_trials=1_000_000,
    max_steps=50,
    bound_mean=0.0,
    bound_std=0.0,
    bound_slope=0.0,
    bound_max=0.75,
    ndt_min=1,
    ndt_max=3,
    lapse_rate=0.0,
    bounds=None,
    rng=None,
    return_steps=True,
):
    """
    Simulate Pigeon task trials (core engine).
    Translated from getPigeon_simulatedData.m.

    Parameters
    ----------
    generative_mean : float or array-like
        Mean step size per SNR level. Multiple values are randomly interleaved
        across trials.
    generative_std : float
        SD of each step (shared across all SNR levels).
    num_trials : int
        Number of trials to simulate.
    max_steps : int
        Maximum steps per trial; trial ends even without a boundary crossing.
    bound_mean : float or array-like
        Mean absolute bound. If length == len(generative_mean), each SNR level
        gets its own bound mean.
    bound_std : float or array-like
        Trial-to-trial SD of bound (0 = fixed bound).
    bound_slope : float
        Multiplicative time-varying slope applied linearly over the first 10
        steps, then held constant.
    bound_max : float
        Screen-edge limit (task geometry: 0.75).
    ndt_min, ndt_max : int or array-like
        Non-decision time range (steps). Drawn uniformly from [ndt_min, ndt_max].
        If arrays, length must equal len(generative_mean) for per-SNR NDTs.
    lapse_rate : float
        Fraction of trials assigned a random choice (ignores DV).
    bounds : array-like or None
        If given, per-trial bound is drawn randomly from this pool (overrides
        bound_mean / bound_std).
    rng : np.random.Generator or None

    Returns
    -------
    choices   : (num_trials,) float array  — 0 or 1 (NaN if no crossing)
    rts       : (num_trials,) int array    — total steps including NDT
    DV        : (num_trials, max_steps+1)  — cumulative step sequences
    bounds_out: (num_trials,) float array  — bound at crossing
    ndts      : (num_trials,) int array    — NDT per trial
    snrs      : (num_trials,) float array  — SNR per trial
    steps     : list of num_trials arrays  — step sequence per trial
    """
    if rng is None:
        rng = np.random.default_rng()

    generative_mean = np.atleast_1d(np.asarray(generative_mean, dtype=float))
    n_means = len(generative_mean)

    # Assign each trial a generative mean (interleave randomly if multiple)
    if n_means > 1:
        i_means = rng.integers(0, n_means, size=num_trials)
        trial_means = generative_mean[i_means]
    else:
        i_means = np.zeros(num_trials, dtype=int)
        trial_means = np.full(num_trials, generative_mean[0])

    # Draw NDTs matching MATLAB: randi(NDTMax-NDTMin+1) - NDTMin + 1
    # This produces values in {2-NDTMin, ..., NDTMax-NDTMin+2-NDTMin} = {2-lo, ..., hi-2lo+2}
    # For NDTMin=2, NDTMax=3 (default): ndts ∈ {0, 1}
    # For NDTMin=1, NDTMax=3: ndts ∈ {1, 2, 3}
    ndt_min_arr = np.atleast_1d(np.asarray(ndt_min, dtype=int))
    ndt_max_arr = np.atleast_1d(np.asarray(ndt_max, dtype=int))
    ndts = np.empty(num_trials, dtype=int)
    if len(ndt_min_arr) == 1 or len(ndt_min_arr) != n_means:
        lo, hi = int(ndt_min_arr[0]), int(ndt_max_arr[0])
        ndts[:] = rng.integers(1, hi - lo + 2, size=num_trials) + (1 - lo)
    else:
        for mm in range(n_means):
            lm = i_means == mm
            lo, hi = int(ndt_min_arr[mm]), int(ndt_max_arr[mm])
            ndts[lm] = rng.integers(1, hi - lo + 2, size=int(lm.sum())) + (1 - lo)

    # Generate DV: starts at 0, each step ~ N(trial_mean, generative_std)
    noise = rng.normal(trial_means[:, None], generative_std,
                       size=(num_trials, max_steps))
    DV = np.concatenate(
        [np.zeros((num_trials, 1)), np.cumsum(noise, axis=1)],
        axis=1,
    )  # shape (num_trials, max_steps + 1)

    # Build per-trial bound values
    bound_mean_arr = np.atleast_1d(np.asarray(bound_mean, dtype=float))
    bound_std_arr  = np.atleast_1d(np.asarray(bound_std,  dtype=float))

    if bounds is not None:
        bounds_arr      = np.asarray(bounds, dtype=float)
        bound_per_trial = bounds_arr[rng.integers(0, len(bounds_arr), size=num_trials)]
    elif len(bound_mean_arr) > 1 and len(bound_mean_arr) == n_means:
        bm = bound_mean_arr[i_means]
        bs = (bound_std_arr[i_means] if len(bound_std_arr) == n_means
              else np.full(num_trials, float(bound_std_arr[0])))
        bound_per_trial = bm.copy()
        variable = bs > 0
        if variable.any():
            bound_per_trial[variable] = np.maximum(
                MIN_BOUND, rng.normal(bm[variable], bs[variable]))
    else:
        bm = float(bound_mean_arr[0])
        bs = float(bound_std_arr[0])
        bound_per_trial = (np.maximum(MIN_BOUND, rng.normal(bm, bs, size=num_trials))
                           if bs > 0 else np.full(num_trials, bm))

    # Build slope vector for time-varying bound (covers columns 0..max_steps-1,
    # the only columns checked for crossings). Skipped when bound_slope == 0.
    if bound_slope != 0:
        n_col = max_steps + 1
        slope_vec = np.concatenate([
            np.linspace(1.0, bound_slope, min(10, n_col)),
            np.full(max(n_col - 10, 0), bound_slope),
        ])[:max_steps]

    # Find first bound crossings — fully vectorized, no Python step loop.
    # crossed_at[i, ss] = True iff |DV[i, ss]| >= bound at step ss.
    # Only checks columns 0..max_steps-1 (matching the original loop range).
    choices_out = np.full(num_trials, np.nan)
    rts_out     = np.full(num_trials, max_steps, dtype=int)
    bounds_out  = np.full(num_trials, np.nan)
    snrs        = trial_means / generative_std
    abs_DV      = np.abs(DV)

    if bound_slope != 0:
        crossed_at = abs_DV[:, :max_steps] >= bound_per_trial[:, None] * slope_vec[None, :]
    else:
        crossed_at = abs_DV[:, :max_steps] >= bound_per_trial[:, None]

    any_crossed = crossed_at.any(axis=1)
    first_ss    = np.argmax(crossed_at, axis=1)  # meaningful only where any_crossed

    # Normal crossings (ss > 0; the common path — requires bound > 0)
    normal = any_crossed & (first_ss > 0)
    if normal.any():
        rows    = np.where(normal)[0]
        ss_vals = first_ss[rows]
        rts_out[rows]     = np.minimum(ss_vals + ndts[rows] + 1, max_steps)
        choices_out[rows] = (DV[rows, ss_vals] > 0).astype(float)
        bounds_out[rows]  = (bound_per_trial[rows] * slope_vec[ss_vals]
                             if bound_slope != 0 else bound_per_trial[rows])

    # ss=0 crossings: only reachable when bound == 0; assign random RT and choice
    ss0 = any_crossed & (first_ss == 0)
    if ss0.any():
        n_ss0 = int(ss0.sum())
        rts_out[ss0]     = rng.integers(1, 3, size=n_ss0)
        choices_out[ss0] = rng.integers(0, 2, size=n_ss0).astype(float)
        bounds_out[ss0]  = (bound_per_trial[ss0] * slope_vec[0]
                            if bound_slope != 0 else bound_per_trial[ss0])

    # Build steps_list: each element is the cumulative DV path for one trial,
    # truncated at RT.  This matches the format of pigeon_steps in the real data
    # and is consumed by get_bounds / infer_bound_mean_per_dt.
    # Pass return_steps=False to skip this allocation (e.g., BiasCorrection.ipynb).
    if return_steps:
        steps_list = [DV[tt, : rts_out[tt]] for tt in range(num_trials)]
    else:
        steps_list = None

    # Lapses: override a random fraction of choices
    if lapse_rate > 0:
        n_lapse   = int(np.ceil(lapse_rate * num_trials))
        lapse_idx = rng.integers(0, num_trials, size=n_lapse)
        choices_out[lapse_idx] = rng.integers(0, 2, size=n_lapse).astype(float)

    return choices_out, rts_out, DV, bounds_out, ndts, snrs, steps_list


def get_simulated_data_table(
    data_table=None,
    specs=None,
    generative_mean=0.05,
    generative_std=0.15,
    num_trials=600,
    max_steps=50,
    bound_mean=0.0,
    bound_std=0.0,
    bound_slope=0.0,
    bound_max=0.75,
    ndt_min=2,
    ndt_max=3,
    lapse_rate=0.0,
    bound_type='given',
    bounds=None,
    correct_bias=BIAS_CORRECTION_FILE,
    rng=None,
):
    """
    Simulate a data table matching the format of get_data_table().
    Translated from getPigeon_simulatedDataTable.m.

    Call modes
    ----------
    data_table mode  (pass data_table)
        Infers per-subject/block generative_mean, NDT, and bound parameters
        from the real data according to bound_type.
    specs mode       (pass specs dict)
        Parametric simulation; all subjects share the same keyword-argument
        parameters.

    Parameters
    ----------
    data_table : pd.DataFrame or None
    specs : dict or None
        Keys: 'num_subjects' (int), 'blocks' (int or list of ints).
    generative_mean : float or array-like
        Mean step size. In specs mode, typically a single positive value
        equal to SNR × generative_std.
    generative_std : float
    num_trials : int
        Maximum trials simulated; actual count is capped at the step budget.
    max_steps : int
        Maximum steps per trial.
    bound_mean : float or array-like
    bound_std : float or array-like
    bound_slope : float
    bound_max : float
    ndt_min, ndt_max : int or array-like
    lapse_rate : float
    bound_type : {'given', 'fixed', 'true', 'var'}
        (data_table mode only)
        'given'  – use the passed bound_mean / bound_std unchanged
        'fixed'  – median |bound| per subject/block; bound_std = 0
        'true'   – sample per-trial bounds from actual data; bound_std = 0
        'var'    – median + std of |bound| per subject/block
        Note: in specs mode, per-SNR bounds can be achieved by passing
        generative_mean and bound_mean as matching-length lists (e.g.
        generative_mean=[lo, hi], bound_mean=[b_lo, b_hi]).  simulate_trials
        then assigns each trial the bound corresponding to its SNR level.
        This is equivalent to MATLAB's boundType='varBySNR'.
    bounds : array-like or None
    correct_bias : str or None
        Path to .mat bias-correction file; applied if the file exists.
    rng : np.random.Generator or None

    Returns
    -------
    pd.DataFrame with columns matching get_data_table() output.
    """
    if rng is None:
        rng = np.random.default_rng()

    # Base simulation kwargs, shared unless overridden per subject/block
    base = dict(
        generative_mean=generative_mean,
        generative_std=generative_std,
        num_trials=num_trials,
        max_steps=max_steps,
        bound_mean=bound_mean,
        bound_std=bound_std,
        bound_slope=bound_slope,
        bound_max=bound_max,
        ndt_min=ndt_min,
        ndt_max=ndt_max,
        lapse_rate=lapse_rate,
        bounds=bounds,
    )

    # ------------------------------------------------------------------
    # Build per-subject × per-block spec dicts
    # ------------------------------------------------------------------
    if data_table is not None:
        subjects = np.sort(data_table['subject_index'].dropna().unique())
        blocks   = np.sort(data_table['block_index'].dropna().unique().astype(int))
        abs_snr  = np.abs(data_table['snr'].to_numpy())
        abs_b    = np.abs(data_table['bound'].to_numpy())
        lg       = data_table['bound'].notna() & (data_table['RT'] >= 0)

        block_specs = [[dict(base) for _ in blocks] for _ in subjects]

        for si, subj in enumerate(subjects):
            ls_sub = data_table['subject_index'] == subj
            for bi, block in enumerate(blocks):
                spec = block_specs[si][bi]
                spec.update(BLOCK_DEFAULTS.get(int(block), BLOCK_DEFAULTS[1]))

                ls_blk   = lg & ls_sub & (data_table['block_index'] == block)
                snr_vals = np.unique(abs_snr[ls_blk & np.isfinite(abs_snr)])
                spec['generative_mean'] = snr_vals * generative_std

                # NDT: derived from congruence on block 1, then copied to later blocks
                if bi == 0:
                    first_idx = ls_sub[ls_sub].index[0]
                    cong = data_table.loc[first_idx, 'congruence']
                    if (cong is not None and isinstance(cong, np.ndarray)
                            and cong.ndim == 2 and cong.shape[1] >= len(snr_vals)):
                        ndt_min_v = np.empty(len(snr_vals), dtype=int)
                        ndt_max_v = np.empty(len(snr_vals), dtype=int)
                        for ri in range(len(snr_vals)):
                            order = np.argsort(cong[:, ri])
                            top2  = order[-2:]
                            # +1 converts 0-indexed Python lag positions to MATLAB's 1-indexed NDT
                            ndt_min_v[ri] = int(top2.min()) + 1
                            ndt_max_v[ri] = int(top2.max()) + 1
                        spec['ndt_min'] = ndt_min_v
                        spec['ndt_max'] = ndt_max_v
                else:
                    spec['ndt_min'] = block_specs[si][0]['ndt_min']
                    spec['ndt_max'] = block_specs[si][0]['ndt_max']

                # Bound parameters according to bound_type
                ab = abs_b[ls_blk]
                ab = ab[np.isfinite(ab)]
                if bound_type == 'fixed':
                    spec['bound_mean'] = float(np.nanmedian(ab))
                    spec['bound_std']  = 0.0
                elif bound_type == 'true':
                    spec['bounds']    = ab
                    spec['bound_std'] = 0.0
                elif bound_type == 'var':
                    spec['bound_mean'] = float(np.nanmedian(ab))
                    spec['bound_std']  = float(np.nanstd(ab))
                # 'given': leave base values unchanged

    elif specs is not None:
        blocks_arr = np.atleast_1d(np.asarray(specs.get('blocks', [1]), dtype=int))
        n_subjects = int(specs.get('num_subjects', 1))
        subjects   = np.arange(1, n_subjects + 1)
        blocks     = blocks_arr

        block_specs = []
        for _ in subjects:
            row = []
            for block in blocks:
                spec = dict(base)
                spec.update(BLOCK_DEFAULTS.get(int(block), BLOCK_DEFAULTS[1]))
                row.append(spec)
            block_specs.append(row)

    else:
        raise ValueError("Provide either data_table or specs.")

    # ------------------------------------------------------------------
    # Run simulations and assemble output rows
    # ------------------------------------------------------------------
    dfs = []

    for si, subj in enumerate(subjects):
        for bi, block in enumerate(blocks):
            spec = block_specs[si][bi]
            bd = {k: spec[k] for k in
                  ('steps_per_block', 'p2p_coins', 'coins_gained', 'coins_lost', 'steps_lost')}
            sim_kw = {k: spec[k] for k in (
                'generative_mean', 'generative_std', 'num_trials', 'max_steps',
                'bound_mean', 'bound_std', 'bound_slope', 'bound_max',
                'ndt_min', 'ndt_max', 'lapse_rate', 'bounds',
            )}

            choices, rts, _DV, _b, _ndts, snrs_sim, steps_sim = simulate_trials(
                **sim_kw, rng=rng,
            )
            n_sim = len(choices)

            # Cumulative steps and coins; truncate when the step budget is reached
            step_counts = np.cumsum(rts + bd['steps_lost'] * (choices == 0))
            coin_counts = np.cumsum(
                -bd['p2p_coins']
                - bd['coins_lost']  * (choices == 0)
                + bd['coins_gained'] * (choices == 1)
            )

            idx         = int(np.searchsorted(step_counts, bd['steps_per_block']))
            trial_count = min(idx + 1, n_sim)

            # Mark last trial as aborted (NaN choice) if budget was exceeded mid-trial
            if trial_count > 0 and step_counts[trial_count - 1] > bd['steps_per_block']:
                choices[trial_count - 1] = np.nan
                if trial_count > 1:
                    coin_counts[trial_count - 1] = coin_counts[trial_count - 2]

            ch = choices[:trial_count]
            sn = snrs_sim[:trial_count]
            st = steps_sim[:trial_count]
            cc = coin_counts[:trial_count]

            # Compute bound and DT from step sequences (same pipeline as real data)
            bound_arr, dt_arr, rt_arr, _, rt_cutoff_arr = get_bounds(
                st, ch, snr=np.abs(sn),
                block_ids=np.full(trial_count, int(block))
            )

            # Build one DataFrame per (subject, block) — faster than appending one dict per trial
            dfs.append(pd.DataFrame({
                'subject_index':  int(subj),
                'block_index':    int(block),
                'trial_number':   np.arange(1, trial_count + 1),
                'bound':          bound_arr,
                'RT':             rt_arr,
                'DT':             dt_arr,
                'choice':         ch,
                'correct':        ch,   # equals choice when generative_mean > 0
                'coin_count':     cc,
                'snr':            sn,
                'steps':          st,
                'rt_cutoff':      rt_cutoff_arr,
                'rt_below_cutoff': rt_arr < rt_cutoff_arr,
            }))

    sim_table = pd.concat(dfs, ignore_index=True) if dfs else pd.DataFrame()

    if len(sim_table) == 0:
        return sim_table

    sim_table['abs_bound'] = np.abs(sim_table['bound'])
    sim_table['wall_hit']  = sim_table['steps'].apply(
        lambda s: (s is not None and len(s) > 0
                   and bool(np.any(np.abs(s[np.isfinite(s)]) >= WALL - 1e-9)))
    )

    if correct_bias and os.path.exists(correct_bias):
        sim_table = apply_bias_corrections(sim_table, correct_bias)

    return sim_table
