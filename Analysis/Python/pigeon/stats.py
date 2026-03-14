"""
Summary statistics utilities for the Pigeon project.
Python translations of getPigeon_performanceSummary, getPigeon_boundSummary.

Both functions apply ``get_good_trial_array`` internally, which enforces the
standard quality filters:

- ``DT ≥ 1``             — excludes trials too short to estimate a bound
- ``trial_number ≥ 10``  — excludes early warm-up trials per block
- ``wall_hit = False``   — excludes trials where the pigeon reached the screen
                            edge (|position| ≥ 0.75), which inflate RT/DT

Input
-----
Both functions expect a DataFrame produced by ``pigeon.data.get_data_table``
(or ``get_simulated_data_table``).  Required columns:
``subject_index``, ``block_index``, ``snr``, ``DT``, ``correct``,
``trial_number``, ``wall_hit``, ``bound``.

Output shapes
-------------
``get_performance_summary`` → ``(n_subjects, n_blocks, n_snrs, 4)``
    Last dim: [accuracy, mean_DT, std_DT, n_trials]

``get_bound_summary`` → ``(n_subjects, n_blocks, n_snrs, n_rt_bins, 3)``
    Last dim: [mean_abs_bound, std_abs_bound, n_trials]
    ``n_rt_bins`` depends on the ``max_rt`` argument:
      - integer    → one bin per integer DT value (1 … max_rt)
      - 'medsplit' → two bins split at the median total RT (not DT); trials with
                     RT ≤ median go into bin 0, RT > median into bin 1
      - other      → single bin (all DTs collapsed)
"""

import numpy as np
import pandas as pd
from .data import get_good_trial_array


def get_performance_summary(data_table, blocks="all", snrs="all"):
    """
    Compute per-subject/block/SNR accuracy and DT stats.
    Translated from getPigeon_performanceSummary.m

    Returns
    -------
    np.ndarray, shape (n_subjects, n_blocks, n_snrs, 4)
        Last dim: [accuracy, mean_DT, std_DT, n_trials]
    """
    subject_indices = np.sort(data_table["subject_index"].dropna().unique())
    n_subjects = len(subject_indices)

    block_indices = np.arange(1, 4) if blocks == "all" else np.atleast_1d(np.asarray(blocks))
    n_blocks = len(block_indices)

    abs_snr = np.abs(data_table["snr"].to_numpy())
    if snrs == "all":
        abs_snrs = np.sort(np.unique(abs_snr[np.isfinite(abs_snr)]))
    else:
        abs_snrs = np.asarray(snrs)
    n_snrs = len(abs_snrs)

    lg = get_good_trial_array(data_table)

    summary = np.zeros((n_subjects, n_blocks, n_snrs, 4))

    for si, subj in enumerate(subject_indices):
        ls = lg & (data_table["subject_index"] == subj)

        for bi, block in enumerate(block_indices):
            lsb = ls & (data_table["block_index"] == block)

            for ni, snr_val in enumerate(abs_snrs):
                lsnr = lsb & np.isclose(abs_snr, snr_val)
                if lsnr.any():
                    subset = data_table[lsnr]
                    summary[si, bi, ni] = [
                        subset["correct"].sum() / len(subset),
                        subset["DT"].mean(),
                        subset["DT"].std(),
                        len(subset),
                    ]
                else:
                    summary[si, bi, ni] = [np.nan, np.nan, np.nan, 0]

    return summary


def get_bound_summary(data_table, correct_only=False, blocks="all",
                      split_by_snr=True, max_rt=10):
    """
    Compute per-subject/block/SNR/RT-bin bound statistics.
    Translated from getPigeon_boundSummary.m

    Returns
    -------
    np.ndarray, shape (n_subjects, n_blocks, n_snrs, n_rt_bins, 3)
        Last dim: [mean_abs_bound, std_abs_bound, n_trials]
    """
    subject_indices = np.sort(data_table["subject_index"].dropna().unique())
    n_subjects = len(subject_indices)
    abs_bounds = np.abs(data_table["bound"].to_numpy())

    block_indices = (
        np.sort(data_table["block_index"].dropna().unique().astype(int))
        if blocks == "all"
        else np.atleast_1d(np.asarray(blocks))
    )
    n_blocks = len(block_indices)

    if split_by_snr:
        abs_snr = np.abs(data_table["snr"].to_numpy())
        snrs = np.sort(np.unique(abs_snr[np.isfinite(abs_snr)]))
        n_snrs = len(snrs)
        snr_masks = np.column_stack([np.isclose(abs_snr, s) for s in snrs])
    else:
        n_snrs = 1
        snr_masks = np.ones((len(data_table), 1), dtype=bool)

    if isinstance(max_rt, int):
        rt_bins = np.arange(1, max_rt + 1)
    elif max_rt == "medsplit":
        rt_bins = np.array([1, 2])
    else:
        rt_bins = np.array([0])  # sentinel for "all"
    n_rts = len(rt_bins)

    summary = np.full((n_subjects, n_blocks, n_snrs, n_rts, 3), np.nan)

    lg = get_good_trial_array(data_table)   # applies min_dt, min_trial_number, wall_hit
    if correct_only:
        lg = lg & (data_table["correct"] == 1)
    else:
        lg = lg & (data_table["correct"] >= 0)
    lg = lg & data_table["bound"].notna() & (data_table["bound"] != 0)

    for si, subj in enumerate(subject_indices):
        lsubj = lg & (data_table["subject_index"] == subj)

        for bi, block in enumerate(block_indices):
            lsb = lsubj & (data_table["block_index"] == block)

            for ni in range(n_snrs):
                lsnr = lsb & snr_masks[:, ni]

                if isinstance(max_rt, int):
                    for ri, rt_val in enumerate(rt_bins):
                        lrt = lsnr & (data_table["DT"] == rt_val)
                        if lrt.any():
                            summary[si, bi, ni, ri] = [
                                np.nanmean(abs_bounds[lrt]),
                                np.nanstd(abs_bounds[lrt]),
                                lrt.sum(),
                            ]
                        else:
                            summary[si, bi, ni, ri] = [np.nan, np.nan, 0]

                elif max_rt == "medsplit":
                    med_rt = np.nanmedian(data_table.loc[lsnr, "RT"])
                    lrt_masks = [
                        lsnr & (data_table["RT"] <= med_rt),
                        lsnr & (data_table["RT"] > med_rt),
                    ]
                    for ri, lrt in enumerate(lrt_masks):
                        if lrt.any():
                            summary[si, bi, ni, ri] = [
                                np.nanmean(abs_bounds[lrt]),
                                np.nanstd(abs_bounds[lrt]),
                                lrt.sum(),
                            ]

                else:
                    if lsnr.any():
                        summary[si, bi, ni, 0] = [
                            np.nanmean(abs_bounds[lsnr]),
                            np.nanstd(abs_bounds[lsnr]),
                            lsnr.sum(),
                        ]

    return summary
