"""
Pigeon — Python analysis package for the Pigeon task.

Modules
-------
data
    Load and pre-process trial data from raw CSV files.  Computes per-trial
    bounds, decision times (DT), response times (RT), and applies optional
    bias correction.

stats
    Summary statistics over the processed data table — per-subject/block/SNR
    accuracy, DT distributions, and bound statistics binned by DT.

simulate
    Generative model of the task.  Simulates a random-walk accumulator with
    an absorbing bound, NDT, and optional lapse rate.  Can reproduce a full
    data-table in the same format as `get_data_table`.

Public API
----------
get_data_table          Load real participant data as a DataFrame.
get_bounds              Infer per-trial bound and DT from step sequences.
get_good_trial_array    Boolean mask for standard trial-quality filters.
get_bound_summary       Mean |bound| per subject/block/SNR/DT bin.
get_performance_summary Accuracy and DT stats per subject/block/SNR.
simulate_trials         Core simulation engine (returns raw arrays).
get_simulated_data_table Full simulated data table matching get_data_table format.
"""

from .data import get_data_table, get_bounds, get_good_trial_array
from .stats import get_bound_summary, get_performance_summary
from .simulate import simulate_trials, get_simulated_data_table
