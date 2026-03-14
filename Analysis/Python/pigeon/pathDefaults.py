"""Path defaults for the Pigeon Code project.

All project-wide file and directory paths are defined here.  Edit this file
to re-point the project at a new data location without touching any notebook
or analysis routine.
"""

from pathlib import Path

# ── Directory anchors ──────────────────────────────────────────────────────────

_HERE = Path(__file__).resolve().parent   # pigeon/
PROJECT_ROOT = _HERE.parent               # Python/

# ── Raw CSV data (external — pigeon-task recordings per subject) ───────────────
# Top-level directory that contains Pigeon_OL/, Pigeon_MX/, etc.
RAW_DATA_DIR = Path('/Users/jigold/GoldWorks/Mirror_jigold/Manuscripts/2023_Pigeon/Data')

# ── Processed / derived data (stored with the project) ────────────────────────
DATA_DIR = PROJECT_ROOT.parent.parent.parent / 'Data' / 'Processed' / 'Python'

BIAS_CORRECTION_FILE = DATA_DIR / 'boundBiasCorrection.mat'   # bound-bias regression coefficients
BOUND_SUMMARIES_FILE = DATA_DIR / 'boundSummaries_wCP.mat'    # CP-task bound summaries (MATLAB)
RR_MATS_FILE         = DATA_DIR / 'RRMats.mat'                # CP reward-rate matrices
RR_MATS_OL_FILE      = DATA_DIR / 'RRMats_OL.mat'             # OL reward-rate matrix

# ── Figures output ─────────────────────────────────────────────────────────────
# Shared Figures/ folder three levels above the Python project root
# (i.e. at the 2025_Pigeon project level)
FIGURES_DIR = PROJECT_ROOT.parent.parent.parent / 'Figures'
