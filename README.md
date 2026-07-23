# obsview_python

A modular pipeline for reading observation diagnostic files and producing
4-panel statistics plots. Reads **ODS** (`.ods`) and **IODA/JEDI** (`.nc4`)
files and plots, per channel: observation count (used vs. unused), mean & RMS
of O-B / O-A, Jo/p, and prescribed vs. estimated errors.

## Running it

Everything currently runs from `run.py` (or `test.py`) — there is no CLI yet.
Set the input file and options directly in the script, then:

```bash
python run.py
```

All CLI flags are unchanged from the original script (`--scale`,
`--obtype`, `--var`, `--tarname`, `--satid`, `--nbins`, `--qc`, `--xGSI`,
`--common`, `--bias`, `--fig`). Run `python run.py --help` for the full
list.


## Pipeline

The code is organized as a linear pipeline. Each stage takes and returns a
plain data object, so the readers (ODS/IODA) are interchangeable and every
downstream stage is format-agnostic:

read → mask → filter → derive → bin → stats → plot

1. **Read** — `ODSReader` / `IODAReader` load a file into an `ObservationData`
   object (raw arrays plus metadata: full channel list `all_lev`, per-variable
   `fill_values`).
2. **Mask & filter** — build a boolean mask (`fill_val_mask`, `qc_pass_mask`,
   `qc_fail_mask`) and apply it with `apply_filter`. Metadata like `all_lev`
   is preserved across filtering.
3. **Derive** — `calc_derived` computes `job`, `joa`, `esigo`, `esigb` *after*
   filtering, so fill values never enter the math (avoids overflow).
4. **Bin** — `create_channel_bins` sorts data by channel and builds a
   `BinnedData` object on the full channel axis.
5. **Stats** — `calculate_stats` returns a `StatisticsData` object
   (per-channel counts, means, RMS, Jo/p, errors).
6. **Plot** — `plot_stats` draws the 4-panel figure and returns the matplotlib
   `Figure` so you can `plt.show()` or `fig.savefig(...)` it.

## Core data objects

- `ObservationData` — per-observation arrays + metadata (`all_lev`,
  `fill_values`).
- `BinnedData` — data sorted by bin, plus `bin_centers`, `bin_indices`,
  `bin_labels`, `bin_heights`.
- `StatisticsData` — per-channel aggregated statistics.

## Latest changes

- Added an **`IODAReader`** class to read `.nc4` files into the same
  `ObservationData` structure used by `ODSReader`, so the rest of the pipeline
  is unchanged between formats.
- Added **QC-fail observations** to the observation-count panel: QC-pass (used)
  and QC-fail (unused) counts share the same channel y-tick, drawn as
  overlapping green-over-red bars. Failed-QC data is counted for the plot only
  and never enters the statistics.

