# obsview_python

Refactored, modular version of `ioda_prs_binned.py`. Reads ODS (`.ods`)
and IODA/JEDI (`.nc4`) diagnostic files and produces 4-panel statistics
plots (obs count, o-b/o-a mean & RMS, Jo/p, prescribed vs. estimated
error), either per pressure bin or per channel.

## Running it

From the directory *containing* `obsview/` (not from inside it):

```bash
python run.py path/to/file.ods --obtype mls55_aura
python run.py path/to/file.nc4 --obtype radiance --xGSI --fig out.png
```

All CLI flags are unchanged from the original script (`--scale`,
`--obtype`, `--var`, `--tarname`, `--satid`, `--nbins`, `--qc`, `--xGSI`,
`--common`, `--bias`, `--fig`). Run `python run.py --help` for the full
list.


## Structure

```
obsview/
├── main.py              CLI entry point: parses args, dispatches to the
│                         right stats function, shows/saves the figure
├── cli.py                argparse setup (--obtype, --scale, --nbins, ...)
├── config.py              obtype -> varname/kt/levlim lookup tables
├── utils.py                 safe_filled() and other small generic helpers
│
├── io/
│   └── readers.py        ioda_from_tarball() (nc4, optionally from a
│                          tarball), is_ods()/file_extension() dispatch
│
├── stats/
│   ├── binning.py         shared bin-construction + mask/accumulate
│   │                      helpers (make_pressure_bins, make_channel_bins,
│   │                      get_pressure_mask, accum_pressure_mask)
│   ├── ods_stats.py        ods_pressure_binned, ods_channel
│   ├── obsview.py        ioda_pressure_binned, ioda_channel
│   └── compare_stats.py      jediXgsi_channel, jediXgsi_pressure_binned
│
├── plotting/
│   └── panels.py          show_nobs_panel, show_resstats_panel,
│                           show_jo_panel, show_sigo_panel, comp_nobs_panel
│
└── coverage/, timeavg/, report/   empty packages, staged for the
    planned features below
```

Each `stats/*.py` function still both *computes* the statistics and
*calls* the plotting panels directly (same as the original script) — I
kept that coupling for this pass to minimize the risk of introducing
bugs during the split. See "Suggested next step" below for why
decoupling those two responsibilities is worth doing before you build
the next few features.

## What changed vs. the original script

- Pure reorganization — function bodies are unchanged except for a
  few no-op cleanups (removing dead/unused variables, replacing a
  count-via-loop with `np.sum(bin_mask)`, that kind of thing).
- `main()`'s big `if/elif` chain mapping `--obtype` to `varname`/`kt` is
  now the `OBTYPE_TO_VARNAME` / `VARNAME_TO_KT` dicts in `config.py`.
  Adding a new obtype is now a one-line dict entry instead of edits in
  three places.
- Verified equivalence: for `ods_pressure_binned`, `ods_channel`,
  `ioda_pressure_binned`, `ioda_channel`, and `jediXgsi_channel`, I
  generated synthetic input files and diffed the refactored output
  against the original script pixel-for-pixel (`np.array_equal` on the
  rendered PNGs) — all identical.

## Known pre-existing issues (not introduced by this refactor, found while verifying it)

These reproduce identically in the original `ioda_prs_binned.py`, confirmed
by running both side by side against the same synthetic data:

1. **`ods_channel` crashes if zero observations are QC-excluded.**
   `count_data.unstack()` only produces an `excluded_count` column when
   at least one row has `qcexcl != 0`; if every observation passes QC,
   `count_full['excluded_count']` raises `KeyError`.
2. **`jediXgsi_pressure_binned` can crash on mismatched valid-obs counts.**
   It computes `bin_indices` once from the JEDI QC mask, then reuses
   those same `bin_indices` for both the JEDI and GSI calls into
   `accum_pressure_mask`. If GSI and JEDI don't have the same set of
   valid observations (their QC masks differ, which is the normal
   case), `accum_pressure_mask`'s internal `valid_mask` for GSI selects
   a different number of observations than `bin_indices` has entries
   for, and indexing fails.
3. **`accum_pressure_mask` references an undefined `scale` variable**
   when `scaleby != 'null'` (it was never passed in as a parameter) —
   only triggers if you pass `--scale` together with `--xGSI`.

None of these block your existing typical usage (default `--scale null`,
mixed QC data), which is presumably why they haven't surfaced yet — but
worth knowing about since you'll likely touch these functions again soon.

## Suggested next step

Right now every `stats/*.py` function ends by calling `plt.figure()` and
the `show_*_panel` functions directly — computing and plotting are still
one function. For the features you mentioned (multi-file input,
time-averaged stats, HTML reports), it's worth introducing a shared
result object, e.g.:

```python
@dataclass
class BinnedStats:
    varname: str
    bin_centers: np.ndarray
    bin_heights: np.ndarray
    sum_nobs: np.ndarray
    mean_ombg: np.ndarray
    rms_ombg: np.ndarray
    # ... etc, one field per stat currently computed
    radiance: bool = False
    metadata: dict | None = None
```

Then each `stats/*.py` function returns a `BinnedStats` instead of
plotting directly, and a single function in `plotting/panels.py` takes a
`BinnedStats` and draws the 4-panel figure. That one change is what
would let you:
- loop over multiple files and get back a `list[BinnedStats]`
- average a `list[BinnedStats]` (one per time) into a time-averaged one
- serialize a `BinnedStats` to a table/JSON for the HTML report, with
  no matplotlib involved at all
- add a coverage-plot module that reuses the same masked/valid arrays
  without duplicating the file-reading logic


