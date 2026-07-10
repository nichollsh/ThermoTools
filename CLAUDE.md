# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## What this is

ThermoTools downloads, parses, and compiles thermodynamic data (heat capacities, saturation
pressures, latent heats, equations of state) for chemical species into per-species NetCDF files,
for use in astrophysics/planetary-science modelling. It is a data pipeline, not an application:
most of the "logic" lives in the `workflow/` Jupyter notebooks, which call into the `thermotools`
Python package for shared functions and reference data.

## Setup

```bash
pip install -Ue .
```

Requires Python >=3.12. Dependencies: matplotlib, netCDF4, numpy, scipy, requests (see `pyproject.toml`).

There is no test suite, linter, or CI configured in this repo. There is no `run` command — the
package's only interface is the notebooks plus library functions imported by them.

## Running the pipeline

The `workflow/` notebooks must be run in numeric order — each stage reads inputs (raw downloads
or previous stage's output) and writes to `generated/`, which later stages and the final
compile step depend on:

1. `01_cp_download.ipynb` — download heat capacity data from JANAF tables (per-element, per-species; slow)
2. `02_cp_parse.ipynb` — parse downloaded heat capacity data
3. `03_lv_parse.ipynb` — parse latent heat data
4. `04_psat_parse.ipynb` — parse saturation curve data
5. `05_vdw_tabulate.ipynb` — tabulate density rho(P,T) via the Van der Waals EOS
6. `06_aqua_tabulate.ipynb` — download and tabulate the AQUA water EOS
7. `07_chabrier_tabulate.ipynb` — download and tabulate the CMS19 hydrogen EOS
8. `99_compile.ipynb` — compile everything into one NetCDF (.nc) file per species

`generated/*/` is gitignored (rebuilt by running the notebooks), except checked-in subfolders
like `generated/sat/dat/` and `generated/sat/plt/` (saturation curve data/plots) and `generated/cp/web/`
(downloaded JANAF text files), which are committed as reference/cached inputs.

## Architecture

- `src/thermotools/__init__.py` — defines the project's root-relative path helpers used
  everywhere else: `get_gendir()` (→ `generated/`) and `get_inpdir()` (→ `src/thermotools/data/`),
  plus `empty_dir()`. All path resolution in the package goes through these two functions rather
  than hardcoded paths, so moving the repo root doesn't break anything.
- `src/thermotools/data/` — checked-in reference/source data, organized by topic
  (`mmw/` molecular weights, `lv/` latent heat/vapour data, `psat/` saturation pressure fits),
  each with a `web/` (raw downloaded/source) and sometimes `dat/` (parsed) subfolder. This is the
  input side; `generated/` is the output side, with the same `web/`/`dat/` split convention.
- `src/thermotools/moles.py` — molecular weight lookups: parses the NIST elements table and
  computes molar mass from a chemical formula string (`count_atoms`, `mmw_from_formula`).
- `src/thermotools/phase.py` — saturation pressure (Psat) models. Two independent sources of
  truth, selected per-gas via the `"prefer"` key in the `gases` dict: Antoine equation fits
  (`antoine()`, coefficients hardcoded per gas/temperature range) and Clausius-Clapeyron
  integration from tabulated latent heat data (`cc_psat()`, reads from `generated/lv/dat/`).
  Separately, `woitke_load_all()` parses Woitke+2017 Table D2 at import time (module-level side
  effect — printed gas list on import) and `woitke_psat()` dispatches to one of ~10 fit-type
  equations (`match gas_data["fit_type"]`) drawn from various literature sources. Gases already
  covered by the Antoine/CC `gases` dict are skipped when loading Woitke data to avoid duplicates.
- `src/thermotools/water.py` — water-specific IAPWS95/Wagner & Pruß (2002) correlations
  (`psat_liquid`, `psat_solid`, `psat_both`, `rho_liq`, `rho_vap`, `delta_vap`) plus a tabulated
  sublimation enthalpy curve (`delta_sub`, PCHIP-interpolated) for temperatures below the triple
  point. Kept separate from `phase.py` because water needs its own dedicated EOS rather than a
  simple fit.
- `src/thermotools/util.py` — generic download/archive helpers (`download`, `untar`, `unzip`,
  `makezip`, `checksum`/`writesum` using BLAKE2b) used by the notebooks to fetch and cache
  external datasets (JANAF, AQUA, CMS19).
- `src/thermotools/plot.py` — shared matplotlib styling and a `pretty_name()` helper that
  LaTeX-ifies chemical formulas (digits → subscripts) for plot labels.
- `misc/` — standalone exploratory notebooks (not part of the main pipeline).

## Working with gas/species data

When adding a new species' saturation pressure fit, decide which mechanism it belongs in:
- Add to the `gases` dict in `phase.py` if you have literature Antoine-equation coefficients or
  want Clausius-Clapeyron integration from latent-heat data — set `"prefer"` to pick which.
- Otherwise it may already be covered by Woitke+2017 Table D2 (`woitke_load_all()` loads anything
  not already in `gases`); check `woitke_gases` before adding a manual fit.

Fit coefficients and citations are hardcoded inline (not in a data file) in `phase.py` — always
cite the literature source in the `"cite"` list when adding a new gas.
