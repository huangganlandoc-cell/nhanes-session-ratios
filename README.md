# Time of blood collection and leukocyte ratio indices: analysis code

Analysis code for the study *Time of blood collection and leukocyte ratio indices: evidence from randomized examination
sessions in the National Health and Nutrition Examination Survey*.

The analysis plan was registered on the Open Science Framework before the outcome data were unblinded
(https://osf.io/x62ra; version 1.0, 26 September 2026). The registered protocol is included in `analysis/prereg/`
together with its SHA-256 checksums.

Archived at Zenodo: https://doi.org/10.5281/zenodo.22993589 (this DOI resolves to the latest version).

## What the code does

- **NHANES III (1988–1994).** Reconstructs the random assignment of households to a morning or an afternoon/evening
  examination session from the session-specific examination weights (WTPFSD6, WTPFMD6) and estimates the
  intention-to-treat effect and the complier average causal effect of session on leukocyte counts and ratio indices
  (granulocyte-to-lymphocyte ratio as the three-part-differential analogue of the neutrophil-to-lymphocyte ratio;
  granulocyte-based systemic immune-inflammation index; platelet-to-lymphocyte ratio). Variance accounts for the
  survey design (Taylor linearization; Fay balanced repeated replication weights).
- **Continuous NHANES 1999–2018 and 2021–2023.** Replication by the session attended (five-part differential: NLR,
  SII, SIRI, PIV, PLR); proportions above cutoffs by marginal standardization with delete-one-PSU jackknife intervals.
- **Mortality** (public-use linked mortality files, follow-up through 31 December 2019). Relative change in the log
  hazard ratio after session standardization, with the session effect re-estimated within each replicate.
- **Prespecified exploratory analyses**: subgroups, restricted cubic spline turning points, specification curve.
- **Independent re-implementation** of the key NHANES III estimates in Python (`analysis/py_check/`), written from the
  protocol without access to the R code, followed by an item-by-item comparison with the R results.

## Requirements

Tested on macOS (Darwin 25.5, Apple silicon) with

- R 4.6.1 with survey 4.5, survival 3.8.6, data.table 1.18.6.1 and jsonlite 2.0.0 (splines and parallel are part of
  base R). The R scripts use forked parallel processing (`parallel::mclapply`) and therefore need macOS or Linux.
  `gzip` must be available on the command line.
- Python 3.14 with the packages in `requirements.txt`.

## Data

Only public NHANES files are used (about 540 MB). `download_data.sh` downloads them from the National Center for
Health Statistics into `$NHANES_DATA_DIR` (default `~/data/nhanes_public`):

- NHANES III laboratory, examination and adult questionnaire files with their SAS layout files
  (https://wwwn.cdc.gov/nchs/data/nhanes3/1a/);
- NHANES III and continuous NHANES 1999–2018 public-use linked mortality files, 2019 release
  (https://ftp.cdc.gov/pub/Health_Statistics/NCHS/datalinkage/linked_mortality/);
- the 95 continuous NHANES files listed in `data_files_continuous.txt` (https://wwwn.cdc.gov/Nchs/Data/Nhanes/Public/).

The files used for the reported results were downloaded on 26 September 2026; `SHA256SUMS_raw_data.txt` lists their
checksums. NCHS occasionally revises public files, so a checksum mismatch indicates a different file version.
No individual-level data are included in this package.

## How to run

```sh
export NHANES_DATA_DIR=~/data/nhanes_public
sh download_data.sh
cd analysis/code
python3 01_build_nhanes3.py        # NHANES III analytic file
python3 02_build_continuous.py     # continuous NHANES analytic file
python3 03_data_notes.py           # sample flow, balance and adherence (no outcome analyses)
python3 04_export_and_scramble.py  # CSV files for R, plus scrambled copies for blinded development
cd ../R
MODE=real sh run.sh                # all R analyses; log in analysis/results/real/run_log.txt
cd ../py_check
python3 key_estimates.py --input ../derived/nhanes3_REAL.csv.gz
cd ../code
python3 05_compare_R_python.py real
python3 06_write_summary.py real   # applies the prespecified decision rules and wording (protocol section 7)
python3 07_figures.py
python3 08_tables.py
```

The R analyses take about 4 minutes on a 12-core computer. Without `MODE=real`, the R pipeline and the Python check run
on the scrambled data (session labels and mortality outcomes permuted), which was used to develop the code before
registration. Runs on the real data require `analysis/prereg/OSF_REGISTRATION.txt` to contain the registration link,
which is included.

## Reproducing the reported numbers

`reference_results/` contains the outputs used for the manuscript: `results_REAL.json` (all estimates),
`key_estimates_REAL.json` (Python check), `compare_R_python_REAL.txt`, `spec_curve.csv` and `summary_REAL.md`.
The analyses are deterministic (bootstrap replicate weights are drawn once with seed 20260926), so a fresh run
reproduces these values exactly; only the recorded run times differ. We verified this by running the package from a clean copy on 2026-09-26: 1,631 numeric values in `results_REAL.json` and the Python check were compared with `reference_results/`, and all were identical.

## Notes

- Code comments, console messages and `summary_REAL.md` are partly in Chinese.
- Figures use the Arial font when it is available.

## License

MIT (see `LICENSE`).
