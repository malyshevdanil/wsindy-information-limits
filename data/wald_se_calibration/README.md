# Paired 5% null standard-error diagnostic

This bundle preserves all five completed `secal_0`–`secal_4` output sets,
including 300 method rows, diagnostics for all 100 observed paths, every
one of the 39,900 bootstrap targets, completion records and execution logs.
`manifest.json` records byte counts and SHA256 hashes. Scientific source
snapshots and configurations are embedded in the corresponding frozen
notebooks under `kaggle/`.

Reproduce the audited full-sample analysis from the repository root:

```sh
python scripts/finish_wald_se_calibration.py --analyze-only
```

This command requires no API token, Kaggle connection or ODE simulation.
It checks the bundle, source snapshots, unchanged protocol and reference
hashes, all original targets and nominal SEs, every bootstrap SE, and all
three rejection decisions on every path.

The paired null outcomes are nominal 23/100, random-design sandwich 2/100,
and full-refit parametric bootstrap 5/100. Both revised SEs pass the frozen
Wilson-compatibility criterion. These are previously examined null paths;
independent size validation and power with the revised SEs remain pending.
