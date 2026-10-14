# Complete endpoint and added-term study outputs

All four endpoint shards (1968 paired records, four branches) and all 1000 complete author WSINDy fits passed the frozen source/configuration audit. The rows here preserve every coefficient result, including zero inclusions.

Frozen protocols are `configs/boundary_confirmation.json` and `configs/support_extra_term.json`; execution sources are embedded in the five frozen notebooks under `kaggle/boundary_confirm_*` and `kaggle/support_extra`. Observed time series can be regenerated from those protocols. Full fit artifacts remain in the downloaded Kaggle output ZIPs; this compact archive contains the numerical results and provenance needed for analysis.

To rerun local analysis, copy each stage directory to `results/kaggle_<stage>/`, then run `scripts/analyze_endpoint_confirmation.py`, `scripts/report_endpoint_confirmation.py` and `scripts/analyze_support_extra_term.py`. Those scripts verify the archived notebook source and require complete samples. Research simulations run on Kaggle.
