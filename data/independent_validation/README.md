# Independent validation outputs

Original-byte outputs from fully completed frozen studies only. The manifest
lists completed stages and records whether the full validation is complete.
A completed subset is not presented as the completed 69-job study.

Each stage retains its full protocol, executed-source hash, completion record,
coefficient rows and diagnostics. Frozen notebook sources are under kaggle/iq_*.
The acceptance bands, observation seeds and bootstrap seeds are fixed in
configs/independent_validation_v2.json. No partial rejection rates are reported.
The analysis script verifies all prescribed seeds and archived-source hashes.
