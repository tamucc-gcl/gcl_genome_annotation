# Regression tests and development diagnostics

See [the testing guide](../docs/testing.md) for current commands, native-tool
scope, focused entries and recovery behavior. Run the Python suite from the
pipeline checkout; run cluster controllers from the consuming project with
its normal Nextflow/Conda/Singularity environment.

Current `test_*.py`, `check_*.groovy`, fixture workflows and stage validators
protect active behavior. Passing their original milestone does not make them
disposable. In particular, several validators import shared controller helpers
from `validate_utr_cluster.py` and `validate_preflight_cluster.py`; later unit
tests also import fixtures from earlier test modules.

The completed custom-AED comparison and one-off review utilities were retired
with their exclusive tests. Shared model/native-score fixtures are now in
`annotation_fixtures.py`; the UTR boundary regression helper is
`utr_label_checks.py`. Current native AED/model checks remain in
`test_ingenannot_aed.py`, `test_annotation_models.py` and `test_final_annotation.py`.
Do not delete an entire regression file to remove only its legacy cases.
See [the consolidation record](../docs/cleanup.md).

Native cluster acceptance uses synthetic inputs and is not proof of full
real-data training in every evidence mode. See [validation limits](../docs/validation.md).

Focused stage launchers live in [validation/entrypoints](../validation/entrypoints).
The validators use those paths and preserve their layout in staged fixtures.
`check_validation_entries.groovy` checks root configuration resolution, includes
and staged helper/environment paths for production and focused launchers.
