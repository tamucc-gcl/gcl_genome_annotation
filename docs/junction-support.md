# RNA splice-junction QC

`main.nf` now runs report-only Portcullis QC (`--run_junction_qc true`) alongside
the other postprocessing work. Each readset is analyzed separately, with its
resolved library strandedness. Native BAMs, CSI indices and the masked reference
are passed directly through Nextflow. No remapping, transcript assembly or BAM
merging is added. The existing BRAKER, UTR and AED inputs remain unchanged.

Portcullis 1.2.4 performs its own preparation, junction analysis and classification.
No predicted annotation is supplied through its reference-rescue option. The
native classifier threshold defaults to 0.5. Low-junction datasets can trigger
Portcullis's own lenient rule-based fallback; this is recorded in the manifest
and the native filter log is retained. The pipeline does not replace its scoring
or invent a second junction filter. `--extra` and BAM filtering are not enabled.

Readset analysis uses 8 CPUs, 48 GB and 48 hours by default, configured centrally
with `portcullis_cpus`, `portcullis_memory`, `portcullis_time` and
`portcullis_threshold`. Physical mate orientation is left UNKNOWN rather than
assuming all supplied paired libraries are FR. Reverse/forward library strand
maps to firststrand/secondstrand; unstranded libraries use unstranded mode.

## Outputs

Under `qc/junctions/<assembly>/`:

- `README.md`: interpretation and report links.
- `transcript_junction_support.tsv`: separate all/coding/UTR-associated counts
  and passing fractions, one row per accepted transcript.
- `junction_support.tsv`: exact model introns, passing/failed readsets, unresolved
  strand observations and raw spliced-alignment counts. Introns use **1-based
  inclusive** coordinates. Native Portcullis TAB coordinates are converted from
  zero-based inclusive coordinates.
- `junction_status.json`: final model/reference identities, readset provenance,
  native table checksums, settings and summary counts. Summary counts refer to
  junction occurrences in transcripts; shared introns can appear more than once.
- `readsets/<readset>/`: native all/pass/fail tables, BED/GFF exports when produced,
  software version, input manifest and preparation/analysis/filter logs. Prepared
  BAM/reference copies are not published.

Exact donor/acceptor boundaries and strand must match. A model junction passes
when at least one readset passes it. A readset rejecting it remains visible even
if another passes. A junction observed on the matching strand but passing nowhere
is `fail`; one absent from those observations is `not_observed`. Observations
whose strand cannot be resolved are marked separately and do not count as
positive strand-specific support. Opposite-strand observations do not support
the model. Counts are alignments, not deduplicated fragments or independent RNA
molecules. Multiple readsets from one library are not independent replicates.

A coding junction has CDS reaching both adjacent exon boundaries. All other
model junctions are labeled `utr_associated`, including mixed CDS/UTR boundaries.
Fractions use only the corresponding category as denominator. Single-exon
transcripts have no junction fraction. RNA-free assemblies report `no_rna` with
unavailable support values rather than a zero support score.

These are complementary confidence measurements, not gene filters or AED
components. Lack of RNA support can reflect tissue/expression sampling. Passing
individual junctions does not validate full transcript chains or UTR assignment
between neighboring genes. Input guards verify accepted genome sequence identity,
BAM dictionary/sort/index and final annotation checksums. Native table integrity
and reference bounds are checked before reporting. Nextflow owns file staging
and caching. Report/model changes do not invalidate native junction analysis;
changing BAMs, strand metadata or classifier settings does.

## Cluster acceptance and real-data pilot

First run the tiny actual-tool test on Crest:

```bash
python gcl_genome_annotation/tests/validate_junctions_cluster.py \
  --outdir testing/junction_contracts_01
```

It checks native Portcullis on positive/negative-strand junctions, multiple
readsets, coding/UTR separation, RNA-free and single-exon cases, main syntax,
resume, report-only/model-only cache reuse and disabled QC. It does not run
mapping, repeats, BRAKER or InGenAnnot. Recovery uses the same output directory
with `--resume`. Tiny native fallback acceptance does not validate full-scale
random-forest behavior; the next step is the real-data pilot.

After acceptance, the usual `main.nf -resume` command includes this QC. To pilot
only junction QC against existing accepted outputs, substitute `junction_qc.nf`
for `main.nf` and add `--junction_inputs_root "${outdir}"`; retain the same root
configuration, profile and work directory. This entry reads published final
annotation and prediction metadata, while integrated main uses native channels.
Moving between the two entries may change task cache identities. Keep the same
entry for subsequent cache tests. `--run_junction_qc false` disables the work.

Implementation and Python tests are verified locally. Actual Nextflow/Portcullis
cluster acceptance and full-data pilot are pending.

## Sources

- [Portcullis paper](https://pmc.ncbi.nlm.nih.gov/articles/PMC6302956/)
- [Native workflow and options](https://portcullis.readthedocs.io/en/latest/using.html)
- [Metrics and pooling behavior](https://portcullis.readthedocs.io/en/latest/metrics.html)
- [Bioconda package recipe](https://github.com/bioconda/bioconda-recipes/blob/master/recipes/portcullis/meta.yaml)
