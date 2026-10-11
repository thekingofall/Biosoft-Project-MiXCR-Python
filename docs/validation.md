# Validation — Biosoft-Project-MiXCR-Python 0.2.0

[中文](validation.zh-CN.md) · [Machine-readable summary](validation-summary.json)

The recorded validation observed zero mismatches in the measured fields across six FASTQ workflows against MiXCR 3.0.13. This is finite differential evidence for the listed inputs, parameters and output fields. It does not establish equivalence for every possible input or every MiXCR feature.

## Reference and environment

The reference executable identifies itself as MiXCR 3.0.13, commit `f61422609eb4743e1e55ad3fb12ab51e34952b62`. Its SHA-256 is `1526a6d36221e86b775e39831b5ebaab6d5877639212c672bef2049a3fd37149`. Source review used the [official v3.0.13 tree](https://github.com/milaboratory/mixcr/tree/v3.0.13), fixed at `045e9afdb83f99a7c39b2e75807736c4632aa344`. The recorded comparison between these revisions found changes in CHANGELOG and POM; the algorithm sources matched. The exported human reference comes from repseqio.v1.6.

The historical run was recorded on 2026-10-12 with Windows 11, Python 3.12.7, NumPy 2.4.6 and Java 21. Python runtime calculations use Python and NumPy. A separate audit completed the 200-pair IGH workflow while enforcing a prohibition on child-process creation and network access; no blocked operations were attempted. Java served as the development-time reference and oracle generator.

## FASTQ-to-clone comparison

Each row starts from FASTQ and includes alignment, clone assembly, clone export and read-to-clone assignment. Paired input counts represent read pairs. The original program ran with matched parameters and a single thread.

| Dataset | Input records | Aligned records | Clones | Alignment / clone / assignment checks | Observed mismatches |
| --- | ---: | ---: | ---: | ---: | ---: |
| Synthetic TRB, single-end | 1,039 reads | 1,036 | 4 | 7,252 / 36 / 2,075 | 0 |
| Synthetic TRB, overlapping paired-end | 1,036 pairs | 1,036 | 4 | 7,252 / 36 / 2,072 | 0 |
| Official TRB, single-end | 100 reads | 100 | 22 | 700 / 198 / 122 | 0 |
| Official IGH, single-end | 200 reads | 104 | 45 | 728 / 405 / 279 | 0 |
| Official TRB, paired-end, VGeneWithP | 100 pairs | 99 | 22 | 693 / 198 / 122 | 0 |
| Official IGH, paired-end | 200 pairs | 169 | 86 | 1,183 / 774 / 337 | 0 |

“Checks” are the comparator's field-check counters. They are not independent biological samples. Synthetic reads reuse four rearrangement templates and include reverse complements, V-region substitutions and indels, high- and low-quality CDR3 mutations, a CDR3 insertion and unmatched reads. The single-end synthetic set contains 17 distinct sequences. Repeated records exercise abundance, mapping and clustering.

Measured alignment fields include alignment success/failure, all retained V/D/J/C gene hits and scores, CDR3 nucleotide sequence, translation and quality. Measured clone fields include these annotations, count and fraction. Assignment checks cover clone membership and mapping type. Python comparisons normalize equally scoring gene names as sets and match baseline clones by CDR3. Thus these baseline checks do not establish identical ordering for every tied hit or uniquely resolve every same-CDR3 clone under arbitrary gene-separation settings. Dedicated extension cases separately check the latter workflow. The recorded Java rebuild also matched the original exports row for row in all six datasets.

Counts follow the original default assembly rule: a clustered read can receive a clone assignment without increasing the retained clone's count. Assigned-read totals and clone-count sums therefore have distinct meanings.

## Additional measured workflows

Seven assembly cases cover CDR1+CDR2+CDR3, CDR3+FR4, VDJRegion, separation by V or V/J/C, and disabled preclustering. All compared sequences, qualities, gene scores, counts, fractions and assignments agree. Three CLI cases include FASTQ analysis, JSON reload and V separation. The V-separation case preserves ten preclustered read labels. Feature extraction additionally agrees for 1,859 paired-IGH feature groups and 506 full-VDJ feature groups, including N-region coordinates. Directed cases cover zero-length junctions and the reference program's different target-selection rules for sequence quality and relative coordinates.

Partial-read assembly has ten single-target and ten paired-target oracle cases. Three chained FASTQ workflows were recorded: a three-read example produces two clones with all source read identifiers retained; the 200-pair IGH workflow produces 86 clones with 951 checks passing; and the three-read partial-to-contig workflow produces two matching contigs. Weighted read counts and per-read assignment expansion are checked after merging.

Full contig command comparisons observed zero measured-field mismatches for synthetic TRB (4 clones), official IGH single-end (45 clones), and official IGH paired-end (86 clones). The paired result contains 89 target sequences and a clone-count sum of 137. Checks cover target sequence and quality, CDR3, counts, all retained V/D/J/C names and scores, and each target alignment's score, reference/query interval and complete M/X/I/D operation path. Contig fixtures also cover 240 point-variant/clustering/trimming cases, 24 read-to-point-to-contig cases and 12 selected command-output clones.

Other historical original-Java checks include 9,000 random/hash/wildcard checks, 8,493 assembly mapping-formula checks, 500 clonal V/J banded alignments, 2,800 multi-segment compatibility cases, 240 seed-mapper cases, 776 gene-extension alignments and 516 paired-merge cases. These counts describe different test layers; they should not be added as independent end-to-end validation datasets.

## Reproduce the published regression suite

From the repository root, install the package and test dependencies, then run:

```sh
python -m pip install ./python pytest pytest-subtests
python -m pytest python/tests
```

The [published tests and saved JSON oracle fixtures](../python/tests/) can be run without Java. Fixtures supply expected outputs to tests; runtime algorithms do not read expected-output fixtures. The retained historical suite log reports **52 tests and 2,214 subtests passed**. This count is a dated result; a fresh checkout's test run determines its current result.

The full historical FASTQ runs, original executable, Java rebuild, raw exports and private command logs are outside this public release. The table and JSON summary are an audited transcription of those results. Public fixtures support regression of the packaged cases; the historical full-run tables cannot be regenerated from this checkout alone. No unpublished evidence path is presented as a downloadable artifact.

## Evidence schema and interpretation

[validation-summary.json](validation-summary.json) records the project/version, reference identity and executable hash, historical environment, six dataset counters, extension and contig results, supported profile, limitations, and hashes of the published test fixtures. All file references in that document are relative to the repository root. `evidence_status` separates historical local-run summaries from packaged regression inputs. A `passed` result means zero observed mismatches in the named measured fields for that case; it is not a probability estimate of correctness.

## Compatibility limits

The validated profile covers the exported human reference, default linear scoring, V-then-J, Max quality aggregation, Opposite paired layout, and the documented feature and gene-separation paths. Complete CDR3 extraction requires contiguous coverage or a successful merge. Contig command validation covers CDR3 clones and default contig parameters.

KAligner2, other species, other paired layouts, arbitrary scoring models, composite features, alternate contig assembling features and custom subCloningRegion require further implementation or verification. Deferred-read candidate ties, exhaustive indel-placement paths in repetitive sequences, and exceptional contig-retention fallbacks remain incompletely sampled. Reports expose known boundaries encountered during execution. Python uses its own JSON interchange format; native `.vdjca`, `.clns` and `.clna` file I/O is outside this Python release. Performance, all-platform behavior and clinical suitability have not been established by these tests.

The port retains the historical MiXCR/MiLib research-license restrictions; see [LICENSE.mixcr](../python/LICENSE.mixcr) and the [dependency notice](../python/mixcrpy/data/NOTICE). Validation does not change the applicable license terms.
