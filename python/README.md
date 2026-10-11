# MiXCR Python reconstruction

[English documentation](https://github.com/thekingofall/Biosoft-Project-MiXCR-Python) · [简体中文](https://github.com/thekingofall/Biosoft-Project-MiXCR-Python/blob/main/README.zh-CN.md)

Python 3.10+ / NumPy implementation of selected human V(D)J alignment, clone assembly, partial-read assembly and contig reconstruction workflows. Compatibility target: MiXCR 3.0.13. Package version: 0.2.0.

Install with `python -m pip install .` from this directory, then run `mixcrpy --help`. From the repository root, use `python -m pip install ./python`.

```bash
mixcrpy analyze reads.fastq output/sample
mixcrpy analyze reads_R1.fastq.gz output/paired --read2 reads_R2.fastq.gz
```

See the repository's bilingual usage guide and validation report for supported profiles, checked datasets, example commands and known limits. Python uses JSON intermediates, TSV tables and FASTA contigs.

MiXCR/MiLib-derived portions retain the research license in LICENSE.mixcr. It permits specified educational, research and nonprofit uses by nonprofit institutions; commercial use requires separate permission. Apache Commons Math-derived code retains its Apache 2.0 license and notices in mixcrpy/data.
