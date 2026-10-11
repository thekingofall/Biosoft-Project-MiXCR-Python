# Usage guide
[English](usage.md) · [简体中文](usage.zh-CN.md) · [Home](../README.md)

## Install
Python 3.10 or later and NumPy are required. From the repository root, create and activate a virtual environment, then run `python -m pip install ./python`. Check the installation with `mixcrpy --version` and `mixcrpy --help`. `python -m mixcrpy` is equivalent and selects the current interpreter explicitly.

## Single and paired reads
```bash
mixcrpy analyze reads.fastq.gz output/single --chains TRB
mixcrpy analyze reads_R1.fastq.gz output/paired --read2 reads_R2.fastq.gz --chains BCR
```

Chain choices: ALL, TCR, BCR, TRA, TRB, TRG, TRD, IGH, IGK, IGL. The default is ALL. FASTQ records require nonempty names, nucleotide sequences and matching Phred+33 quality scores.

Paired files must have equal record counts and matching identifiers. CASAVA mate descriptions and /1, /2 suffixes are accepted. If your files deliberately pair by record order and names differ, use `--pair-by-order`; establish their ordering first. Original descriptions contribute to deterministic wildcard handling. `--no-merge` retains separate targets. The supported paired layout is Opposite.

## Separate alignment and assembly
```bash
mixcrpy align reads.fastq output/aligned
mixcrpy assemble output/aligned.alignments.json output/clones
```

The default assembly feature is CDR3. Alignment JSON records its chain and V-reference profile for later steps. `--v-feature VGeneWithP` selects the alternate V-gene feature; VRegionWithP is the default.

## Multiple features and gene separation
```bash
mixcrpy analyze full_reads.fastq output/cdr123 --features CDR1 CDR2 CDR3
mixcrpy analyze full_reads.fastq output/vdj --features VDJRegion --separate-by V J C
```

These features require sufficient sequence coverage. `--separate-by` makes the chosen gene assignments part of clone separation. The tested profiles include CDR3, CDR1+CDR2+CDR3, CDR3+FR4 and VDJRegion.

## Partial-read workflow
```bash
mixcrpy align partial.fastq output/raw --allow-partial-alignments
mixcrpy assemble-partial output/raw.alignments.json output/overlaps
mixcrpy assemble output/overlaps.alignments.json output/clones
mixcrpy assemble-contigs output/overlaps.alignments.json output/clones.clones.json output/contigs
```

Paired data can add `--read2 partial_R2.fastq` to alignment. `--allow-no-cdr3-part-alignments` retains qualifying partial alignments without a CDR3-adjacent anchor. `assemble-partial --overlapped-only` emits accepted overlaps. The resulting alignment document preserves contributing read IDs; use this same document for both assembly and contigs.

## Contigs
```bash
mixcrpy assemble-contigs output/sample.alignments.json output/sample.clones.json output/contigs
```

Use matching alignment and clone documents from the same run and reference profile. Contigs require CDR3 clones and use the default FullSeqAssembler configuration. Outputs are `.contigs.fasta`, `.contigs.tsv`, `.contigs.json` and `.report.json`. Separate branches and targets may exist within a clone.

## Compare with the original
The optional developer script requires your own authorized copy of the original MiXCR 3.0.13 JAR and Java. The repository contains captured expectations for its portable tests; it does not bundle the executable.

```bash
python scripts/compare_example.py --jar /path/to/mixcr-3.0.13.jar --fastq examples/trb.fastq --output results/original-check
```

Use `--java /path/to/java` if Java is outside PATH. The script checks the original JAR SHA-256, runs alignment, assembly and exports for both programs, and saves commands, logs and the comparison. Choose a new output directory. Its comparison covers the fields declared in the script; the historical broader validation is described in [Validation](validation.md).

## Reproduce the Python checks
```bash
python -m pip install pytest
python -m pytest python/tests -q
```

The tests use included JSON expectations generated from the original program. Java is unnecessary for these checks. To regenerate the synthetic examples, run `python scripts/make_example.py`.

## Practical limits
Input errors stop the command. Existing outputs require `--force` to replace. Use one output prefix per run. A workflow interrupted during writing can leave partial files; rerun with a fresh prefix.

The CLI loads alignment records into memory. Start with a representative subset and assess memory and run time before a large analysis. Timing equivalence and clinical performance have not been established. Reference data and parameters affect results: compare against 3.0.13 using the same profile. See [scope and validation](validation.md) before interpreting differences.

Native `.vdjca`, `.clns` and `.clna` files require the original tool; this package reads its own JSON intermediates. Latest MiXCR features and license terms are summarized [separately](latest-version.md).
