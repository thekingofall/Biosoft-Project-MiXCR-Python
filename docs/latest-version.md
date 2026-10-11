# Official MiXCR v4.7.0 and this Python port

[中文](latest-version.zh-CN.md)

MiXCR v4.7.0 is the release marked **Latest** on the official GitHub release page when checked on 12 October 2026. It was published on 7 August 2024. This page describes that fixed version; future releases may change the latest-version designation. [Official release](https://github.com/milaboratory/mixcr/releases/tag/v4.7.0).

The official release emphasizes consistent assembly features across samples, improved single-cell processing and gene assignment, and corrections to downstream analyses.

| Area | Official v4.7.0 changes |
| --- | --- |
| Assembly | Explicit feature selection where the protocol leaves it unspecified; separate hit-filter settings for clone and contig assembly. |
| Single-cell and barcodes | Cell isolation fixes, per-chain consensus assembly, optional tag-correction bypass and MiTool preassembly in the 10x VDJ preset. |
| Downstream analysis | SHM-tree topology and export fixes; `groupClones` renamed to `assembleCells`. |
| References and inputs | Rabbit/sheep IG references, human reference corrections, long-read and BAM fixes, and additional kit presets. |

These are upstream release-note claims. Full option names, migration details and the complete change list are available in the [v4.7.0 notes](https://github.com/milaboratory/mixcr/releases/tag/v4.7.0).

**The Python implementation targets MiXCR 3.0.13 behavior.** Its research scope includes human V(D)J alignment, overlapping and unmerged paired reads, feature extraction, clone assembly, partial-read assembly and default CDR3 contig reconstruction. Its tests compare scientific output with the original 3.0.13 program. Supported profiles and data formats are described in the project README.

Single-cell barcode correction, `assembleCells`, allele inference, SHM-tree construction, MiTool integration, newer kit presets and v4.7.0 reference-library updates are outside this port's supported interface. Results validated against 3.0.13 establish compatibility for those tested workflows. They do not establish equivalence to v4.7.0. The Python package uses its own JSON intermediates, TSV tables and FASTA contigs; native MiXCR binary-format compatibility remains outside its interface.

The v4.7.0 archive is available from the [official download](https://github.com/milaboratory/mixcr/releases/download/v4.7.0/mixcr-4.7.0.zip). The downloaded artifact was checked with the following values; the hash is a project-computed checksum, not a vendor signature.

| Artifact | Verified value |
| --- | --- |
| File | `mixcr-4.7.0.zip` |
| Size | 106,402,425 bytes |
| SHA-256 | `4fde83750f1d1308ba6b7431bcafe98297b361da7e17efc4d072ed935045cac7` |

The [v4.7.0 license](https://github.com/milaboratory/mixcr/blob/v4.7.0/LICENSE) defines eligible academic use and grants conditional download, installation and use rights. It includes restrictions on modification, redistribution, derivatives, reverse engineering and disclosure of performance analyses. Its definition of academic use also includes funding conditions; academic affiliation alone does not describe every requirement. Consult the complete text and the licensor for use outside its grant.

This repository's research Python implementation derives from the older 3.0.13 sources and retains their applicable notices and terms. Its distribution contains no v4.7.0 executable or implementation. Access to the official v4.7.0 download remains governed by that release's license. MiXCR is a MiLaboratories product; this reconstruction has no official endorsement.
