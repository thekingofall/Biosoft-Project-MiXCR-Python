# Workflow illustration

`workflow.png` is the shared English-language figure used by both project READMEs. `workflow.pdf` embeds the same image for download. Captions and the surrounding explanation follow each README's language.

## Visual specification

Use a white background, large Arial-style labels, dark outlines and restrained panel borders. The palette follows the project's reference illustration: jade teal `#00877C`, slate blue `#477EA9`, muted amethyst `#8F79B4`, amber ochre `#D69A32`, and coral accent `#C86E64`. V is blue, D is ochre, J is teal, and C is lavender throughout the diagram. These are target prompt colors; raster shading introduces additional shades.

The style reference is the [fetal immune atlas graphical abstract](https://github.com/thekingofall/scFetallife_Blood/blob/main/docs/figures/graphical_abstract.png). The artwork here depicts this software's own workflow and uses an independently generated composition.

## Scientific interpretation

The main lane shows FASTQ input, alignment, clone assembly and inspectable outputs. Partial-read overlap and contig reconstruction are optional. D segments apply to the appropriate receptor loci. Nucleotide contigs use A/C/G/T; the CDR3 annotations in the clone schematic use amino-acid symbols. Gene identities and sequences are illustrative. The counts 12 and 8 are the verified result of the bundled 20-read synthetic example.

The lower lane describes an independent comparison with original MiXCR 3.0.13. The Python runtime uses Python and NumPy. Validation applies to the finite cases in the [validation report](../validation.md).

## Image provenance

Created on 12 October 2026 with Codex's built-in image-generation tool, followed by a targeted correction of the nucleotide contig lettering. The tool did not expose a model identifier. The final PNG was visually checked against the documented workflow; the PDF is a format export of that image. The generation specification is saved in [workflow-prompt.md](workflow-prompt.md).
