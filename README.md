# 2-bit-encoded-VCF-format

This repository implements a self-contained binary format for genomic variation data. Unlike conventional VCF + FASTA storage, this format integrates the reference genome, variant records, and metadata into a single file, making it compact, portable, and fully restorable.

Motivation
The Variant Call Format (VCF) is human-readable but inefficient for large-scale analysis:
Storing VCF together with a FASTA reference wastes space.
VCF alone cannot reconstruct a full genome.
Our 2-bit encoded format solves these issues by:
Embedding the entire reference sequence inside the compressed file.
Encoding SNPs, insertions, and deletions with fixed-length binary records.
Preserving ambiguous bases with a 1-bit mask.
Attaching compact metadata (allele frequency, depth, genotype).
This achieves ~87% compression compared to VCF + FASTA while supporting full genome restoration
.
File Layout
Each compressed file (*.hex) follows this structure:
Block	Description
Reference Block	Full reference genome, 2-bit encoded (A=00, C=01, G=10, T=11).
Mask Block	1-bit per base, marking ambiguous sites (N).
Variant Block	Fixed-length records for SNPs, insertions, deletions.
Metadata Block	Compact per-variant fields: AF, DP, GT.
Variant Encoding Structure
Each variant begins with a 1-byte type flag:
0x00: SNP (Single Nucleotide Polymorphism)
0x01: Insertion
0x02: Deletion

SNP (7 bytes)
1 byte: Type flag (0x00)
4 bytes: Position (POS)
1 byte: Reference base (REF)
1 byte: Alternate base (ALT)

Insertion (23 bytes)
1 byte: Type flag (0x01)
4 bytes: Position (POS)
2 bytes: Length of inserted sequence (LEN)
16 bytes: Inserted sequence (ALT), encoded in 2-bit format

Deletion (23 bytes)
1 byte: Type flag (0x02)
4 bytes: Position (POS)
2 bytes: Length of deleted sequence (LEN)
16 bytes: Deleted sequence (REF), encoded in 2-bit format

Metadata Encoding
Each variant is mapped one-to-one with metadata fields, appended in the same order:
1 byte: Allele Frequency (AF), scaled to 0–255
1 byte: Read Depth (DP), capped at 255
1 byte: Genotype (GT):
0 = Homozygous reference (0/0)
1 = Heterozygous (0/1)
2 = Missing / ambiguous
3 = Homozygous alternate (1/1)
The metadata block is preceded by a 4-byte META marker and a 4-byte length indicator.

Workflow
All steps are wrapped into compression.py, but conceptually the pipeline is:

Reference + Mask Encoding
Input: FASTA
Output: reference_with_mask.hex

Variant Encoding
Input: VCF
Output: variants_extracted.bin

Metadata Encoding
Input: VCF INFO / FORMAT fields
Output: metadata.bin

Final Assembly
Input: Reference + Mask + Variants + Metadata
Output: complete_version_with_ref_N.hex

Restoration
Input: Final .hex file
Output: Reconstructed FASTA sequence


Results (HG00157, Chr11)
Compression
VCF + FASTA: 573.85 MB
2-bit encoded format: 74.50 MB
→ 87.0% reduction
Variant Inclusion
Variants in VCF: 3,881,791
Encoded: 3,870,547
Genome Restoration
Alignment (10 Mbp region): 96.3% sequence similarity
Remaining 3.7% differences reflect true biological variation, not encoding errors


Future Work
Support longer insertions and structural variants.
Extend to multi-sample encoding.
Improve compatibility with standard pipelines (VCF/CRAM).
Explore integration with alternative references (e.g. T2T pangenome).
.
