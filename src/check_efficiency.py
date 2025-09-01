import os

base_dir = "/Users/jayjung/Comp571/final project"
vcf = os.path.join(base_dir, "HG00157.chr11.vcf")
fasta = os.path.join(base_dir, "chr11.fasta")
hexfile = os.path.join(base_dir,  "chr11_fasta_with_ref_N_masking.hex")

vcf_size = os.path.getsize(vcf)
fasta_size = os.path.getsize(fasta)
hex_size = os.path.getsize(hexfile)

original_size = vcf_size + fasta_size
compression_ratio = 1 - (hex_size / original_size)
compression_percent = round(compression_ratio * 100, 2)

print(f"VCF+FASTA: {original_size / (1024*1024):.2f} MB")
print(f"Compressed HEX: {hex_size / (1024*1024):.2f} MB")
print(f"압축률: {compression_percent}% reduced")