import os

def calculate_compression_ratio(vcf_path, fasta_path, compressed_path):
    vcf_size = os.path.getsize(vcf_path)
    fasta_size = os.path.getsize(fasta_path)
    compressed_size = os.path.getsize(compressed_path)

    original_size = vcf_size + fasta_size

    compression_percent = (1 - compressed_size / original_size) * 100
    compression_ratio = original_size / compressed_size

    print(" Results")
    print("────────────────────────────")
    print(f" VCF size:        {vcf_size / (1024 * 1024):.2f} MB")
    print(f" FASTA size:      {fasta_size / (1024 * 1024):.2f} MB")
    print(f" original size:    {original_size / (1024 * 1024):.2f} MB")
    print(f" compressed file size: {compressed_size / (1024 * 1024):.2f} MB")
    print("────────────────────────────")
    print(f" Compression percent:          {compression_percent:.2f}% 감소")
    print(f" compression_ratio:       {compression_ratio:.2f}배")


if __name__ == '__main__':
    base_dir = os.path.dirname(os.path.abspath(__file__))

    vcf_file = os.path.join(base_dir, 'HG00157.chr11.vcf')
    fasta_file = os.path.join(base_dir, 'chr11.fasta')
    compressed_file = os.path.join(base_dir, 'combined_final.hex')

    calculate_compression_ratio(vcf_file, fasta_file, compressed_file)
