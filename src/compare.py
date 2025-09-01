from Bio import SeqIO

# 2bit encoding
base_to_bits = {'A': '00', 'C': '01', 'G': '10', 'T': '11'}

def encode_2bit(seq):
    return ''.join([base_to_bits.get(base.upper(), '00') for base in seq])

# 1. Load reference sequence
ref_seq = None
for record in SeqIO.parse("chr11.fasta", "fasta"):
    ref_seq = str(record.seq)
    break

# 2. Read variant data from binary file
variant_dict = {}
with open("chr11_combined_binary.hex", "rb") as f:
    ref_len = int.from_bytes(f.read(4), 'big')
    variant_count = int.from_bytes(f.read(4), 'big')
    f.read(ref_len)  # Skip reference sequence

    for _ in range(variant_count):
        pos = int.from_bytes(f.read(4), 'big')
        diff_bytes = f.read(2)
        variant_dict[pos] = diff_bytes

# 3. Compare with VCF for validation
mismatch_count = 0
checked = 0

with open("HG00157.chr11.vcf", "r") as vcf:
    for line in vcf:
        if line.startswith("#"):
            continue
        tokens = line.strip().split("\t")
        pos = int(tokens[1])
        ref = tokens[3]
        alt = tokens[4]
        gt = tokens[9].split(":")[0]

        # Check if the variant is valid
        if alt == '.' or '1' not in gt:
            continue

        # Verify XOR difference based on position
        start = pos - 1
        ref_4mer = ref_seq[start:start+4]
        alt_4mer = list(ref_4mer)
        alt_4mer[0] = alt[0]  # 간단화된 방식
        alt_4mer = ''.join(alt_4mer)

        ref_bin = encode_2bit(ref_4mer)
        alt_bin = encode_2bit(alt_4mer)
        expected_xor = int(ref_bin, 2) ^ int(alt_bin, 2)
        expected_bytes = expected_xor.to_bytes(2, 'big')

        if pos in variant_dict:
            if variant_dict[pos] != expected_bytes:
                print(f"Mismatch at {pos}: expected {expected_bytes.hex()}, got {variant_dict[pos].hex()}")
                mismatch_count += 1
        else:
            print(f"Missing variant at {pos}")
            mismatch_count += 1
        
        checked += 1
        if checked >= 100:  # 100개만 확인
            break

print(f"\n검사 완료: {checked}개 중 {mismatch_count}개 mismatch")
