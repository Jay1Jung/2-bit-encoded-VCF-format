import os
import struct

def encode_base2bit(base):
    return {'A': 0b00, 'C': 0b01, 'G': 0b10, 'T': 0b11}.get(base, 0b00)

def encode_seq_variable(seq):
    bits = 0
    for base in seq:
        bits <<= 2
        bits |= encode_base2bit(base)
    byte_len = (len(seq) * 2 + 7) // 8  # 2 bit per base
    return bits.to_bytes(byte_len, 'big')

# Step 1: Extract all VCF mutation without filtering GT
def extract_vt_variants(vcf_path):
    variants = []
    base_dir = os.path.dirname(os.path.abspath(__file__))
    full_vcf_path = os.path.join(base_dir, vcf_path)

    with open(full_vcf_path, 'r') as f:
        for line in f:
            if line.startswith('#'):
                continue

            fields = line.strip().split('\t')
            pos = int(fields[1])
            ref = fields[3]
            alt_field = fields[4]
            info = fields[7]

            if 'VT=SNP' in info:
                for alt in alt_field.split(','):
                    variants.append((pos, ('SNP', ref, alt)))
            elif 'VT=INDEL' in info:
                for alt in alt_field.split(','):
                    variants.append((pos, ('INDEL', ref, alt)))
    return variants

# Step 2: Convert mutation into binary record 
def convert_variants_to_binary_records(vt_variants):
    binary_records = []

    for pos, (vt_type, ref, alt) in vt_variants:
        if vt_type == 'SNP' and len(ref) == 1 and len(alt) == 1:
            record = (
                b'\x00' +                             # SNP 타입
                struct.pack('>I', pos) +              # POS (4B)
                bytes([encode_base2bit(ref)]) +       # REF (1B)
                bytes([encode_base2bit(alt)])         # ALT (1B)
            )
            binary_records.append(record)

        elif vt_type == 'INDEL':
            if len(ref) > len(alt) and ref.startswith(alt):  # Deletion
                del_seq = ref[len(alt):]
                encoded = encode_seq_variable(del_seq)
                record = (
                    b'\x02' +
                    struct.pack('>I', pos) +
                    struct.pack('>H', len(del_seq)) +
                    encoded
                )
                binary_records.append(record)

            elif len(ref) < len(alt) and alt.startswith(ref):  # Insertion
                ins_seq = alt[len(ref):]
                encoded = encode_seq_variable(ins_seq)
                record = (
                    b'\x01' +
                    struct.pack('>I', pos) +
                    struct.pack('>H', len(ins_seq)) +
                    encoded
                )
                binary_records.append(record)

    print("🧾 앞 15개의 변이 binary record (hex):")
    for rec in binary_records[:15]:
        print(rec.hex())

    return binary_records

# Step 3: Store as binary file
def write_binary_records(output_path, binary_records):
    with open(output_path, 'wb') as f:
        for record in binary_records:
            f.write(record)
    print(f"Binary file stored to: {output_path}")

if __name__ == '__main__':
    vcf_file = "HG00157.chr11.vcf"
    output_bin_file = "variants_extracted.bin"

    vt_variants = extract_vt_variants(vcf_file)
    binary_records = convert_variants_to_binary_records(vt_variants)
    write_binary_records(output_bin_file, binary_records)

    print(f" Number of finally encoded binary: {len(binary_records)}")
