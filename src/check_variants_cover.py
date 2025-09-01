import os
import struct

def encode_base2bit(base):
    return {'A': 0b00, 'C': 0b01, 'G': 0b10, 'T': 0b11}.get(base, 0b00)

def encode_seq_variable(seq):
    bits = 0
    for base in seq:
        bits <<= 2
        bits |= encode_base2bit(base)
    byte_len = (len(seq) * 2 + 7) // 8
    return bits.to_bytes(byte_len, 'big')

# Step 1: extract every variants from vcf
def extract_vt_variants(vcf_path):
    variants = []
    base_dir = os.path.dirname(os.path.abspath(__file__))
    full_vcf_path = os.path.join(base_dir, vcf_path)

    with open(full_vcf_path, 'r') as f:
        for line in f:
            if line.startswith('#'):
                continue

            fields = line.strip().split('\t')
            if len(fields) < 8:
                continue

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

# Step 2: convert to binary record
def convert_variants_to_binary_records(vt_variants):
    binary_records = []
    for pos, (vt_type, ref, alt) in vt_variants:
        if vt_type == 'SNP' and len(ref) == 1 and len(alt) == 1:
            record = (
                b'\x00' +
                struct.pack('>I', pos) +
                bytes([encode_base2bit(ref)]) +
                bytes([encode_base2bit(alt)])
            )
            binary_records.append((pos, record))

        elif vt_type == 'INDEL':
            if len(ref) > len(alt) and ref.startswith(alt):  # Deletion
                del_seq = ref[len(alt):]
                if len(del_seq) == 0:
                    continue
                encoded = encode_seq_variable(del_seq)
                record = (
                    b'\x02' +
                    struct.pack('>I', pos) +
                    struct.pack('>H', len(del_seq)) +
                    encoded
                )
                binary_records.append((pos, record))

            elif len(ref) < len(alt) and alt.startswith(ref):  # Insertion
                ins_seq = alt[len(ref):]
                if len(ins_seq) == 0:
                    continue
                encoded = encode_seq_variable(ins_seq)
                record = (
                    b'\x01' +
                    struct.pack('>I', pos) +
                    struct.pack('>H', len(ins_seq)) +
                    encoded
                )
                binary_records.append((pos, record))

    print(" Muatation from first 15 binary record (hex):")
    for _, rec in binary_records[:15]:
        print(rec.hex())

    return binary_records

# Step 3: save as binary 
def write_binary_records(output_path, binary_records):
    with open(output_path, 'wb') as f:
        for _, record in binary_records:
            f.write(record)
    print(f" saved as binary file: {output_path}")

# Step 4: Extract POS from binary file
def extract_positions_from_bin(bin_path):
    positions = set()
    with open(bin_path, 'rb') as f:
        while True:
            header = f.read(1)
            if not header:
                break
            variant_type = header[0]
            pos = struct.unpack('>I', f.read(4))[0]
            positions.add(pos)
            if variant_type == 0x00:
                f.read(2)  # REF+ALT
            elif variant_type in (0x01, 0x02):  # INS or DEL
                length = struct.unpack('>H', f.read(2))[0]
                byte_len = (length * 2 + 7) // 8
                f.read(byte_len)
            else:
                break
    return positions


if __name__ == '__main__':
    vcf_file = "HG00157.chr11.vcf"
    bin_file = "variants_extracted.bin"

    vt_variants = extract_vt_variants(vcf_file)
    binary_records = convert_variants_to_binary_records(vt_variants)
    write_binary_records(bin_file, binary_records)

   
    vcf_pos = set(pos for pos, _ in vt_variants)
    bin_pos = extract_positions_from_bin(bin_file)

    common = vcf_pos & bin_pos
    vcf_only = sorted(vcf_pos - bin_pos)
    bin_only = sorted(bin_pos - vcf_pos)

    print(f"\n total mutation: {len(binary_records)}")
    print(f" shared mutation: {len(common)}")
    print(f" POS only in VCF: {len(vcf_only)}")
    print(f" POS only in bin: {len(bin_only)}")

    print("\n first 10 - POS only in VCF :", vcf_only[:10])
    print("first 10 -  POS only in Bin:", bin_only[:10])
