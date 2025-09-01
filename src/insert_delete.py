
import struct
import os

base_dir = os.path.dirname(__file__)
vcf_path = os.path.join(base_dir, "HG00157.chr11.vcf")
out_path = os.path.join(base_dir, "chr11_indels.bin")


def encode_base2bit(base):
    return {'A': 0b00, 'C': 0b01, 'G': 0b10, 'T': 0b11}[base]

def encode_sequence_fixed16(seq):
    bits = 0
    for base in seq:
        bits <<= 2
        bits |= encode_base2bit(base)
    byte_len = (len(seq) * 2 + 7) // 8
    return bits.to_bytes(byte_len, 'big').ljust(16, b'\x00')  # 고정 16바이트

def extract_indels_using_vt(vcf_path, output_bin):
    with open(vcf_path, 'r') as vcf, open(output_bin, 'wb') as out:
        for line in vcf:
            if line.startswith('#'):
                continue
            parts = line.strip().split('\t')
            pos = int(parts[1])
            ref = parts[3]
            alt = parts[4]
            info = parts[7]

            if 'VT=INDEL' not in info:
                continue  # SNP, SV 등은 제외

            for alt_allele in alt.split(','):
                
                # Deletion
                if len(ref) > len(alt_allele):
                    del_len = len(ref) - len(alt_allele)
                    out.write(b'\x00')  # type = deletion
                    out.write(struct.pack('>I', pos))
                    out.write(struct.pack('B', del_len))
                # Insertion
                elif len(ref) < len(alt_allele):
                    insert_seq = alt_allele[len(ref):]
                    if any(b not in 'ACGT' for b in insert_seq):
                        print(f"[SKIP] Invalid base in insertion at {pos}: {insert_seq}")
                        continue
                    if len(insert_seq) > 64:
                        print(f"[TRUNCATE] Insertion too long at {pos}, cutting to 64 bases: {insert_seq[:64]}")
                        insert_seq = insert_seq[:64]
                    encoded = encode_sequence_fixed16(insert_seq)
                    out.write(b'\x01')  # type = insertion
                    out.write(struct.pack('>I', pos))
                    out.write(struct.pack('B', len(insert_seq)))
                    out.write(encoded)  # always 16 bytes
                

extract_indels_using_vt(vcf_path, out_path)