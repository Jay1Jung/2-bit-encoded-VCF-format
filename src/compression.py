import os
import struct
from Bio import SeqIO

def encode_base2bit(base):
    if base not in 'ACGT':
        return 0b00  # fallback to A (used only when masked)
    return {'A': 0b00, 'C': 0b01, 'G': 0b10, 'T': 0b11}[base]

def encode_seq_fixed16(seq):
    bits = 0
    for base in seq:
        bits <<= 2
        bits |= encode_base2bit(base)
    byte_len = (len(seq) * 2 + 7) // 8
    return bits.to_bytes(byte_len, 'big').ljust(16, b'\x00')

def encode_af(val):
    try:
        return int(round(min(max(float(val), 0.0), 1.0) * 255))
    except:
        return 0

def encode_dp(val):
    try:
        return min(int(val), 255)
    except:
        return 0

def encode_gt(gt_str):
    gt = gt_str.replace('|', '/').split('/')
    if gt == ['0', '0']:
        return 0
    elif '0' in gt and '1' in gt:
        return 1
    elif gt == ['1', '1']:
        return 3
    return 0

def encode_reference_with_mask(fasta_path):
    record = next(SeqIO.parse(fasta_path, "fasta"))
    seq = str(record.seq).upper()
    bitstream = []
    mask_bits = []

    bits = 0
    count = 0
    for base in seq:
        if base not in 'ACGT':
            mask_bits.append(1)
            base = 'A'
        else:
            mask_bits.append(0)

        bits = (bits << 2) | encode_base2bit(base)
        count += 1
        if count % 4 == 0:
            bitstream.append(bits & 0xFF)
            bits = 0

    if count % 4 != 0:
        bits <<= (2 * (4 - (count % 4)))
        bitstream.append(bits & 0xFF)

    mask_bytes = bytearray()
    for i in range(0, len(mask_bits), 8):
        chunk = mask_bits[i:i+8]
        byte = 0
        for bit in chunk:
            byte = (byte << 1) | bit
        byte <<= (8 - len(chunk))  # 패딩
        mask_bytes.append(byte)

    return bytes(bitstream), mask_bytes

def generate_ref_hex_with_mask(vcf_filename, fasta_filename, output_hex_filename):
    base_dir = os.path.dirname(os.path.abspath(__file__))
    vcf_path = os.path.join(base_dir, vcf_filename)
    fasta_path = os.path.join(base_dir, fasta_filename)
    variant_bin = os.path.join(base_dir, "ref_variants.bin")
    meta_bin = os.path.join(base_dir, "meta_variants.bin")

    variant_count = 0

    with open(vcf_path, 'r') as vcf, \
         open(variant_bin, 'wb') as out_var, \
         open(meta_bin, 'wb') as out_meta:

        for line in vcf:
            if line.startswith('#'):
                continue

            parts = line.strip().split('\t')
            pos = int(parts[1])
            ref = parts[3]
            alt_field = parts[4]
            info = parts[7]
            format_tags = parts[8].split(':')
            sample_values = parts[9].split(':')

            info_dict = dict(i.split('=') for i in info.split(';') if '=' in i)
            af = encode_af(info_dict.get('AF', 0))
            dp = encode_dp(info_dict.get('DP', 0))
            gt_str = sample_values[format_tags.index('GT')] if 'GT' in format_tags else '0/0'
            gt_indices = [int(g) for g in gt_str.replace('|', '/').split('/') if g.isdigit()]

            alt_list = alt_field.split(',')
            if all(alt == ref for alt in alt_list):
                continue

            # SNP 처리
            if len(ref) == 1 and all(len(alt) == 1 for alt in alt_list):
                for idx in set(gt_indices):
                    if idx == 0 or idx > len(alt_list):
                        continue
                    alt = alt_list[idx - 1]
                    out_var.write(b'\x00')
                    out_var.write(struct.pack('>I', pos))
                    out_var.write(bytes([encode_base2bit(ref)]))
                    out_var.write(bytes([encode_base2bit(alt)]))
                    out_meta.write(bytes([af, dp, encode_gt(gt_str)]))
                    variant_count += 1
                continue

            # INDEL 처리
            for idx in set(gt_indices):
                if idx == 0 or idx > len(alt_list):
                    continue
                alt = alt_list[idx - 1]

                if len(ref) > len(alt):  # Deletion
                    del_seq = ref[len(alt):]
                    del_len = len(del_seq)
                    encoded = encode_seq_fixed16(del_seq)
                    out_var.write(b'\x02')
                    out_var.write(struct.pack('>I', pos))
                    out_var.write(struct.pack('B', del_len))
                    out_var.write(encoded)

                elif len(ref) < len(alt):  # Insertion
                    insert_seq = alt[len(ref):]
                    if any(b not in 'ACGT' for b in insert_seq):
                        continue
                    insert_seq = insert_seq[:64]
                    encoded = encode_seq_fixed16(insert_seq)
                    out_var.write(b'\x01')
                    out_var.write(struct.pack('>I', pos))
                    out_var.write(struct.pack('B', len(insert_seq)))
                    out_var.write(encoded)

                else:
                    continue

                out_meta.write(bytes([af, dp, encode_gt(gt_str)]))
                variant_count += 1

    ref_block, mask_block = encode_reference_with_mask(fasta_path)

    with open(output_hex_filename, 'wb') as out_hex, \
         open(variant_bin, 'rb') as in_var, \
         open(meta_bin, 'rb') as in_meta:

        var_data = in_var.read()
        meta_data = in_meta.read()

        out_hex.write(ref_block)
        out_hex.write(mask_block)
        out_hex.write(var_data)
        out_hex.write(b'META')
        out_hex.write(struct.pack('>I', len(meta_data)))
        out_hex.write(meta_data)

    print(f" HEX created to: {output_hex_filename}")
    print(f" number of all mutations: {variant_count}개")
    return output_hex_filename

# 실행 예시
generate_ref_hex_with_mask("HG00157.chr11.vcf", "chr11.fasta", "chr11_fasta_with_ref_N_masking.hex")
