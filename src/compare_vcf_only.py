import os
import struct

IDX_TO_BASE = ['A', 'C', 'G', 'T']

#VCF parser: extract only VT = SNP or VT = INDEL
def parse_vcf(vcf_path):
    vcf_variants = {}
    with open(vcf_path, 'r') as f:
        for line in f:
            if line.startswith('#'):
                continue
            parts = line.strip().split('\t')
            pos = int(parts[1])
            ref = parts[3]
            alt_field = parts[4]
            info = parts[7]

            # Define the type of mutation from VT field
            if 'VT=SNP' in info:
                for alt in alt_field.split(','):
                    vcf_variants[pos] = ('SNP', ref, alt)
            elif 'VT=INDEL' in info:
                for alt in alt_field.split(','):
                    delta = len(ref) - len(alt)
                    if delta > 0:
                        vcf_variants[pos] = ('DEL', None, f"-{delta}bp")
                    elif delta < 0:
                        vcf_variants[pos] = ('INS', None, f"+{abs(delta)}bp")
    return vcf_variants

# HEX parser
def parse_hex_variants(hex_path, ref_len):
    with open(hex_path, 'rb') as f:
        hex_data = f.read()

    meta_index = hex_data.find(b'META')
    ref_bytes = (ref_len + 3) // 4
    mask_bytes = (ref_len + 7) // 8
    variant_block = hex_data[ref_bytes + mask_bytes:meta_index]

    hex_variants = {}
    i = 0
    while i + 5 < len(variant_block):
        variant_type = variant_block[i]
        pos = struct.unpack('>I', variant_block[i + 1:i + 5])[0]
        if variant_type == 0x00 and i + 7 <= len(variant_block):  # SNP
            ref_idx = variant_block[i + 5]
            alt_idx = variant_block[i + 6]
            ref = IDX_TO_BASE[ref_idx] if ref_idx < 4 else '?'
            alt = IDX_TO_BASE[alt_idx] if alt_idx < 4 else '?'
            hex_variants[pos] = ('SNP', ref, alt)
            i += 7
        elif variant_type == 0x01 and i + 22 <= len(variant_block):  # INS
            length = variant_block[i + 5]
            hex_variants[pos] = ('INS', None, f"+{length}bp")
            i += 22
        elif variant_type == 0x02 and i + 22 <= len(variant_block):  # DEL
            length = variant_block[i + 5]
            hex_variants[pos] = ('DEL', None, f"-{length}bp")
            i += 22
        else:
            break
    return hex_variants

# Compare
def compare_variants(vcf_dict, hex_dict):
    common = 0
    mismatched = []
    missing_in_hex = []

    for pos in vcf_dict:
        if pos in hex_dict:
            vcf_type, vcf_ref, vcf_alt = vcf_dict[pos]
            hex_type, hex_ref, hex_alt = hex_dict[pos]

            if vcf_type == hex_type:
                if vcf_type == 'SNP':
                    if vcf_ref == hex_ref and vcf_alt == hex_alt:
                        common += 1
                    else:
                        mismatched.append((pos, vcf_dict[pos], hex_dict[pos]))
                else:  # INS/DEL
                    if vcf_alt == hex_alt:
                        common += 1
                    else:
                        mismatched.append((pos, vcf_dict[pos], hex_dict[pos]))
            else:
                mismatched.append((pos, vcf_dict[pos], hex_dict[pos]))
        else:
            missing_in_hex.append((pos, vcf_dict[pos]))
    print("\nVariant comparison results:")
    print(f" Number of matching variants: {common}")
    print(f" Number of variants missing in HEX: {len(missing_in_hex)}")
    print(f" Number of format or ALT mismatches: {len(mismatched)}\n")

    if mismatched:
        print("Example - Format or ALT mismatch:")
        for x in mismatched[:5]:
            print(x)
    if missing_in_hex:
        print("\nExample - Variants missing in HEX:")
        for x in missing_in_hex[:5]:
            print(x)

# 실행부
if __name__ == "__main__":
    base_dir = os.path.dirname(os.path.abspath(__file__))
    vcf_path = os.path.join(base_dir, "HG00157.chr11.vcf")
    hex_path = os.path.join(base_dir, "chr11_fasta_with_ref_N_masking.hex")
    ref_len = 134534784  # chr11 길이

    vcf_variants = parse_vcf(vcf_path)
    hex_variants = parse_hex_variants(hex_path, ref_len)
    compare_variants(vcf_variants, hex_variants)
