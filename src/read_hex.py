import struct

# 새로운 hex 파일 경로
hex_path = "/Users/jayjung/Comp571/final project/chr11_fasta_with_ref_N_masking.hex"

def parse_hex_footer_and_summary(path):
    with open(path, 'rb') as f:
        data = f.read()

    marker = b'META'
    idx = data.rfind(marker)
    if idx == -1:
        return " META marker not found. Not properly formatted."

    length_start = idx + 4
    meta_len = struct.unpack('>I', data[length_start:length_start+4])[0]
    meta_start = length_start + 4
    meta_data = data[meta_start:meta_start+meta_len]

    variant_data = data[:idx]
    variant_count = 0
    i = 0
    while i < len(variant_data):
        tag = variant_data[i]
        if tag == 0x00:  # SNP
            i += 1 + 4 + 16
        elif tag == 0x01:  # Insertion
            i += 1 + 4 + 1 + 16
        elif tag == 0x02:  # Deletion
            i += 1 + 4 + 1
        else:
            break  # Unknown tag
        variant_count += 1

    meta_entries = [meta_data[i:i+3] for i in range(0, len(meta_data), 3)]
    sample_meta = meta_entries[:5]

    return {
        "Variants (parsed)": variant_count,
        "META entries": len(meta_entries),
        "Sample META (AF, DP, GT)": [tuple(e) for e in sample_meta]
    }

print(parse_hex_footer_and_summary(hex_path))
