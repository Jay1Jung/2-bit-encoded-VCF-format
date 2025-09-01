import os
import struct

def parse_final_fasta_free_hex(filename):
    base_dir = os.path.dirname(os.path.abspath(__file__))
    path = os.path.join(base_dir, filename)

    with open(path, 'rb') as f:
        data = f.read()

    result = {}

    # find META section
    marker = b'META'
    idx = data.rfind(marker)
    if idx == -1:
        result["error"] = "META marker not found"
        return result

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
            i += 1 + 4 + 1 + 1
        elif tag == 0x01:  # Insertion
            i += 1 + 4 + 1 + 16
        elif tag == 0x02:  # Deletion
            i += 1 + 4 + 1 + 16
        else:
            break
        variant_count += 1

    meta_entries = [meta_data[j:j+3] for j in range(0, len(meta_data), 3)]
    sample_meta = meta_entries[:5]

    result["Variants (parsed)"] = variant_count
    result["META entries"] = len(meta_entries)
    result["Sample META (AF, DP, GT)"] = [tuple(entry) for entry in sample_meta]

    return result


# 실행
res = parse_final_fasta_free_hex("chr11_fasta_free_final.hex")
for k, v in res.items():
    print(f"{k}: {v}")
