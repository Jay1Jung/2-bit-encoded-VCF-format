import os
import struct

# --- 인코딩 함수들 ---

def encode_af(value):
    try:
        val = float(value)
        return int(round(min(max(val, 0.0), 1.0) * 255))
    except:
        return 0

def encode_dp(value):
    try:
        return min(int(value), 255)
    except:
        return 0

def encode_gt(gt_str):
    gt = gt_str.replace('|', '/').split('/')

    # When GT is Lost or malformed
    if len(gt) != 2 or not all(g in {'0', '1'} for g in gt):
        return 2

    if gt == ['0', '0']:
        return 0
    elif '0' in gt and '1' in gt:
        return 1
    elif gt == ['1', '1']:
        return 3
    return 2  # fallback

# --- Extract metadata ---

def extract_variant_metadata(vcf_path):
    metadata = []

    base_dir = os.path.dirname(os.path.abspath(__file__))
    full_path = os.path.join(base_dir, vcf_path)

    with open(full_path, 'r') as f:
        for line in f:
            if line.startswith('#'):
                continue

            parts = line.strip().split('\t')
            if len(parts) < 10:
                continue

            info = parts[7]
            format_tags = parts[8].split(':')
            sample_fields = parts[9].split(':')

            # Filter VT=SNP or VT=INDEL 
            if 'VT=SNP' not in info and 'VT=INDEL' not in info:
                continue

            # Extract AF, DP 
            info_dict = dict(kv.split('=') for kv in info.split(';') if '=' in kv)
            af = encode_af(info_dict.get('AF', 0))
            dp = encode_dp(info_dict.get('DP', 0))

            # Extract GT 
            try:
                gt_idx = format_tags.index("GT")
                gt_str = sample_fields[gt_idx]
                gt = encode_gt(gt_str)
            except:
                gt = 2  # Treat as error when there is no GT field

            metadata.append(bytes([af, dp, gt]))

    return metadata

# ---  Storing Metadata binary file  ---

def write_metadata_block(output_path, metadata_list):
    with open(output_path, 'ab') as f:  # append mode
        f.write(b'META')  # 4-byte marker
        f.write(struct.pack('>I', len(metadata_list) * 3))  # 4-byte length
        for meta in metadata_list:
            f.write(meta)
    print(f" Saved metadata block: {len(metadata_list)}개 변이에 대해 3B씩 저장됨")

# --- 실행 예시 ---

if __name__ == '__main__':
    vcf_file = "HG00157.chr11.vcf"
    output_file = "variants_extracted.bin"  # 기존 BIN 파일에 이어붙이는 걸 권장

    metadata = extract_variant_metadata(vcf_file)
    write_metadata_block(output_file, metadata)
