import os

def combine_blocks(reference_path, variant_path, output_path):
    with open(reference_path, 'rb') as ref_file:
        ref_data = ref_file.read()

    with open(variant_path, 'rb') as var_file:
        var_data = var_file.read()

    with open(output_path, 'wb') as out_file:
        out_file.write(ref_data)     #  Reference + Mask
        out_file.write(var_data)     #  Variants + Metadata

    print(f" Done creating final BIN file → {output_path}")
    print(f"Reference+Mask Size: {len(ref_data)} bytes")
    print(f" Variant+Meta size: {len(var_data)} bytes")
    print(f" Size of the whole file: {len(ref_data) + len(var_data)} bytes")

# 실행부
if __name__ == '__main__':
    base_dir = os.path.dirname(os.path.abspath(__file__))

    reference_file = os.path.join(base_dir, 'reference_with_mask.hex')
    variant_file   = os.path.join(base_dir, 'variants_extracted.bin')
    output_file    = os.path.join(base_dir, 'combined_final.hex')

    combine_blocks(reference_file, variant_file, output_file)
