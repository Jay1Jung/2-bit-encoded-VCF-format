import struct
import os

BASE_TABLE = ['A', 'C', 'G', 'T']

def decode_2bit_sequence(data, length):
    bases = bytearray()
    bit_idx = 0
    for i in range(length):
        byte = data[bit_idx // 8]
        shift = 6 - (bit_idx % 8)
        val = (byte >> shift) & 0b11
        bases.append(ord(BASE_TABLE[val]))
        bit_idx += 2
    return bases

def apply_variants(sequence: bytearray, variants):
    # Align in Reverse order 
    for v in sorted(variants, key=lambda x: x['pos'], reverse=True):
        if v['type'] == 'snp':
            sequence[v['pos']] = ord(v['alt'])
        elif v['type'] == 'ins':
            ins_seq = bytearray(v['alt'], 'ascii')
            sequence[v['pos']:v['pos']] = ins_seq  # insert
        elif v['type'] == 'del':
            del sequence[v['pos']:v['pos'] + v['len']]
    return sequence.decode('ascii')

def read_combined_file(hex_path, ref_length):
    with open(hex_path, 'rb') as f:
        # 1. Read 2-bit encoded reference
        ref_bytes = f.read((ref_length * 2 + 7) // 8)
        ref_seq = decode_2bit_sequence(ref_bytes, ref_length)

        # 2. Read 1-bit mask 
        mask_bytes = f.read((ref_length + 7) // 8)
        for i in range(ref_length):
            if (mask_bytes[i // 8] >> (7 - (i % 8))) & 1:
                ref_seq[i] = ord('N')

        # 3. Read variant records
        variants = []
        while True:
            type_byte = f.read(1)
            if not type_byte:
                break
            vtype = type_byte[0]

            # invalid variant type → terminate
            if vtype not in (0x00, 0x01, 0x02):
                print(f"Unknown variant type: {vtype:#x}, stopping.")
                break

            try:
                pos_bytes = f.read(4)
                if len(pos_bytes) < 4:
                    print("Unexpected EOF while reading position.")
                    break
                pos = struct.unpack('>I', pos_bytes)[0]

                if vtype == 0x00:  # SNP
                    ref = f.read(1)
                    alt = f.read(1)
                    if len(ref) < 1 or len(alt) < 1:
                        break
                    variants.append({'type': 'snp', 'pos': pos, 'alt': alt.decode()})

                elif vtype == 0x01:  # Insertion
                    len_bytes = f.read(2)
                    if len(len_bytes) < 2:
                        break
                    length = struct.unpack('>H', len_bytes)[0]
                    seq_bytes = f.read((length * 2 + 7) // 8)
                    if len(seq_bytes) < (length * 2 + 7) // 8:
                        break
                    ins_seq = decode_2bit_sequence(seq_bytes, length)
                    variants.append({'type': 'ins', 'pos': pos, 'alt': ins_seq.decode('ascii')})

                elif vtype == 0x02:  # Deletion
                    len_bytes = f.read(2)
                    if len(len_bytes) < 2:
                        break
                    length = struct.unpack('>H', len_bytes)[0]
                    del_bytes = f.read((length * 2 + 7) // 8)
                    if len(del_bytes) < (length * 2 + 7) // 8:
                        break
                    variants.append({'type': 'del', 'pos': pos, 'len': length})

            except Exception as e:
                print("Error while reading variant:", e)
                break

        print(f"Total variants read: {len(variants)}")
        return apply_variants(ref_seq, variants)

def write_fasta(sequence_str, output_path, header=">restored_chr11"):
    with open(output_path, 'w') as f:
        f.write(header + '\n')
        for i in range(0, len(sequence_str), 60):
            f.write(sequence_str[i:i+60] + '\n')

if __name__ == "__main__":
    REF_LENGTH = 135086622  

    base_dir = os.path.dirname(os.path.abspath(__file__))
    hex_path = os.path.join(base_dir, "combined_final.hex")
    output_path = os.path.join(base_dir, "restored_chr11.fasta")

    genome_seq = read_combined_file(hex_path, REF_LENGTH)
    write_fasta(genome_seq, output_path)
    print(f" FASTA restored and saved to: {output_path}")
