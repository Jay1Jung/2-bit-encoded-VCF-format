import os
from Bio import SeqIO

# ACGT -> encoded as 2bits(when there is N, it is converted to 00) 
def encode_base2bit(base):
    return {
        'A': 0b00,
        'C': 0b01,
        'G': 0b10,
        'T': 0b11
    }.get(base, 0b00)  

# Convert FASTA → reference block + mask block 
def encode_reference_and_mask(fasta_path):
    record = next(SeqIO.parse(fasta_path, "fasta"))
    seq = str(record.seq).upper()

    bitstream = bytearray()
    mask_bits = []

    bits = 0
    count = 0

    for base in seq:
        is_N = base not in 'ACGT'
        mask_bits.append(1 if is_N else 0)

        bits = (bits << 2) | encode_base2bit(base)
        count += 1

        if count % 4 == 0:
            bitstream.append(bits & 0xFF)
            bits = 0

    if count % 4 != 0:
        bits <<= 2 * (4 - (count % 4))  # padding
        bitstream.append(bits & 0xFF)

    # Compressing mask bit
    mask_bytes = bytearray()
    for i in range(0, len(mask_bits), 8):
        chunk = mask_bits[i:i+8]
        byte = 0
        for bit in chunk:
            byte = (byte << 1) | bit
        byte <<= (8 - len(chunk))  # padding
        mask_bytes.append(byte)

    return bitstream, mask_bytes

def write_reference_and_mask(output_path, ref_bytes, mask_bytes):
    with open(output_path, 'wb') as f:
        f.write(ref_bytes)
        f.write(mask_bytes)

    print(f" Reference + Mask saved at: → {output_path}")
    print(f" Reference bytes: {len(ref_bytes)}")
    print(f" Mask bytes: {len(mask_bytes)}")


if __name__ == '__main__':
    base_dir = os.path.dirname(os.path.abspath(__file__))


    fasta_file = os.path.join(base_dir, "chr11.fasta") 
    output_file = os.path.join(base_dir, "reference_with_mask.hex")

    print(f"📂 FASTA location: {fasta_file}")

    ref_bytes, mask_bytes = encode_reference_and_mask(fasta_file)
    write_reference_and_mask(output_file, ref_bytes, mask_bytes)
