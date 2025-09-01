import os
hex_path = "/Users/jayjung/Comp571/final project/chr11_combined_binary.hex"

def decode_2bit_fixed16(bits, length):
    base_map = ['A', 'C', 'G', 'T']
    result = ""
    for i in reversed(range(length)):
        base = (bits >> (i * 2)) & 0b11
        result += base_map[base]
    return result

def read_last_indels(hex_path, num_records=5):
    with open(hex_path, 'rb') as f:
        f.seek(0, os.SEEK_END)
        file_size = f.tell()

        f.seek(max(0, file_size - 5000), os.SEEK_SET)  # 끝에서 5KB 읽기
        data = f.read()

    i = 0
    count = 0
    while i < len(data) and count < num_records:
        type_byte = data[i]
        pos = int.from_bytes(data[i+1:i+5], 'big')
        length = data[i+5]

        if type_byte == 0x00:
            print(f"[Deletion] Pos={pos}, Len={length}")
            i += 6
        elif type_byte == 0x01:
            encoded_bytes = data[i+6:i+22]
            bits = int.from_bytes(encoded_bytes, 'big')
            seq = decode_2bit_fixed16(bits, length)
            print(f"[Insertion] Pos={pos}, Len={length}, Seq={seq}")
            i += 22
        else:
            print(f"[?] Unknown type {type_byte} at offset {i}")
            break

        count += 1

read_last_indels(hex_path, num_records=2000)