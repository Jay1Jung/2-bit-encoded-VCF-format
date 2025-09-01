import os
def check_duplicate_indels(hex_path):
    seen = set()
    duplicates = []

    with open(hex_path, 'rb') as f:
        f.seek(0, os.SEEK_END)
        f.seek(max(0, f.tell() - 20000), os.SEEK_SET)
        data = f.read()

    i = 0
    while i < len(data):
        t = data[i]
        pos = int.from_bytes(data[i+1:i+5], 'big')
        length = data[i+5]

        if t == 0x00:
            key = ('del', pos, length)
            i += 6
        elif t == 0x01:
            encoded = data[i+6:i+22]
            key = ('ins', pos, length, encoded)
            i += 22
        else:
            break

        if key in seen:
            duplicates.append(key)
        else:
            seen.add(key)

    if duplicates:
        print(" Identified overlapping ins/del:")
        for d in duplicates:
            print(f"  {d}")
    else:
        print(" No identified overlaaping ins/del.")


check_duplicate_indels("/Users/jayjung/Comp571/final project/chr11_combined_binary.hex")