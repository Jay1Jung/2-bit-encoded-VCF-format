import os

def append_bin_to_hex_with_os(hex_filename, bin_filename):
    base_dir = os.path.dirname(__file__)
    
    hex_path = os.path.join(base_dir, hex_filename)
    bin_path = os.path.join(base_dir, bin_filename)

    if not os.path.exists(hex_path):
        print(f" '{hex_path}' File DNE.")
        return
    if not os.path.exists(bin_path):
        print(f" '{bin_path}' File DNE.")
        return

    
    with open(hex_path, 'ab') as hex_file, open(bin_path, 'rb') as bin_file:
        data = bin_file.read()
        hex_file.write(data)
        print(f" {len(data)} bytes appended from '{bin_filename}' to '{hex_filename}'")

#example
append_bin_to_hex_with_os("chr11_combined_binary.hex", "chr11_indels.bin")