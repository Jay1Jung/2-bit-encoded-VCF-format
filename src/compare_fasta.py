import os
import edlib  # pip install edlib

def read_fasta_sequence(fasta_path):
    seq = []
    with open(fasta_path, 'r') as f:
        for line in f:
            if line.startswith('>'):
                continue
            seq.append(line.strip().upper())
    return ''.join(seq)

def normalize_sequence(seq):
    return seq.upper().replace('.', 'N')

def compare_sequences(seq1, seq2, show_mismatch_examples=True, max_examples=10):
    min_len = min(len(seq1), len(seq2))
    matches = 0
    mismatches = []

    for i in range(min_len):
        if seq1[i] == seq2[i]:
            matches += 1
        else:
            mismatches.append((i, seq1[i], seq2[i]))
            if show_mismatch_examples and len(mismatches) >= max_examples:
                break

    total = max(len(seq1), len(seq2))
    match_rate = matches / total * 100

    print(f" FASTA 비교 결과:")
    print(f"• Total Length: {total}")
    print(f"• Number of matched base: {matches}")
    print(f"• Match rate: {match_rate:.2f}%")
    print(f"• Number of mistmatches: {total - matches}")

    if mismatches:
        print("\n mismatches (max 10):")
        for idx, ref, alt in mismatches:
            print(f"  - 위치 {idx}: ref={ref} vs restored={alt}")

def preview_alignment(seq1, seq2, start=72960, window=100):
    print(f"\n Sequence comparison [{start}:{start+window}]")
    print("ref:      ", seq1[start:start+window])
    print("restored: ", seq2[start:start+window])
    for i in range(window):
        print(" " if seq1[start+i] == seq2[start+i] else "^", end='')
    print("\n")

def align_sequences(seq1, seq2):
    print("\n🔧 edlib aligning (1Mbp) ...")
    result = edlib.align(seq1, seq2, mode='NW', task='distance')
    print(f" edit distance: {result['editDistance']}")
    similarity = 100 * (1 - result['editDistance'] / max(len(seq1), len(seq2)))
    print(f" Alignment-based similarity: {similarity:.4f}%")

if __name__ == "__main__":
    base_dir = os.path.dirname(os.path.abspath(__file__))

    original_path = os.path.join(base_dir, "chr11.fasta")
    restored_path = os.path.join(base_dir, "restored_chr11.fasta")

    print(" Reading FASTA file")
    ref_seq = normalize_sequence(read_fasta_sequence(original_path))
    restored_seq = normalize_sequence(read_fasta_sequence(restored_path))

    min_len = min(len(ref_seq), len(restored_seq))
    ref_seq = ref_seq[:min_len]
    restored_seq = restored_seq[:min_len]

    print("🔬 Comparing with FASTA")
    compare_sequences(ref_seq, restored_seq)

    # Visual comparison: manually inspecting mismatches at specific positions
    preview_alignment(ref_seq, restored_seq, start=72960, window=100)

    # Alignment-based accuracy (quick evaluation on a 1 Mbp segment)
    align_sequences(ref_seq[:134741342], restored_seq[:134741342])
