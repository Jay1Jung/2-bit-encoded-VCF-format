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

    print(f"🧪 FASTA 비교 결과:")
    print(f"• 총 길이: {total}")
    print(f"• 일치한 염기 수: {matches}")
    print(f"• 일치율: {match_rate:.2f}%")
    print(f"• 불일치 수: {total - matches}")

    if mismatches:
        print("\n🔍 불일치 예시 (최대 10개):")
        for idx, ref, alt in mismatches:
            print(f"  - 위치 {idx}: ref={ref} vs restored={alt}")

def preview_alignment(seq1, seq2, start=72960, window=100):
    print(f"\n🔬 시퀀스 시각 비교 [{start}:{start+window}]")
    print("ref:      ", seq1[start:start+window])
    print("restored: ", seq2[start:start+window])
    for i in range(window):
        print(" " if seq1[start+i] == seq2[start+i] else "^", end='')
    print("\n")

def align_sequences(seq1, seq2):
    print("\n🔧 edlib 정렬 중 (1Mbp) ...")
    result = edlib.align(seq1, seq2, mode='NW', task='distance')
    print(f"✅ edit distance: {result['editDistance']}")
    similarity = 100 * (1 - result['editDistance'] / max(len(seq1), len(seq2)))
    print(f"🧮 정렬 기반 유사도: {similarity:.4f}%")

if __name__ == "__main__":
    base_dir = os.path.dirname(os.path.abspath(__file__))

    original_path = os.path.join(base_dir, "chr11.fasta")
    restored_path = os.path.join(base_dir, "restored_chr11.fasta")

    print("📥 FASTA 파일 읽는 중...")
    ref_seq = normalize_sequence(read_fasta_sequence(original_path))
    restored_seq = normalize_sequence(read_fasta_sequence(restored_path))

    # 길이 맞추기
    min_len = min(len(ref_seq), len(restored_seq))
    ref_seq = ref_seq[:min_len]
    restored_seq = restored_seq[:min_len]

    print("🔬 FASTA 직접 비교 중...")
    compare_sequences(ref_seq, restored_seq)

    # (1) 시각 비교: 특정 위치에서 눈으로 mismatch 확인
    preview_alignment(ref_seq, restored_seq, start=72960, window=100)

    # (2) alignment 기반 정확도 (1Mbp만 잘라서 빠르게)
    align_sequences(ref_seq[:134741342], restored_seq[:134741342])
