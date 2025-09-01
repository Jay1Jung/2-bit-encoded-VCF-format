import random
import pandas as pd


original_path = "/Users/jayjung/Comp571/final project/chr11.fasta"
restored_path = "/Users/jayjung/Comp571/final project/check_efficiency.py"


sample_count = 5
sample_length = 1000


def load_clean_sequence(path):
    with open(path, 'r') as f:
        lines = [line.strip() for line in f if not line.startswith('>')]
    return ''.join(lines).upper().replace("N", "")


original_seq = load_clean_sequence(original_path)
restored_seq = load_clean_sequence(restored_path)


max_len = min(len(original_seq), len(restored_seq)) - sample_length


samples = []
for _ in range(sample_count):
    start = random.randint(0, max_len)
    orig_sample = original_seq[start:start + sample_length]
    rest_sample = restored_seq[start:start + sample_length]
    match_count = sum(1 for a, b in zip(orig_sample, rest_sample) if a == b)
    match_ratio = match_count / sample_length * 100
    samples.append({
        "샘플 시작 위치": start,
        "일치한 염기 수": match_count,
        "일치율 (%)": match_ratio
    })


df = pd.DataFrame(samples)
print(df.to_string(index=False))
