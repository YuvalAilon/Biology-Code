
def kmers_of(sequence, k):
    if len(sequence) < k:
        raise Exception("k must be less than or equal to the length of the sequence")

    kmers = []
    for i in range(len(sequence) - (k - 1)):
        kmers.append(sequence[i: i + k])

    return kmers

def l_r_k_minus_1_mers(kmer):
    return [
        kmer[0: len(kmer) - 1],
        kmer[1: len(kmer)]
    ]