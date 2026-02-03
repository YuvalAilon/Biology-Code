new_indel = -5 # Not yet used
extend_indel = -1
indel_char = "-"

similarity_matrix = {
    frozenset({"A", "A"}): 1,
    frozenset({"A", "G"}): -1,
    frozenset({"A", "T"}): -1,
    frozenset({"A", "C"}): -1,

    frozenset({"G", "G"}): 1,
    frozenset({"G", "T"}): -1,
    frozenset({"G", "C"}): -1,

    frozenset({"T", "T"}): 1,
    frozenset({"T", "C"}): -1,

    frozenset({"C", "C"}): 1,
}

diagonal = "↖"
horizontal = "←"
vertical = "↑"


# Assumes both sequences are the same length and sequence1 has no gaps
def get_diff_score(sequence1, sequence2):
    diff_score = 0
    for i in range(len(sequence1)):
        if sequence2[i] == indel_char:
            diff_score += new_indel if i == 0 or sequence2[i - 1] != indel_char else extend_indel
        else:
            diff_score += similarity_matrix.get(frozenset({sequence1[i], sequence2[i]}))
    return diff_score


def needleman_wunsch(sequence1, sequence2):

    len_seq1 = len(sequence1)
    len_seq2 = len(sequence2)

    score_array = [[0 for _ in range(len_seq1 + 1)] for _ in range(len_seq2 + 1)]
    direction_array = [["" for _ in range(len_seq1 + 1)] for _ in range(len_seq2 + 1)]
    print(direction_array)
    for i in range(1, len_seq1 + 1):
        score_array[0][i] = score_array[0][i - 1] + extend_indel
    for i in range(1, len_seq2 + 1):
        score_array[i][0] = score_array[i - 1][0] + extend_indel

    for N in range(1, len_seq1 + 1):
        for M in range(1, len_seq2 + 1):
            max_score = max(
                score_array[M - 1][N - 1] + + similarity_matrix.get(frozenset({sequence1[N - 1], sequence2[M - 1]})),
                score_array[M][N - 1] + extend_indel,
                score_array[M - 1][N] + extend_indel
            )
            score_array[M][N] = max_score
            if max_score == score_array[M - 1][N - 1] + + similarity_matrix.get(frozenset({sequence1[N - 1], sequence2[M - 1]})):
                direction_array[M][N] = direction_array[M][N] + diagonal
            if max_score == score_array[M][N - 1] + extend_indel:
                direction_array[M][N] = direction_array[M][N] + horizontal
            if max_score == score_array[M - 1][N] + extend_indel:
                direction_array[M][N] = direction_array[M][N] + vertical


    print("    " + sequence1)
    for i in range(0, len(score_array)):
        if (i == 0):
            print("  " + str(score_array[i]))
        else:
            print(sequence2[i-1] + " " + str(score_array[i]))

    print("    " + sequence1)
    for i in range(0, len(direction_array)):
        if (i == 0):
            print("  " + str(direction_array[i]))
        else:
            print(sequence2[i - 1] + " " + str(direction_array[i]))

    reconstruct(sequence1, sequence2, direction_array)
    print("Score: " + str(score_array[M][N]))

def reconstruct(sequence1, sequence2, direction_array):
    sequence1_gaps = ""
    sequence2_gaps = ""

    len_seq1 = len(sequence1)
    len_seq2 = len(sequence2)

    currentN = len_seq1
    currentM = len_seq2

    while(currentM > 0 or currentN > 0):
            if diagonal in direction_array[currentM][currentN]:
                sequence1_gaps = prepend(sequence1_gaps, sequence1[currentN - 1])
                sequence2_gaps = prepend(sequence2_gaps, sequence2[currentM - 1])
                currentN -= 1
                currentM -= 1
            elif horizontal in direction_array[currentM][currentN] or currentM == 0:
                sequence1_gaps = prepend(sequence1_gaps, sequence1[currentN - 1])
                sequence2_gaps = prepend(sequence2_gaps, indel_char)
                currentN -= 1
            elif vertical in direction_array[currentM][currentN] or currentN == 0:
                sequence1_gaps = prepend(sequence1_gaps, indel_char)
                sequence2_gaps = prepend(sequence2_gaps, sequence2[currentM - 1])
                currentM -= 1
            print(sequence1_gaps)
            print(sequence2_gaps)
            print("M:" + str(currentM) + " | N:" + str(currentN))

def prepend(string, pre):
    return pre + string

needleman_wunsch(
    "CGATGCTAGCTAGCTAGCTA",
    "CGCGCGCGCTAAAAGCTAGC"
)
