from enum import Enum

new_indel = -5
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

space_used = "■"


class reconstruct_state(Enum):
    VERTICAL = 1
    HORIZONTAL = 2
    MATCH = 3


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

    match_array = [[0 for _ in range(len_seq1 + 1)] for _ in range(len_seq2 + 1)]
    horizontal_indel_array = [[0 for _ in range(len_seq1 + 1)] for _ in range(len_seq2 + 1)]
    vertical_indel_array = [[0 for _ in range(len_seq1 + 1)] for _ in range(len_seq2 + 1)]

    direction_array = [["" for _ in range(len_seq1 + 1)] for _ in range(len_seq2 + 1)]

    match_array[1][0] = new_indel + extend_indel
    match_array[0][1] = new_indel + extend_indel

    for i in range(1, len_seq1 + 1):
        if i != 1:
            match_array[0][i] = match_array[0][i - 1] + extend_indel
        vertical_indel_array[0][i] = float('-inf')
        direction_array[0][i] = horizontal
    for i in range(1, len_seq2 + 1):
        if i != 1:
            match_array[i][0] = match_array[i - 1][0] + extend_indel
        horizontal_indel_array[i][0] = float('-inf')
        direction_array[i][0] = vertical

    direction_array[0][0] = diagonal

    for N in range(1, len_seq1 + 1):
        for M in range(1, len_seq2 + 1):
            horizontal_max = max_horizontal(horizontal_indel_array, match_array, M, N)
            horizontal_indel_array[M][N] = horizontal_max

            vertical_max = max_vertical(vertical_indel_array, match_array, M, N)
            vertical_indel_array[M][N] = vertical_max

            match_max = max_match(
                match_array, vertical_indel_array, horizontal_indel_array, M, N, sequence1, sequence2
            )
            match_array[M][N] = match_max

            max_score = max(horizontal_max, vertical_max, match_array[M - 1][N - 1] + similarity_matrix.get(frozenset(
                {sequence1[N - 1], sequence2[M - 1]}
            )))

            if max_score == horizontal_max:
                direction_array[M][N] += horizontal
            if max_score == vertical_max:
                direction_array[M][N] += vertical
            if max_score == match_array[M - 1][N - 1] + similarity_matrix.get(
                    frozenset({sequence1[N - 1], sequence2[M - 1]})):
                direction_array[M][N] += diagonal

    reconstruct(sequence1, sequence2, direction_array, vertical_indel_array, horizontal_indel_array, match_array)

    print("Score: " + str(match_array[M][N]) + "\n")

    print(pretty_print_array(direction_array, format_square, sequence1, sequence2))



def max_vertical(vertical_array, match_array, M, N):
    return max(
        match_array[M - 1][N] + new_indel + extend_indel,
        vertical_array[M - 1][N] + extend_indel
    )


def max_match(match_array, vertical_array, horizontal_array, M, N, sequence1, sequence2):
    return max(
        match_array[M - 1][N - 1] + similarity_matrix.get(frozenset({sequence1[N - 1], sequence2[M - 1]})),
        horizontal_array[M][N],
        vertical_array[M][N]
    )


def max_horizontal(horizontal_array, match_array, M, N):
    return max(
        match_array[M][N - 1] + new_indel + extend_indel,
        horizontal_array[M][N - 1] + extend_indel
    )


def reconstruct(sequence1, sequence2, direction_array, vertical_array, horizontal_array, match_array):
    sequence1_gaps = ""
    sequence2_gaps = ""

    len_seq1 = len(sequence1)
    len_seq2 = len(sequence2)

    currentN = len_seq1
    currentM = len_seq2

    state = reconstruct_state.MATCH

    while (currentM > 0 or currentN > 0):
        direction_array[currentM][currentN] += space_used
        directions = direction_array[currentM][currentN]
        match state:
            case reconstruct_state.MATCH:
                if diagonal in directions:
                    sequence1_gaps = prepend(sequence1_gaps, sequence1[currentN - 1])
                    sequence2_gaps = prepend(sequence2_gaps, sequence2[currentM - 1])
                    currentN -= 1
                    currentM -= 1
                elif horizontal in directions:
                    state = reconstruct_state.HORIZONTAL
                elif vertical in directions:
                    state = reconstruct_state.VERTICAL
            case reconstruct_state.HORIZONTAL:
                sequence1_gaps = prepend(sequence1_gaps, sequence1[currentN - 1])
                sequence2_gaps = prepend(sequence2_gaps, indel_char)

                if horizontal_array[currentM][currentN] == (
                        match_array[currentM][currentN - 1] + new_indel + extend_indel):
                    state = reconstruct_state.MATCH

                currentN -= 1

            case reconstruct_state.VERTICAL:
                sequence1_gaps = prepend(sequence1_gaps, indel_char)
                sequence2_gaps = prepend(sequence2_gaps, sequence2[currentM - 1])

                if vertical_array[currentM][currentN] == (
                        match_array[currentM - 1][currentN] + new_indel + extend_indel):
                    state = reconstruct_state.MATCH

                currentM -= 1
    print("*" * (max(len_seq1, len_seq2) + 4))
    print_match(sequence1_gaps, sequence2_gaps)
    print("*" * (max(len_seq1, len_seq2) + 4))


def print_match(sequence1, sequence2):
    match_string = ""
    for i in range(0, len(sequence1)):
        if sequence1[i] == sequence2[i] and sequence1[i] != indel_char:
            match_string += "|"
        elif sequence1[i] != sequence2[i] and sequence1[i] != indel_char and sequence2[i] != indel_char:
            match_string += "*"
        else:
            match_string += " "
    print(sequence1)
    print(match_string)
    print(sequence2)


def prepend(string, pre):
    return pre + string


def pretty_print_array(array, formatter, sequence1, sequence2):
    matrix = ""
    for i in range(len(array)):
        is_top = i == 0
        is_bottom = i == len(array) - 1
        matrix += print_matrix_row(list(map(formatter, array[i])), is_top, is_bottom, sequence1, sequence2[i-1]) + "\n"
    return matrix


def format_square(square):
    grid_square = ["", ""]
    grid_square[0] = diagonal if diagonal in square else " "
    grid_square[0] += vertical if vertical in square else " "

    grid_square[1] = horizontal if horizontal in square else " "
    grid_square[1] += space_used if space_used in square else " "

    return grid_square


def print_matrix_row(row, is_top, is_bottom, sequence1, sequence2_letter):
    side_letter_gap = "  "
    print_string = ""
    if is_top:
        print_string += side_letter_gap + (" " * (len(row[0][0]) + 1))
        for char in sequence1:
            print_string += (" " * (len(row[0][0]) + 1)) + char
        print_string += "\n"
    if not is_top:
        print_string += sequence2_letter + " "
    print_string += "│" if not is_top else side_letter_gap + divider_line("┬", len(row), len(row[0][0]) + 1) + "\n" + side_letter_gap + "│"
    for i in range(len(row[0])):
        for square in row:
            print_string += square[i] + " │"
        print_string += "\n"
        if i != len(row[0]) - 1:
            print_string += side_letter_gap + "│"
    bottom_char = "┴" if is_bottom else "┼"
    print_string += side_letter_gap + divider_line(bottom_char, len(row), len(row[0][0]) + 1)
    return print_string


def divider_line(divider_char, length, cell_width):
    return (divider_char + "─" * (cell_width)) * length + divider_char


needleman_wunsch(
    "CAAATG",
    "CATG"
)
