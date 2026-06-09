# pylint: disable=missing-module-docstring
# pylint: disable=missing-function-docstring
import difflib
import numpy as np
from numba import jit
from ipamd.public.models.sequence import Sequence
from ipamd.public.models.data import String, Scalar

def func(seq1: Sequence, seq2: Sequence, algorithm='needleman-wunsch', match_score=1, mismatch_score=-1, gap_score=-1):
    """
    Calculate the difference between two sequences
    :param seq1: the first sequence
    :param seq2: the second sequence
    :param algorithm: the algorithm to use
    :param match_score: the score for a match
    :param mismatch_score: the score for a mismatch
    :param gap_score: the score for a gap
    :return: a list of data objects
    """
    match algorithm:
        case 'needleman-wunsch':
            compute_func = needleman_wunsch
        case 'smith-waterman':
            compute_func = smith_waterman
        case 'diff':
            compute_func = diff
        case _:
            raise ValueError(f"Invalid algorithm: {algorithm}")
    seq1_text, seq2_text, score = compute_func(
        seq1.sequence,
        seq2.sequence,
        match_score,
        mismatch_score,
        gap_score
        )
    return [
        String(title='reference', data=seq1_text),
        String(title='target', data=seq2_text),
        Scalar(title='score', data=score)
    ]

def diff(seq1_text, seq2_text, match_score, mismatch_score, gap_score):
    matcher = difflib.SequenceMatcher(None, seq1_text, seq2_text)
    align1 = ""
    align2 = ""
    for opcode in matcher.get_opcodes():
        type_, i1, i2, j1, j2 = opcode
        if type_ == 'equal':
            align1 += seq1_text[i1: i2]
            align2 += seq2_text[j1: j2]
        elif type_ == 'replace':
            len1 = i2 - i1
            len2 = j2 - j1
            delta_len = abs(len1 - len2)
            if len1 > len2:
                target = 2
            elif len1 < len2:
                target = 1
            else:
                target = 0
            align1 += f"{seq1_text[i1:i2]}{'-' * delta_len if target == 1 else ''}"
            align2 += f"{seq2_text[j1:j2]}{'-' * delta_len if target == 2 else ''}"
        elif type_ == "insert":
            align1 += "-" * (j2 - j1)
            align2 += f"{seq2_text[j1:j2]}"
        elif type_ == "delete":
            align1 += f"{seq1_text[i1:i2]}"
            align2 += "-" * (i2 - i1)
    score = 0
    for ci, cj in zip(align1, align2):
        if ci == cj:
            score += match_score
        elif ci == '-' or cj == '-':
            score += gap_score
        else:
            score += mismatch_score
    return align1, align2, score

@jit
def needleman_wunsch(seq1_text, seq2_text, match_score, mismatch_score, gap_score):
    """
    Needleman-Wunsch algorithm.
    """
    m = len(seq1_text)
    n = len(seq2_text)
    score_matrix = np.zeros((m + 1, n + 1), dtype=np.int64)
    # initialize the first row and column
    for i in range(m + 1):
        score_matrix[i][0] = gap_score * i
    for j in range(n + 1):
        score_matrix[0][j] = gap_score * j
    # fill the score matrix
    for i in range(1, n + 1):
        for j in range(1, m + 1):
            # calculate the score for the current position
            score = match_score if seq1_text[i - 1] == seq2_text[j - 1] else mismatch_score
            current_match_score = score_matrix[i - 1][j - 1] + score
            delete_score = score_matrix[i - 1][j] + gap_score
            insert_score = score_matrix[i][j - 1] + gap_score
            score_matrix[i][j] = max(current_match_score, delete_score, insert_score)
    # trace back to get the alignment
    align1 = ""
    align2 = ""
    align_score = score_matrix[m][n]
    i = m
    j = n
    while i > 0 and j > 0:
        current_score = score_matrix[i][j]
        diag_score = score_matrix[i - 1][j - 1]
        up_score = score_matrix[i - 1][j]

        current_match_score = \
            match_score\
            if seq1_text[i - 1] == seq2_text[j - 1]\
            else mismatch_score

        if current_score == diag_score + current_match_score:
            align1 = seq1_text[i - 1] + align1
            align2 = seq2_text[j - 1] + align2
            i -= 1
            j -= 1
        elif current_score == up_score + gap_score:
            align1 = seq1_text[ i - 1] + align1
            align2 = "-" + align2
            i -= 1
        else:
            align1 = "-" + align1
            align2 = seq2_text[j - 1] + align2
            j -= 1

    while i > 0:
        align1 = seq1_text[i - 1] + align1
        align2 = "-" + align2
        i -= 1
    while j > 0:
        align1 = "-" + align1
        align2 = seq2_text[j - 1] + align2
        j -= 1

    return align1, align2, align_score

@jit
def smith_waterman(seq1_text, seq2_text, match_score, mismatch_score, gap_score):
    """
    Smith-Waterman algorithm
    """
    m = len(seq1_text)
    n = len(seq2_text)
    score_matrix = np.zeros((m + 1, n + 1), dtype=np.int64)

    for i in range(m + 1):
        score_matrix[i][0] = 0
    for j in range(n + 1):
        score_matrix[0][j] = 0

    max_score = 0
    max_pos = (m, n)

    for i in range(1, m + 1):
        for j in range(1, n + 1):
            score = match_score if seq1_text[i - 1] == seq2_text[j - 1] else mismatch_score
            current_match_score = score_matrix[i - 1][j - 1] + score
            delete_score = score_matrix[i - 1][j] + gap_score
            insert_score = score_matrix[i][j - 1] + gap_score

            score_matrix[i][j] = max(current_match_score, delete_score, insert_score, 0)

            if score_matrix[i][j] > max_score:
                max_score = score_matrix[i][j]
                max_pos = (i, j)

    align1 = ""
    align2 = ""
    align_score = max_score
    i, j = max_pos

    while i > 0 and j > 0 and score_matrix[i][j] > 0:
        current_score = score_matrix[i][j]
        diag_score = score_matrix[i - 1][j - 1]
        up_score = score_matrix[i - 1][j]

        current_match_score = (
            match_score
            if seq1_text[i - 1] == seq2_text[j - 1]
            else mismatch_score
        )

        if current_score == diag_score + current_match_score:
            align1 = seq1_text[i - 1] + align1
            align2 = seq2_text[j - 1] + align2
            i -= 1
            j -= 1
        elif current_score == up_score + gap_score:
            align1 = seq1_text[i - 1] + align1
            align2 = "-" + align2
            i -= 1
        else:
            align1 = "-" + align1
            align2 = seq2_text[j - 1] + align2
            j -= 1

    return align1, align2, align_score
