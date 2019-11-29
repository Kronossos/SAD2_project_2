import re
import numpy as np
import sys
from collections import Counter
import pandas

from HMM_errors import *


class SequenceReader:
    """

    """

    def __init__(self, path_to_file, amino_acids=True, *args, **kwargs):

        alignment_list = []

        with open(path_to_file) as seq_file:
            for line in seq_file:
                alignment_list.append(list(line.strip()))

        if not alignment_list:
            raise EmptyFile("File {} is empty!".format(path_to_file))

        self.alignment_list = alignment_list
        self.alignment_len = len(alignment_list)

        self.position_list = list(map(list, zip(*alignment_list)))
        self.position_len = len(self.position_list)

        self.alignment_alphabet = (
            "A", "R", "N", "D", "C", "E", "Q", "G", "H", "I", "L", "K", "M", "F", "P", "S", "T", "W", "Y",
            "V") if amino_acids else ("A", "C", "G", "T", "U")

        self.calculate_emissions(*args, **kwargs)
        self.calculate_transitions(*args, **kwargs)

    def calculate_transitions(self, *args, **kwargs):






        for i, letter in enumerate(self.alignment_list):
            before_seq = SequenceReader.clear_list(letter[0:self.true_seq[0]])






        for x, y in SequenceReader.pairwise(self.true_seq):
            for sequence in self.alignment_list:
                pass

    def calculate_emissions(self, deletion_sign="-", treshold=0.5):

        true_seq = []

        mE = []
        iE = []

        full_occurrence_count = {}

        for i, position_i in enumerate(self.position_list):

            current_occurrence_count = Counter(position_i)  # UPDATE HERE!

            deletion_count = current_occurrence_count.get(deletion_sign, 0)

            if deletion_count / self.alignment_len > treshold:
                mE.append(self.build_count_column(current_occurrence_count))
                iE.append(self.build_count_column(full_occurrence_count))
                true_seq.append(i)
                full_occurrence_count = {}
            else:
                full_occurrence_count = Counter(full_occurrence_count) + Counter(current_occurrence_count)
        else:
            iE.append(self.build_count_column(full_occurrence_count))

        self.L = len(true_seq)
        self.true_seq = true_seq
        self.mE = mE
        self.iE = iE

    def build_count_column(self, occurrence_count):
        match_column = []
        for letter in self.alignment_alphabet:
            match_column.append(occurrence_count.get(letter, 0))
        return match_column

    @staticmethod
    def pairwise(iterable):
        "s -> (s0,s1), (s1,s2), (s2, s3), ..."

        from itertools import tee

        a, b = tee(iterable)
        next(b, None)
        return zip(a, b)

    @staticmethod
    def clear_list(list_to_clear, element_to_delete = "-"):
        return [x for x in list_to_clear if x != element_to_delete]
