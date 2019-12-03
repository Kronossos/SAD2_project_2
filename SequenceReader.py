import numpy as np
from collections import Counter
import pandas

from HMM_errors import *


class SequenceReader:
    """
    This class builds HMM from  protein or DNA sequence alignment.


    """

    def __init__(self, path_to_file, amino_acids=True, deletion_sign="-"):

        alignment_list = []

        with open(path_to_file) as seq_file:
            for line in seq_file:
                alignment_list.append(list(line.strip()))

        if not alignment_list:
            raise EmptyFile("File {} is empty!".format(path_to_file))

        self.alignment_list = alignment_list
        self.alignment_len = len(alignment_list)

        self.position_list = SequenceReader.transpose(alignment_list)
        self.position_len = len(self.position_list)

        self.deletion_sign = deletion_sign

        self.alignment_alphabet = (
            "A", "R", "N", "D", "C", "E", "Q", "G", "H", "I", "L", "K", "M", "F", "P", "S", "T", "W", "Y",
            "V") if amino_acids else ("A", "C", "G", "T")

        self.true_seq, self.match_emission, self.insert_emission = self.calculate_emissions()
        self.L = len(self.true_seq)
        self.trans = self.calculate_transitions()




    def calculate_transitions(self):
        match_transition = {"M": [0] * (self.L + 1), "D": [0] * (self.L + 1), "I": [0] * (self.L + 1)}
        insert_transition = {"M": [0] * (self.L + 1), "D": [0] * (self.L + 1), "I": [0] * (self.L + 1)}
        delete_transition = {"M": [float('nan')] + ([0] * self.L), "D": [float('nan')] + ([0] * self.L),
                             "I": [float('nan')] + ([0] * self.L)}

        trans = {"M": match_transition, "D": delete_transition, "I": insert_transition}

        for alignment_i, sequence in enumerate(self.alignment_list):
            last_seq_position = 0
            first_transition_index = "M"

            for true_seq_i, x in enumerate(self.true_seq):
                self.process_transitions(sequence, trans, last_seq_position, x, true_seq_i, first_transition_index)
                first_transition_index = "D" if sequence[x] == self.deletion_sign else "M"
                last_seq_position = x + 1
            else:
                if last_seq_position != 0:
                    true_seq_i += 1
                    self.process_transitions(sequence, trans, last_seq_position, len(sequence), true_seq_i,
                                             first_transition_index)


        return SequenceReader.divide_dict(trans)

    def process_insertion(self, insertion, insert_transition, sequence_sign, true_seq_i):
        last_transition_index = "D" if sequence_sign == self.deletion_sign else "M"
        if insertion:
            insert_transition["I"][true_seq_i] += len(insertion) - 1
            insert_transition[last_transition_index][true_seq_i] += 1
            last_transition_index = "I"
        return last_transition_index

    def process_transitions(self, sequence, trans, last_seq_position, x, true_seq_i, first_transition_index):
        insertion = SequenceReader.clear_list(sequence[last_seq_position:x])
        sign_in_sequence = sequence[x] if x < len(sequence) else True
        last_transition_index = self.process_insertion(insertion, trans["I"], sign_in_sequence, true_seq_i)
        trans[first_transition_index][last_transition_index][true_seq_i] += 1

    def calculate_emissions(self, deletion_sign="-", treshold=0.5):
        true_seq = []

        match_emission = []
        insert_emission = []

        full_occurrence_count = {}

        for i, position_i in enumerate(self.position_list):

            current_occurrence_count = Counter(position_i)  # UPDATE HERE!

            deletion_count = current_occurrence_count.get(deletion_sign, 0)

            if deletion_count / self.alignment_len < treshold:
                match_emission.append(SequenceReader.divide_list(self.build_count_column(current_occurrence_count)))
                insert_emission.append(SequenceReader.divide_list(self.build_count_column(full_occurrence_count)))
                true_seq.append(i)
                full_occurrence_count = {}
            else:
                full_occurrence_count = Counter(full_occurrence_count) + Counter(current_occurrence_count)
        else:
            insert_emission.append(SequenceReader.divide_list(self.build_count_column(full_occurrence_count)))

        self.L = len(true_seq)

        return true_seq, match_emission, insert_emission

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
    def transpose(iterable):
        return list(map(list, zip(*iterable)))

    @staticmethod
    def clear_list(list_to_clear, element_to_delete="-"):
        return [x for x in list_to_clear if x != element_to_delete]

    @staticmethod
    def divide_list(list_to_divide):
        list_sum = sum(list_to_divide)
        if not list_sum:
            return list_to_divide
        return [x / list_sum for x in list_to_divide]

    @staticmethod
    def divide_dict(dict_to_divide):
        for first, second in dict_to_divide.items():
            fraction = [sum(x) for x in zip(*second.values())]
            for key_of_list, list_to_divide in second.items():
                list_to_divide = [x / y if y != 0 else x for x, y in zip(list_to_divide, fraction)]
                dict_to_divide[first][key_of_list] = list_to_divide
        return dict_to_divide

def main():
    test = SequenceReader("test.txt", amino_acids=False)
    test2 = SequenceReader("test2.txt", amino_acids=False)


if __name__ == "__main__":
    # execute only if run as a script
    main()


test = SequenceReader("test.txt", amino_acids=False)
test2 = SequenceReader("test2.txt", amino_acids=False)