import numpy as np
from collections import Counter
import pandas as pd
from math import log10

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
        self.amino_acids = amino_acids

        self.true_seq, self.match_emission, self.insert_emission = self.calculate_emissions()

        self.match_emission = pd.DataFrame(SequenceReader.transpose(self.match_emission))
        self.match_emission.index = self.alignment_alphabet
        self.insert_emission = pd.DataFrame(SequenceReader.transpose(self.insert_emission))
        self.insert_emission.index = self.alignment_alphabet

        self.L = len(self.true_seq)
        self.trans = self.calculate_transitions()

    def calculate_transitions(self):
        match_transition = {"M": [0] * (self.L + 1), "D": [0] * (self.L + 1), "I": [0] * (self.L + 1)}
        insert_transition = {"M": [0] * (self.L + 1), "D": [0] * (self.L + 1), "I": [0] * (self.L + 1)}
        # delete_transition = {"M": [float('nan')] + ([0] * self.L), "D": [float('nan')] + ([0] * self.L),
        #                      "I": [float('nan')] + ([0] * self.L)}
        delete_transition = {"M": [float(0)] + ([0] * self.L), "D": [float(0)] + ([0] * self.L),
                             "I": [float(0)] + ([0] * self.L)}

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

        # match_emission.insert(0,[float("NaN")]*len(self.alignment_alphabet))
        match_emission.insert(0, [float(0)] * len(self.alignment_alphabet))

        return true_seq, match_emission, insert_emission

    def build_count_column(self, occurrence_count):
        match_column = []
        for letter in self.alignment_alphabet:
            match_column.append(occurrence_count.get(letter, 0))
        return match_column

    def forward(self, sequence):

        FM = np.zeros(shape=(len(sequence) + 1, self.L + 1), dtype=np.float64)
        FI = np.zeros(shape=(len(sequence) + 1, self.L + 1), dtype=np.float64)
        FD = np.zeros(shape=(len(sequence) + 1, self.L + 1), dtype=np.float64)

        # row_len = seqeunce len = i
        # column_len = model len = j

        FI[0, 0] = np.log(self.insert_emission.loc[sequence[0]][0])
        # FM[1, 1] = np.log(self.match_emission.loc[sequence[1]][1]) + np.log(
        #     self.trans["M"]["M"][0] * np.exp(FM[0, 0]) + self.trans["I"]["M"][0] * np.exp(FI[0, 0])
        # )
        #
        # FI[1, 1] = np.log(self.insert_emission.loc[sequence[1]][1]) + np.log(
        #     self.trans["M"]["I"][0] * np.exp(FM[0, 1]) + self.trans["I"]["I"][0] * np.exp(FI[0, 1])
        # )
        #
        # FD[1, 1] = np.log(self.trans["M"]["D"][0] * np.exp(FM[1, 0]) + self.trans["I"]["D"][0] * np.exp(FI[1, 0]))

        for i, sign in enumerate(sequence):
            i += 1
            for j in range(1, self.L + 1):
                FM[i, j] = np.log(self.match_emission.loc[sign][j] / (1 / len(self.alignment_alphabet))) + np.log(
                    (self.trans["M"]["M"][j - 1] * np.exp(FM[i - 1][j - 1])) +
                    (self.trans["I"]["M"][j - 1] * np.exp(FI[i - 1][j - 1])) +
                    (self.trans["D"]["M"][j - 1] * np.exp(FD[i - 1][j - 1]))
                )
                FI[i, j] = np.log(self.insert_emission.loc[sign][j] / (1 / len(self.alignment_alphabet))) + np.log(
                    (self.trans["M"]["I"][j - 1] * np.exp(FM[i - 1][j])) +
                    (self.trans["I"]["I"][j - 1] * np.exp(FI[i - 1][j])) +
                    (self.trans["D"]["I"][j - 1] * np.exp(FD[i - 1][j]))
                )
                FD[i, j] = np.log(
                    (self.trans["M"]["D"][j - 1] * np.exp(FM[i][j - 1])) +
                    (self.trans["I"]["D"][j - 1] * np.exp(FI[i][j - 1])) +
                    (self.trans["D"]["D"][j - 1] * np.exp(FD[i][j - 1]))
                )
        return FM[i][j]

    @staticmethod
    def compare(sequence, HMM_model_1, HMM_model_2=False):
        if not HMM_model_2:
            HMM_model_2 = SequenceReader.random_model(HMM_model_1.position_len, HMM_model_1.alignment_len,
                                                      HMM_model_1.amino_acids, HMM_model_1.deletion_sign)
        return log10(HMM_model_1.forward(sequence) / HMM_model_2.forward(sequence))

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
        list_to_divide = [x + 1 for x in list_to_divide]
        list_sum = sum(list_to_divide)
        return [x / list_sum for x in list_to_divide]

    @staticmethod
    def divide_dict(dict_to_divide):
        for first, second in dict_to_divide.items():
            fraction = [sum(x) + 1 * 3 for x in zip(*second.values())]
            for key_of_list, list_to_divide in second.items():
                list_to_divide = [(x + 1) / y for x, y in zip(list_to_divide, fraction)]
                dict_to_divide[first][key_of_list] = list_to_divide
        return dict_to_divide


def main():
    atp = SequenceReader("ATPases.txt")
    gtp = SequenceReader("GTP_binding_proteins.txt")

    print("ATPase model len: {}\n GTP model len: {}".format(atp.L, gtp.L))
    print("ATPase T, mE, iE shapes: {},{},{}\n GTP T, mE, iE shapes: {},{},{}".format((9, len(atp.trans["M"]["M"])),
                                                                                      atp.match_emission.shape,
                                                                                      atp.insert_emission.shape,
                                                                                      (9, len(gtp.trans["M"]["M"])),
                                                                                      gtp.match_emission.shape,
                                                                                      gtp.insert_emission.shape))

    print("----------Position 50 of ATPases HMM model.------------")
    print("Match emission: {}".format(atp.match_emission.loc[:, 50]))
    print("Insert emission: {}".format(atp.insert_emission.loc[:, 50]))
    for x1, y1 in atp.trans.items():
        for x2, y2 in y1.items():
            print("{}-{}: {}".format(x1, x2, y2[50]))

    print("----------Position 50 of GTP binding proteins HMM model.------------")
    print("Match emission: {}".format(gtp.match_emission.loc[:, 50]))
    print("Insert emission: {}".format(gtp.insert_emission.loc[:, 50]))
    for x1, y1 in gtp.trans.items():
        for x2, y2 in y1.items():
            print("{}-{}: {}".format(x1, x2, y2[50]))

    print("----------Position 49 of ATPases HMM model.------------")
    print("Match emission: {}".format(atp.match_emission.loc[:, 49]))
    print("Insert emission: {}".format(atp.insert_emission.loc[:, 49]))
    for x1, y1 in atp.trans.items():
        for x2, y2 in y1.items():
            print("{}-{}: {}".format(x1, x2, y2[49]))

    print("----------Position 49 of GTP binding proteins HMM model.------------")
    print("Match emission: {}".format(gtp.match_emission.loc[:, 49]))
    print("Insert emission: {}".format(gtp.insert_emission.loc[:, 49]))
    for x1, y1 in gtp.trans.items():
        for x2, y2 in y1.items():
            print("{}-{}: {}".format(x1, x2, y2[49]))

    y = []

    import matplotlib.pyplot as plt

    with open("Unclassified_proteins.txt") as file:
        for line in file:
            atp_score = atp.forward(line.strip())
            gtp_score = gtp.forward(line.strip())
            print()
            print(line.strip())
            print("ATP score: {}".format(atp_score))
            print("GTP score: {}".format(gtp_score))
            y.append(atp_score - gtp_score)

    fig = plt.figure(figsize=(10, 6))

    y_plus = [x if x >= 0 else 0 for x in y]
    y_minus = [x if x < 0 else 0 for x in y]

    x = [x for x in range(len(y))]
    plt.bar(x, y, tick_label=x, width=0.5)

    fig = plt.figure(figsize=(10, 6))
    ax = fig.add_subplot(111)
    ax.bar(x, y_minus, tick_label=x, width=0.5, color='r', label='GTP binding proteins')
    ax.bar(x, y_plus, tick_label=x, width=0.5, color='b', label='ATPases')
    plt.title('Family membership')
    ax.legend()
    plt.show()


if __name__ == "__main__":
    # execute only if run as a script
    main()
