from Bio import SeqIO, SeqFeature, SeqRecord
from collections import Counter
import csv
from random import random
from pandas import DataFrame
import matplotlib.pyplot as plt


def get_bar_chart(data: dict[int, int], bins, title: str) -> None:
    # Get the frequency of each slice of bin
    bin_size = max(data.keys())//bins + 10
    # aproximate the value of the bin to the nearest integer divisible by 10
    bin_size = int(bin_size/10)*10
    # Create a list of bins
    bin_list = [bin for bin in range(0, max(data.keys()), bin_size)]
    # Create a list of frequencies
    frequency_list = [0 for _ in range(len(bin_list))]
    # Fill the frequency list
    for key, value in data.items():
        index = key//bin_size
        try:
            frequency_list[index] += value
        except IndexError:
            frequency_list[-1] += value
    # Create a dataframe
    frequency_list = [value/sum(frequency_list)
                      for value in frequency_list]  # to be in percentage
    df = DataFrame({"bins": bin_list, "frequency": frequency_list})
    # Create the bar chart
    plt.bar(df["bins"], df["frequency"], width=bin_size,
            color="#3686C9", edgecolor='black', align='edge')
    plt.title(title)
    # plt.xlabel("Genomic position")
    # plt.ylabel("Frequency of SNPs")
    # plt.text((21555+ 266) // 2, -0.025, 'ORF1AB')
    # plt.text((21563+25384) // 2 , -0.025, 'S')
    # plt.text((25393+26220) // 2 , -0.045, 'ORF3a', rotation=90)
    # plt.text((26245+26472) // 2 , -0.045, 'E', rotation=90)
    # plt.text((26523+27191) // 2 , -0.045, 'M', rotation=90)
    # plt.text((27202+27387) // 2 , -0.045, 'ORF6', rotation=90)
    # plt.text((27394+27759) // 2 , -0.045, 'ORF7a', rotation=90)
    # plt.text((27894+28259) // 2 , -0.045, 'ORF8', rotation=90)
    # plt.text((28274+29533) // 2 , -0.045, 'N', rotation=90)
    # plt.text((29558+29674) // 2 , -0.045, 'ORF10', rotation=90)
    # plt.xlim(21563, 25384)
    plt.ylim(0, 1)
    plt.savefig(
        f"{title}.png")
    plt.show()


def get_dictionary(filename: str) -> dict[int, int]:
    with open(filename, "r", encoding="utf8") as handle:
        reader = csv.reader(handle)
        dic = {}
        for row in reader:
            key, value = row[0].split("\t")
            dic[int(key)] = int(value)
    return dic


def rollete_choice(percentage_list: list[float]) -> int:
    random_number = random()
    for idx, percentage in enumerate(percentage_list):
        if random_number <= percentage:
            return idx
    return len(percentage_list) - 1


def rollete(data: dict[int, int], roll_times: int):
    """
    The rollete function takes a dictionary of values and their counts,
    and returns a list of the chosen values. The function uses the rollete method to
    choose which value to pick next based on its percentage chance of being picked.
    The function will continue until it has rolled all the desired values.

    :param data: dict[int: Represent the position of the value in the rollete,
    and int is used to represent how many times that value has been chosen
    :param int]: Specify the type of data that is going to be returned
    :param roll_times: int: Determine how many times the function will roll
    :return: A list of values that were chosen
    :doc-author: Trelent
    """
    chosen_values: list[int] = []
    percentage_list: list[float] = []
    positions_list: list[int] = []
    current_percentage: float = 0
    total_count = sum(data.values())
    for _ in range(roll_times):
        # print(_)
        for key, value in data.items():
            current_percentage += value / total_count
            percentage_list.append(current_percentage)
            positions_list.append(key)
        chosen_index = rollete_choice(percentage_list)
        # print(chosen_index, positions_list)
        chosen_key = positions_list[chosen_index]
        chosen_values.append(chosen_key)
        total_count -= data[chosen_key]
        del data[chosen_key]
        percentage_list = []
        positions_list = []
        current_percentage = 0
    return chosen_values


def generate_rollete(positions_file, roll_times: int, positions_dict: dict[int, int]):
    """
    The generate_rollete function takes a file of positions and rolls the rollete that many times.
    The function returns a list of strings, each string is one roll.

    :param positions_file: Get the dictionary of words from the file
    :param roll_times: int: Specify how many times the rollete will be rolled
    :return: A list of strings
    :doc-author: Trelent
    """

    return rollete(positions_dict, roll_times)


def generate_new_reference(old_reference_file, positions_change):
    with open(old_reference_file, "r", encoding="UTF8") as handler:
        old_ref = SeqIO.parse(handler, "fasta")
        old_ref = list(old_ref)[0].seq
        old_ref = str(old_ref)
        old_ref = [*old_ref]
    for position in positions_change:
        print(position)
        if old_ref[position] == "A":
            old_ref[position] = "T"
        elif old_ref[position] == "T":
            old_ref[position] = "A"
        elif old_ref[position] == "C":
            old_ref[position] = "G"
        elif old_ref[position] == "G":
            old_ref[position] = "C"
    new_reference = "".join(old_ref)
    return new_reference


def main():
    positions_file = "../montly_information/positions_2022_11_2022_12.csv"
    number_of_snps = 50
    positions_dict = get_dictionary(positions_file)
    total_positions = Counter({})
    # for _ in range(100000):
    #     print(_)
    #     total_positions += Counter(rollete(positions_dict.copy(),
    #                                number_of_snps))
    # get_bar_chart(dict(total_positions), 50,
    #               "Artifical_Data_50bins_100000_samples")
    get_bar_chart(positions_dict, 50,
                  "Pos2022_2023")



if __name__ == "__main__":
    main()
