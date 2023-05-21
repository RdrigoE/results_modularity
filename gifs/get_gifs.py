#!/usr/bin/env python

import matplotlib.pyplot as plt
import imageio
from pandas import DataFrame
import csv
import streamlit as st
from pathlib import Path


def get_positions(positions: tuple[int, int]):
    origin_folder = Path("montly_information")
    for file in origin_folder.iterdir():
        if file.is_file and file.suffix == ".csv":
            with open(file, "r") as f:
                lines = f.readlines()
                save_positions = []
                for line in lines:
                    position, _ = line.split("\t")
                    position = int(position)
                    if positions[1] >= position >= positions[0]:
                        save_positions.append(line)
                with open(f"gifs/positions/{file.name}", "w") as f:
                    f.writelines(save_positions)


months = [
    (12, 2019, 4, 2020),
    (5, 2020, 7, 2020),
    (8, 2020, 10, 2020),
    (11, 2020, 1, 2021),
    (2, 2021, 4, 2021),
    (5, 2021, 7, 2021),
    (8, 2021, 10, 2021),
    (11, 2021, 1, 2022),
    (2, 2022, 4, 2022),
    (5, 2022, 7, 2022),
    (8, 2022, 10, 2022),
    (11, 2022, 12, 2022),
]


def get_dictionary(filename: str) -> dict[int, int]:
    with open(filename, "r", encoding="utf8") as handle:
        reader = csv.reader(handle)
        dic = {}
        for row in reader:
            key, value = row[0].split("\t")
            dic[int(key)] = int(value)

    return dic


# @st.cache_data
def get_bar_chart(data: dict[int, int], bins: int, title: str, positions: tuple[int, int], month, lim) -> None:
    # Get the frequency of each slice of bin
    # bin_size = max(data.keys())//bins*10 + 10
    # aproximate the value of the bin to the nearest integer divisible by 10
    amplitude = positions[1] - positions[0]
    print(data.keys())
    # print(amplitude, bin_size)
    # bin_size = int(bin_size/amplitude)
    bin_size = int(amplitude * bins)
    print("New", bin_size)
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
    plt.xlim(positions[0], positions[1])
    plt.ylim(lim[0], lim[1])
    plt.savefig(
        f"gifs/images/barplot_{month[1]}_{month[0]}_{month[3]}_{month[2]}.png")
    plt.clf()


def create_gif(gene: str, bins: int, positions: tuple[int, int], duration: float, ylim: tuple[float, float]):

    get_positions(positions)
    all_info = {}
    for month in months:
        dic = get_dictionary(
            f"gifs/positions/positions_{month[1]}_{month[0]}_{month[3]}_{month[2]}.csv")
        all_info[f"{month[1]}_{month[0]}_{month[3]}_{month[2]}"] = dic
    all_info.keys()

    for month in months:
        month_string = f"{month[1]}_{month[0]}_{month[3]}_{month[2]}"
        get_bar_chart(all_info[month_string], bins,
                      f"Percentage of SNPs from {month[1]}-{month[0]} to {month[3]}-{month[2]} at {gene}", positions, month, ylim)
    images = []
    filenames = []
    for month in months:
        filenames.append(
            f"gifs/images/barplot_{month[1]}_{month[0]}_{month[3]}_{month[2]}.png")

    for filename in filenames:
        images.append(imageio.imread(filename))
    print(duration) 

    durations = [duration*500] + [duration*100 ] * (len(months) - 1)
    imageio.mimsave(f'gifs/images/{gene}_{bins}_{duration}.gif',
                    images, duration=durations, loop=0)


def get_relative_frequency(data: dict[int, int], positionStart: int, positionEnd: int, all_info: dict[str, dict[int, int]]):
    inside_position = 0
    outside_position = 0
    for key, value in data.items():
        if positionStart <= key <= positionEnd:
            inside_position += value
        else:
            outside_position += value
    total = inside_position + outside_position
    return (inside_position / total, outside_position / total)


# In[20]:
#@st.cache_data
def get_relative_plot(positions: tuple[int, int], gene: str):
    all_info: dict[str, dict[int, int]] = {}
    for month in months:
        dic = get_dictionary(
            f"montly_information/positions_{month[1]}_{month[0]}_{month[3]}_{month[2]}.csv")
        all_info[f"{month[1]}_{month[0]}_{month[3]}_{month[2]}"] = dic

    inside_l, outside_l = [], []
    for _, value in all_info.items():
        inside, outside = get_relative_frequency(
            value, positions[0], positions[1], all_info)
        inside_l.append(inside)
        outside_l.append(outside)
    plt.title(
        f"Relative frequency of SNPs in the region of {positions[0]}-{positions[1]} ({gene})")
    plt.ylabel("Relative SNPS frequency")
    plt.xlabel("Months")
    plt.plot(inside_l, label="Inside", marker="o")
    plt.grid(True)
    plt.xticks(labels=all_info.keys(), rotation=90,
               ticks=range(len(all_info.keys())))
    plt.savefig(f"gifs/images/relative_frequency_{gene}.png",
                dpi=300, bbox_inches='tight')
    plt.clf()


if __name__ == "__main__":
    positions = (20000, 24000)
    bins = 100
    get_positions(positions)
    print("hello")
    create_gif("GENE_S", bins, positions, 0.5, (0.0, 0.5))
