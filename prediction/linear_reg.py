#!/usr/bin/env python
# coding: utf-8

import matplotlib.pyplot as plt
import numpy as np
from sklearn.linear_model import LinearRegression
from pandas import DataFrame


# In[932]:


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

data = []
names = []
dfs = []
for month in months:
    with open(f"montly_information/number_of_snps_per_sample_{month[1]}_{month[0]}_{month[3]}_{month[2]}.txt", "r", encoding="utf8") as handle:
        data.append(sorted(list(map(lambda x: x.strip(), handle.readlines()))))
        names.append(f"{month[1]}_{month[0]}_{month[3]}_{month[2]}")
        dfs.append(
            DataFrame(data[-1], columns=[f"{month[1]}_{month[0]}_{month[3]}_{month[2]}"]))

final_df = dfs[0]
for df in dfs[1:]:
    final_df = final_df.join(df)

median_list = []
for col in final_df:
    lists = []
    for x in final_df[col]:
        if isinstance(x, str) and isinstance(int(x), int):
            lists.append(int(x))
    lists = sorted(lists)
    if len(lists):
        median_list.append(lists[int(len(lists)//2)])

x = np.array([1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12]).reshape((-1, 1))
y = np.array(median_list)
model = LinearRegression().fit(x, y)
new_model = LinearRegression().fit(x, y.reshape((-1, 1)))
y_pred = model.predict(x)
predict_range = range(1, 81, 4)
b = list(predict_range)
new_x = np.array([predict_range]).reshape((-1, 1))
new_y = model.predict(new_x)
a = [x for x in range(1, 13)]
# for x in new_y:
#     median_list.append(x)
#
list_years = []
year = 2020
for times in predict_range:
    list_years.append(year)
    year += 1
list_years = list(map(str, list_years))
a = [2020, 2020, 2020, 2020, 2021, 2021, 2021, 2021, 2022, 2022, 2022, 2023]
a = list(map(str, a))

plt.plot(list_years, new_y)
plt.plot(a, median_list)
plt.legend(["predicted", "real"])
plt.ylabel("Number of SNPs in SARS-CoV-2 Genome")
plt.title("SNPs frequency prediction in SARS-CoV-2")
plt.xlabel("years")
plt.xticks(rotation=90)
plt.savefig("prediction.png")
