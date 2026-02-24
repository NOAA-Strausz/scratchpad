#!/usr/bin/env python


import pandas as pd
import matplotlib.pyplot as plt

names=["time", "O2conc", "O2sat"]

df = pd.read_csv("optode_682_zero_cleaned.csv", header=0, names=names, parse_dates=["time"])
df.set_index("time", inplace=True)
df.O2conc.plot(ylabel = "O2 Concentration", legend=True)
df.O2sat.plot(secondary_y=True, ylabel="O2 Percent Sat", legend=True, title="Optode zero bath response")
plt.show()
