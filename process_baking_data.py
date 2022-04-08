import csv
import numpy as np
import os
import re


def process_data(folder="Solution_baking"):
    points = []
    data = []

    # extract high temp data
    strings = []
    for f in os.listdir(folder):
        strings.append(os.path.join("/" + folder, f))

    for s in strings:
        match_number = re.compile('-?\ *[0-9]+\.?[0-9]*(?:[Ee]\ *-?\ *[0-9]+)?')
        e = re.findall(match_number, s)
        points.append(float(e[0])*10**float(e[1]))
        data.append({})
        data[-1]["T_baking"] = points[-1]
        filename = folder + "/T={:.2e}/derived_quantities.csv".format(points[-1])

        t = []
        inventory = []
        with open(filename, 'r') as csvfile:
            plots = csv.reader(csvfile, delimiter=',')
            next(plots)
            for row in plots:
                t.append(float(row[0]))
                inventory.append(
                    2*(float(row[-1]) +
                        float(row[-2]) +
                        float(row[-3])))

        data[-1]["t"] = t
        data[-1]["inventory"] = inventory
    return data


if __name__ == "__main__":
    process_data()
