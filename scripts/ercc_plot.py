import logging
import sys
import math
import numpy as np
import matplotlib.pyplot as plt

def ercc_plot(args):
    ercc_counts = {}
    dirpath = args[0]

    # sort the ercc_expected_concentration file just in case, such that our ercc_coverage.txt file can be compared simultaneously
    with open(args[1]) as f:
        for line in f:
            if line == "\n":
                continue
            if line == "":
                break
            curr = line.strip().split("\t")
            if len(curr) < 2:
                continue
            if float(curr[1]) == 0:
                continue
            ercc_counts[curr[0]] = [math.log(float(curr[1]), 10)]

    # now do the same but for the actual read counts
    with open(dirpath + "/ercc_coverage.txt") as f:
        for line in f:
            if line == "\n":
                continue
            if line == "":
                break
            curr = line.strip().split("\t")
            if len(curr) < 2:
                continue
            if float(curr[1]) == 0:
                continue
            if curr[0] in ercc_counts.keys():
                ercc_counts[curr[0]].append(math.log(float(curr[1]), 10))

    keys_to_remove = []

    for key in ercc_counts.keys():    
        if len(ercc_counts[key]) == 1:
            keys_to_remove.append(key)

    for key in keys_to_remove:
        del ercc_counts[key]

    # check if we have at least one value
    if len(ercc_counts) > 0:
        data = np.asarray(list(ercc_counts.values()))
        x_values = data[:,0]
        y_values = data[:,1]
        #plt.plot(data[:,0], data[:, 1])
        a, b = np.polyfit(x_values, y_values, 1)
        plt.scatter(x_values, y_values, color='red', edgecolor='black', linewidth=0.4)

        best_fit = [a*x_val+b for x_val in x_values]
        mean_y = y_values.mean()

        # step 3
        differences_regression_line = []
        for i in range(0, len(x_values)):
            differences_regression_line.append((a * x_values[i] + b) - y_values[i])

        regression_sum = 0
        for i in differences_regression_line:
            regression_sum += (i*i)

        # step 4
        differences_mean_line = [mean_y - y_val for y_val in y_values]
        mean_sum = 0
        for i in differences_mean_line:
            mean_sum += (i*i)

        # step 5
        r_squared = (mean_sum - regression_sum)/mean_sum

        plt.plot([min(x_values), max(x_values)], [min(best_fit), max(best_fit)], color='black', linestyle=':', linewidth=0.8)
        plt.text(min(x_values), max(y_values), "R² = " + str(r_squared))

        # we have got the above from: https://www.statology.org/line-of-best-fit-python/
        plt.xlabel("Log10 ERCC spike-in concentration")
        plt.ylabel("Log10 coverage per gene")
        plt.savefig(dirpath + '/ercc_plot.png', dpi=300) # dpi to control resolution

if __name__ == "__main__":
    logging.basicConfig(level=logging.INFO)
    logger = logging.getLogger(__name__)
    logger.info("Running ERCC plot script")
    args = sys.argv[1:]
    ercc_plot(args)