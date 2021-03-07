import statistics as stats
import sys
import numpy as np

def main():
    len_array = []
    for line in sys.stdin:
        attributes = line.strip().split()
        start = int(attributes[1]) 
        end = int(attributes[2])
        len_array.append(end - start)
    if len(len_array) == 0:
        print("None")
        return
    len_array.sort(reverse=True)
    total_len = sum(len_array)
    n_blocks = len(len_array)
    running_sum = 0
    nx_blocks = []
    nx_i = 0
    for i in range(n_blocks):
        running_sum += len_array[i]
        if running_sum >= (total_len * nx_i * 10 / 100):
            nx_i += 1
            nx_blocks.append(len_array[i])
    median_blocks = stats.median(len_array)
    mean_blocks = stats.mean(len_array)
    print("Stats:\n\t#: {}\n\tTotal Length: {}\n\tN50: {}\n\tMedian: {}\n\tMean: {:.2f}\n\tMax: {}\n\tMin: {}".format(n_blocks, total_len, nx_blocks[5], median_blocks, mean_blocks, len_array[0], len_array[-1]))
    print("\tNx:")
    for i in range(11):
        print("\t\t{}%: {:d}".format(i * 10,nx_blocks[i]))
    print("\tPercentiles:")
    percentiles = np.percentile(len_array,q=np.arange(0,110,10))
    for i in range(11):
        print("\t\t{}%: {:.1f}".format(i * 10,percentiles[i]))

main()
