import os
import sys

count_dict = {}
with open("current_alignment_index.txt") as fh:
    for line in fh:
        line = line.strip()
        gca = line.split("/")[2]
        if gca not in count_dict:
            count_dict[gca] = line
        else:
            count_dict[gca] = "2"


for gca in count_dict:
    if count_dict[gca] != "2":
        print(count_dict[gca])