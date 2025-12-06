#!/usr/bin/env python3

import numpy as np
from math import sqrt, asin

def cohen_h(p1, p2):
    return 2 * (asin(sqrt(p1)) - asin(sqrt(p2)))

revio = {
    'DEL': [172297, 182272],
    'G->A': [159515, 324824],
    'G->C': [49844, 84115],
    'G->T': [98421, 176128],
    'INS': [155037, 174652],
    'T->A': [145146, 264718],
    'T->C': [154175, 178627],
    'T->G': [95559, 128216]
}

total_dZTP = sum([c[0] for c in revio.values()])
total_dATP = sum([c[1] for c in revio.values()])

for category, counts in revio.items():
    p1 = counts[0] / (total_dZTP/1.695)
    p2 = counts[1] / total_dATP
    h = cohen_h(p1, p2)
    print(f"revio\t{category}\t{h}")

seq2 = {
    'DEL': [290683, 226663],
    'G->A': [229210, 511250],
    'G->C': [119299, 388428],
    'G->T': [186184, 654645],
    'INS': [570135, 215526],
    'T->A': [296293, 1060852],
    'T->C': [244830, 470236],
    'T->G': [189755, 529041]
}


total_dZTP = sum([c[0] for c in seq2.values()])
total_dATP = sum([c[1] for c in seq2.values()])

for category, counts in seq2.items():
    p1 = counts[0] / (total_dZTP/1.615)
    p2 = counts[1] / total_dATP
    h = cohen_h(p1, p2)
    print(f"seq2\t{category}\t{h}")