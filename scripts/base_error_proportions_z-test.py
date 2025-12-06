#!/usr/bin/env python3

from statsmodels.stats.proportion import proportions_ztest

# Data for each error type category between dZTP and dATP on the Revio platform
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

# Assuming `nobs_dZTP` and `nobs_dATP` represent the total number of observations in each condition
# (e.g., total number of reads or total opportunities for mutations to occur).
# These should be specific to your experiment.
nobs_dZTP = sum([172297, 159515, 49844, 98421, 155037, 145146, 154175, 95559])  # Total observations for dZTP
nobs_dATP = sum([182272, 324824, 84115, 176128, 174652, 264718, 178627, 128216])  # Total observations for dATP

for category, counts in revio.items():
    count = counts  # count contains [count_dZTP, count_dATP]
    nobs = [nobs_dZTP, nobs_dATP]
    stat, p_value = proportions_ztest(count, nobs)
    print(f"revio\t{category}\t{stat}\t{p_value}")

# Data for each error type category between dZTP and dATP on Sequel 2 platform
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

# Assuming `nobs_dZTP_seq2` and `nobs_dATP_seq2` represent the total number of observations in each condition
# (e.g., total number of reads or total opportunities for mutations to occur).
nobs_dZTP_seq2 = sum([290683, 229210, 119299, 186184, 570135, 296293, 244830, 189755])  # Total observations for dZTP
nobs_dATP_seq2 = sum([226663, 511250, 388428, 654645, 215526, 1060852, 470236, 529041])  # Total observations for dATP

for category, counts in seq2.items():
    count = counts  # count contains [count_dZTP, count_dATP]
    nobs = [nobs_dZTP_seq2, nobs_dATP_seq2]
    stat, p_value = proportions_ztest(count, nobs)
    print(f"Sequel 2\t{category}\t{stat}\t{p_value}")