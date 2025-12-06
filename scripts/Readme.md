# **Analysis Scripts for Z-Calling**

This directory contains the custom Python scripts used for the development, benchmarking, and evaluation of Z-Calling as described in the Supplementary Notes.

## **1\. Quality Value (QV) & Error Analysis**

Used in **Supplementary Note 2** to assess base quality and error distributions.

* **QV\_statistics\_and\_making\_graphs.py**  
  * **Description:** Performs Mann-Whitney U tests to compare Quality Value (QV) distributions between canonical A and Z conditions across different platforms and pass numbers (NP). It also adjusts p-values using the Benjamini-Hochberg procedure.  
* **VCF.stats.py**  
  * **Description:** Processes VCF files generated from CCS reads to count occurrences of specific base substitution types (e.g., G-\>A, T-\>C) and indels (DEL, INS), filtering out non-sense alleles.  
* **base\_error\_proportions\_z-test.py**  
  * **Description:** Performs a Z-test for proportions (proportions\_ztest) to statistically compare the frequency of specific base error types between dZTP and dATP amplicon datasets.  
* **base\_error\_proportions\_Cohen\_h.py**  
  * **Description:** Calculates Cohen’s *h* effect size to quantify the magnitude of difference in error proportions between conditions.

## **2\. Kinetic Signal Analysis**

Used in **Supplementary Note 3** to explore polymerase kinetics (IPD and PW).

* **extract\_kinetics\_signals.py**  
  * **Description:** Parses BAM files to extract Pulse Width (PW) and Inter-Pulse Duration (IPD) values for specific k-mer contexts surrounding target bases. Supports both specific motif targeting and global extraction.  
* **tsne\_pca\_analysis.py**  
  * **Description:** Performs dimensionality reduction (PCA and t-SNE) on extracted kinetic signals to visualize clustering patterns between A and Z bases. Includes data down-sampling for balanced comparisons.  
* **ZA\_signal\_ttest\_and\_medianratio.py**  
  * **Description:** Conducts two-tailed t-tests (with BH correction) and calculates median ratios for IPD and PW signals at each position within a k-mer context to identify significant kinetic differences.  
* **ZA\_signal\_feature\_importance.py**  
  * **Description:** Trains a Random Forest Classifier on kinetic features to calculate feature importance scores, identifying which signal positions (e.g., PW at \+1 position) are most critical for discrimination.  
* **kinetics\_stat.py**  
  * **Description:** Analyzes the relationship between the median PW/IPD ratio and false positive rates. Used to determine the filtering thresholds (top 1% and ratio \< 0.3) described in **Supplementary Note 5**.

## **3\. Model Training & Benchmarking**

Used in **Supplementary Note 3 & 4** for deep learning model development.

* **torch\_mlp\_shap\_gpu.py**  
  * **Description:** Trains a PyTorch MLP model and performs SHAP (SHapley Additive exPlanations) analysis to interpret feature importance, revealing the biological mechanism of Z-base detection.  
* **benchmark\_models\_gpu.py**  
  * **Description:** A comprehensive benchmarking script that trains and evaluates five different architectures (RF, MLP, CNN, BiLSTM, GCN) on identical datasets using dynamic early stopping to determine the optimal model architecture.  
* **mlp\_layer\_benchmark.py**  
  * **Description:** Systematically tests MLP architectures with varying depths (1 to 5 hidden layers) across different motifs to determine the optimal network depth (found to be 2 layers).
