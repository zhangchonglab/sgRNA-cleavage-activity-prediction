# sgRNA-cleavage-activity-prediction

## What is this?
This python script collection is used to predict the activity scores of sgRNAs (CRISPR/Cas9). This program is especially useful for biologists and bioengineerers interested in using CRISPR/Cas9 genome editing technology in bacteria, because the sgRNA activity dataset used to train this model is experimentally determined in model bacteria *E. coli*, which is the first such dataset obtained in bacterial host to the best of our knowledge. The paper comprehensively describing this algorithm and relevant experiments is currently under review, we are going to post it once its publication. Please cite the paper if this program is useful to your work.

## General description of the algorithm
The utility of this software is fairly simple. Given a fasta file containing the target DNA sequences (N4N20NGGN3), the program will predict activity scores for each one, with output in .csv format.

## How to use it?
### Step 1：Installation
1. Install Python version 2.7
2. Install biopython version 1.66 or above
3. Install matplotlib version 2.0.2 or above
4. Install numpy version 1.13.1 or above
5. Install scikit-learn version **0.19.0 (This is important)**
6. Install pandas version 0.18.1 or above

### Step 2：Prepare the necessary files.
All these files (or subdirectories) should be organized under a common working directory together with the all .py scripts.
Please check the example files post at GitHub, which are described as below.
