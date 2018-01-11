# sgRNA-cleavage-activity-prediction

## What is this?
This python script collection is used to predict the activity scores of sgRNAs (CRISPR/Cas9). This program is especially useful for biologists and bioengineerers interested in using CRISPR/Cas9 genome editing technology in bacteria, because the sgRNA activity dataset used to train this model is experimentally determined in model bacteria *E. coli*, which is the first such dataset obtained in bacterial host to the best of our knowledge. The paper comprehensively describing this algorithm and relevant experiments is currently under review, we are going to post it once its publication. Please cite the paper if this program is useful to your work.

## General description of the algorithm
The utility of this software is fairly simple. Given a fasta file containing the target DNA sequences (N<sub>4</sub>N<sub>20</sub>NGGN<sub>3</sub>), the program will predict activity scores for each one, with output in .csv format.

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

#### File 1: fasta file of DNA target sequences (see example_sgRNA.fasta)
The fasta file contains the target DNA sequences with N<sub>4</sub>N<sub>20</sub>NGGN<sub>3</sub> format (N = A, T, C, G).

#### File 2: configure file (see example_configure.txt)
The configure file is used to set all the necessary parameters and tell the program where to find necessary files. This file is in a two-column format using colon (:) as delimiter. Each line starts with one word (name of one parameter) separated with the following (setting of this parameter) by a colon delimiter. We describe each parameter as below.

**targetFasta**: The location of the abovementioned fasta file of DNA target sequences.

**model**: This software provides two options to choose a model to predict sgRNA activity, "Cas9" and "eSpCas9". For the details, please read our paper. Briefly, "Cas9" refers to the wild type Cas9 nuclease, which is currently most widely used in CRISPR/Cas genome editing. By contrast, "eSpCas9" is a off-target-reducing mutant of "Cas9" with three amino acid mutations.  

**normalization**: Because we use Z score of the original dataset to train the model, hence the output for each sgRNA is a real number typically located between (0, 50). To normalize the result, we provide here an option to linearly project the result onto (0, 1). "Yes" or "No" should be specified here. 

**prefix**: prefix used for naming of all output files, keep it simple. For example, "myprediction" is fine.

Below is **an example configure file with default parameters**.

parameter|value
---------|-----
[configdesign]
targetFasta:|example_target.fasta
model:|Cas9
normalization:|Yes
prefix:|myprediction
