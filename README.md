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

### Step 4：Run the pipeline
Open the command line window (for example, terminal in Macbook), cd to the working directory and run the analysis pipeline.

cd path_to_your_working_directory

python sgRNA_desgin_main.py configure.txt

We also post a toy example together with the scripts and the example_configure.txt has been edit to make it compatible. For this test, cd to the working directory, type in: 

python sgRNA_desgin_main.py example_configure.txt

Check [here](./image/successful_running.png) for the output during a successful running of the abovementioned test.

For a typical Macbook (for example, 2.6 GHz processor and 8 GB memory), the example test can be finalized within 30 minutes. The rate-limiting step is off-target site identification across the genome. For a typical Macbook, we expect a processing speed of ~ 50 sgRNAs designed per minute. For a genome-scale library with 50 k members (10 sgRNAs per gene assuming 5 k genes encoded by a genome), the laptop thus needs roughly 20 hours to finalize the design process.

## Output description
The output files will be organized in the subdirectory whose name is specified by the 'prefix' option in configure file under the working directory (prefiex/). We term this subdirectory 'result directory' thereafter.

You can find many files under the result directory.

[your result directory after running the test](./image/resultdir_after_example_running.png)

Below is the description. For the mathematical processing, see our paper. **All .csv flat files use tab as delimiter unless mentioned**

**prefix.N20.fasta.txt**: sgRNA library sequences (N20) in .fasta format, including negative control sgRNAs if it is specified in the configure file. **This file, after adding designed flanking nucleotides, can be subjected to microarray based oligomer synthesis** to prepare the sgRNA library. 

**prefix.N20NGG.fasta.txt**: The reverse complementary of target region of each sgRNA including the PAM site in .fasta format.

**prefix.sgRNA_statistics.txt**: position information within relevant gene coding region and the GC content of each designed sgRNA in .csv format.

sgRNAID|sgRNA_position_in_gene |GCcontent
-------|-----------------------|---------
insI1b4284_32|0.028|39.13
...|...|...

**prefix.cluster.txt**: Clusters of genes with homologs (see Step 2, "Coping with multiple copy issue" section). The file has no header line and uses tab as delimiter. It is consisted of two columns: cluster name and the genes in each cluster. Genes in one cluster are separated by comma.

araFb1901|araFb1901
-------|--------------------
yjhXb4566|yjhXb4566
tufAb3339|tufAb3339,tufBb3980
...|...

**prefix.gene_statistics.txt**: file in .csv format with one header line contains the length and the designed sgRNA numbers of each gene.

gene_name|gene_length|sgRNA_number_in_gene
---------|-----------|--------------------
ydcCb1460|1137|0
insI1b4284|1152|10
pinRb1374|591|9

**prefix.fasta.txt**: target gene sequences in .fasta format. It is similar to the input gene fasta file. For genome-wide sgRNA library design, we use the combination of gene and synonym (.ptt or .rnt annotation file) to rename each gene. Hence, we give the refined .fasta file with new gene names for some convenience in following usage.

**N20_library.csv**: the sgRNA library file is at .csv formate **containing one header line**, in which there are three columns in order of id, sequence and gene respectively. **Use comma as delimiter**. If negative control (NC) sgRNAs are within this synthetic library (specified in configure file), name them NCx and assign '0' at 'gene' column of these sgRNAs. This file is used as input for the data processing subpackage after pooled screening experiment and NGS.

id|sequence|gene
--|--------|----
sgRNA1|ATCCCCCCCCCCGGGGG|recA
NC1|TGTGTGTGTGTGTGTGTGTG|0
...|...|...

**sgRNA_position.txt**: Flat file of sgRNA position (relative location of sgRNA in the coding region) information in gene **without header line**. The file contains three columns in order of gene name, sgRNAid and the relative position of sgRNA in the gene. Actually, you can also find this file as output of our library design subpackage. This file is also used as input for the data processing subpackage after pooled screening experiment and NGS.

rsmE|rsmE_9|0.012
----|------|-----
rsmE|rsmE_10|0.014
0|NC1|0
...|...|...

Two .png figures use histogram to summarize the basic information of [number of sgRNA per gene](./image/The_distrution_of_sgRNA_number_per_gene.png) and [position of sgRNAs in the coding region of relevant genes](./image/The_distrution_of_sgRNA_position_per_gene.png).
