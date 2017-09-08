# Sequence-Order Frequency Matrix
## Introduction
Protein remote homology detection and fold recognition are two critical tasks for the studies of protein structure and function. Currently, the profile-based methods achieve the state-of-the-art performance in these fields. However, the widely used sequence profiles, like Position-Specific Frequency Matrix (PSFM) and Position-Specific Scoring Matrix (PSSM), ignore the sequence-order effects along protein sequence. In this study, we propose a novel profile, called Sequence-Order Frequency Matrix (SOFM), which can extract the amino acid patterns with sequence-order information from Multiple Sequence Alignment (MSA). The SOFM is applied to protein remote homology detection and fold recognition by combining with two profile vectorizing approaches Top-n-grams and Smith-Waterman algorithm, based on which two computational predictors called SOFM-Top and SOFM-SW are proposed, respectively. Experimental results show SOFM can achieve more information content than traditional profiles and these two predictors outperform other related state-of-the-art methods, indicating that SOFM is a more effective protein profile than traditional profiles. Because profiles have been widely used in constructing computational predictors for protein studies, SOFM will have many potential applications.  

If you use this library, please cite the following papaers:  
1. Junjie Chen, Mingyue Guo, Xiaolong Wang, BinLiu*. Protein remote homology detection and fold recognition based on Sequence-Order Frequency Matrix[J]. *IEEE/ACM Transactions on Computational Biology and Bioinformatics*, 2017.   
2. Junjie Chen, Mingyue Guo, Xiaolong Wang, BinLiu*. SOFM-Top: Protein Remote Homology Detection and Fold Recognition Based on Sequence-Order Frequency Matrix[C]//International Conference on Intelligent Computing. Springer, Cham, 2017: 469-480.

## Install
### 1. Main dependence  
* [python2 >= 2.7.6](https://www.python.org/)
* [ncbi-blast >= 2.2.30+](https://blast.ncbi.nlm.nih.gov/Blast.cgi)
* [Numpy >= 1.3.1](http://www.numpy.org/)
* [BioPython >= 1.69](http://biopython.org/)

### 2. Configure  
You have to configure the paths of ncbi-blast and blast database in the PATH.conf
```Bash  
# ncbi blast verison ncbi-blast-2.2.30+
[PSIBLAST]
NCBI_BLAST_BIN = /path/to/ncbi-blast-2.2.30+/bin/
BLAST_DATABASE = /path/tp/BLAST_DB/db_name
```  
The command of PSIBLAST used in SOFM is   
```  
{NCBI_BLAST_BIN}/psiblast -query input -db {BLAST_DATABASE} -out xmlfile -outfmt 5 -num_iterations 3 -evalue 0.001 -num_threads 5
```  

### 3. Test  
You can test it by running an example:  
```Bash  
python SOFM.py -i example/d119l__.fasta -k 3
```

## Usage
```
python SOFM.py -h  

-------------------
usage: SOFM.py [-h] -i INPUT [-o OUTPUT] [-k KMER] [--version]

Sequence-Order Frequency Matrix (SOFM) is a novel protein profile, which can
achieve more information content than traditional profiles.

optional arguments:
  -h, --help                    show this help message and exit
  -i INPUT, --input INPUT       input a single protein sequence in FASTA format
  -o OUTPUT, --output OUTPUT    output SOFM file, default filename: {input}.sofm{k}
  -k KMER, --kmer KMER          length of substrings, default value: k=3
  --version                     show program's version number and exit 
```  
   

