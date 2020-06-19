# SNU Team IDASH2019

### Introduction

This is a solution from the SNU team for HE-based genotype imputation. It consists of four directories: `/ModHEaaN`, `/data_origin`, `/plain` and `/encrypted`.
`/ModHEaaN` includes a library which implements the CKKS homomorphic encryption scheme (without RNS). In `/plain` directory, there are some python files that generates model files in name `New_gen_model_W_*.py`, and these models are saved in `/encrypted` directory. The original data should be downloaded from external link.
You can find our main code for HE-based genotype imputation in `/encrypted`. 

## Download Original Data

At first, we need two groups of large data. Please follow steps below:
1. Download original data files from [this link](https://drive.google.com/drive/folders1EVFLogAoqAajHxCBlen4vzy2Y0JbkpbU?usp=sharing) and save in folder `/data_origin`. 
1. Download modified data files from [this link](https://drive.google.com/drive/folders/15JNx48B-dUDoIr1eVNqegj2fmIMB9K9Y?usp=sharing) and save in folder `/plain/data_mod`.


## how to Run

The following is the way to build and run our solution:

### Build plain model
Plain model here

### Build ModHEaaN
ModHeean Here

### Build HE-imputation executable by
```bash
$./encrypted/impute_dnn
$ cmake CMakeLists.txt
$ make all
```
### Run
Choose two parameters, `window_size` and ``. `window_size` is the number of nearby (known) SNPs used for each target SNP, and `number_of_targetSNP` is literally the number of target SNPs. The choices of `window_size` and `number_of_targetSNP` are among (8 / 16 / 24 / ... / 72) and (20 / 40 / 80), respectively. The executable with command line arguments by
```bash
$ ./enc_impute window_size number_of_targetSNP
```
For instance, The argument below runs imputation of 80k target SNPS, with window size 40.
```bash
$ ./enc_impute 40 80
```

* Note : The running time of (encrypted) imputation grows with `window_size`, but the quality of imputation (accuracy) does not. For the sample data, a moderate window size (about 40) shows the best result.
* WARN: mode `populations` is only compatible with '80'. The command for poluation mod is as follow.
```bash
$ ./enc_impute populations 80
```

### Measure MicroAUC
If you succeed to run our enc_impute solution (Step 2), then the score output file will be saved in `/encrypted/impute_dnn` in name `score_window(window_size)__(number_of_targetSNP)k.csv`. Note that the real genotypes of test data is saved in `/plain` named as `real_number_of_targetSNP)k.csv`.

e.g. After running `$ ./enc_impute 40 80` in Step 2, then
```bash
$ cd ./plain
$ python3 evaluation.py -i ../encrypted/impute_dnn/score_window40_80k.csv -t real_80k.csv -o output.png
```

e.g. After running `$ ./enc_impute populations 80` in Step 2, then
```bash
$ cd ./plain
$ python3 evaluation.py -i ../encrypted/impute_dnn/score_EUR_80k.csv -t real_EUR_80k.csv -o output_EUR.png
$ python3 evaluation.py -i ../encrypted/impute_dnn/score_AMR_80k.csv -t real_AMR_80k.csv -o output_AMR.png
$ python3 evaluation.py -i ../encrypted/impute_dnn/score_AFR_80k.csv -t real_AFR_80k.csv -o output_AFR.png
```
