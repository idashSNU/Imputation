# Privacy-preserving Genotype Imputation based on Homomorphic Encryption

### Introduction

This is a solution from the SNU team for privacy-preserving genotype imputation based on Homomorphic Encryption (HE). Basically it consists of four directories: `/ModHEaaN`, `/data_origin`, `/plain` and `/encrypted`.
`/ModHEaaN` includes an HE library [ModHEaaN](https://github.com/idashSNU/Imputation/tree/master/ModHEaaN) which implements an approximate HE scheme [CKKS](https://eprint.iacr.org/2016/421.pdf). In `/plain` directory, there are some python files that generates model files in name `New_gen_model_W_*.py`, and these models are saved in `/encrypted` directory. The original data should be downloaded from external link.
You can find our main code for HE-based genotype imputation in `/encrypted`. 

## Download Original Data

At first, we need two groups of large data. Please follow steps below:
1. Download original data files from [here](https://drive.google.com/drive/folders/1EVFLogAoqAajHxCBlen4vzy2Y0JbkpbU?usp=sharing) and save in folder `/data_origin`. 
1. Download modified data files from [here](https://drive.google.com/drive/folders/15JNx48B-dUDoIr1eVNqegj2fmIMB9K9Y?usp=sharing) and save in folder `/plain/data_mod`.


## how to Run

The following is the way to build and run our solution:

### Build plain model
We generate models for several datasets represented by "population". For each dataset, one can choose different models determined by "window_size", which denotes the number of adjacent tag SNPs for each target SNP. Experiments based on the given dataset shows that the choice "window_size = 40" provides the best accuracy of genotype imputation.
Note that the larger "window_size" implies the higher computational cost. We also provide the multi-processing option for the acceleration so that one can dynamically choose the number of processes based on his/her computer environment.

To generate the imputation models, command the following:
```bash
$ cd ./plain
$ python New_gen_model_W.py -p <population> -w <window_size> -n <number_of_processes>
```
* population: Total, EUR, AMR, AFR 
* window_size: 8, 16, 24, 32, 40, 48, 56, 64, 72 
* number_of_processes: 1, 2, 4, 8, 16 

The generated plain models will be saved in the directory `/encrypted/<population>_DNNmodels/DNNmodels_<window_size>_c`.

### Build ModHEaaN
For the encryption of test data, we use the [ModHEaaN](https://github.com/idashSNU/Imputation/tree/master/ModHEaaN) library, which is a light-version implementation of the approximate HE scheme [CKKS](https://eprint.iacr.org/2016/421.pdf). Contrary to the original implementation of the [CKKS](https://eprint.iacr.org/2016/421.pdf) scheme [HEAAN](https://github.com/snucrypto/HEAAN), our [ModHEaaN](https://github.com/idashSNU/Imputation/tree/master/ModHEaaN) library does not have any dependency on multi-precision libraries GMP and NTL, and only supports homomorphic addition and constant multiplication (hence bootstrapping disabled).

To build [ModHEaaN](https://github.com/idashSNU/Imputation/tree/master/ModHEaaN), command the following:
```bash
$ ./ModHEaaN/heaan
$ cmake CMakeLists.txt
$ make all
```

### Build HE-imputation executable by
```bash
$ ./encrypted/impute_dnn
$ cmake CMakeLists.txt
$ make all
```
The executable file is generated as `enc_impute` in the directory `./encrypted/impute_dnn`. 

### Run
Choose two parameters, `window_size` and `number_of_targetSNP`. `window_size` is the number of adjacent tag SNPs for each target SNP, and `number_of_targetSNP` is literally the number of target SNPs. The choices of `window_size` and `number_of_targetSNP` are among (8 / 16 / 24 / ... / 72) and (20 / 40 / 80), respectively. 
```bash
$ cd ./encrypted/impute_dnn
$ ./enc_impute <window_size> <number_of_targetSNP>
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
