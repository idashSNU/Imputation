# Privacy-preserving Genotype Imputation based on Homomorphic Encryption

This is a solution from the SNU team for privacy-preserving genotype imputation based on Homomorphic Encryption (HE). Basically it consists of four directories: `/ModHEaaN`, `/data_origin`, `/plain` and `/encrypted`. To run our HE-based genotype imputation solution, follow the instruction described below in order.


## Download Datasets

At first, we need to download four datasets "Total", "EUR", "AMR" and "AFR". Please follow the steps below:
1. Download original data files from [here](https://drive.google.com/drive/folders/1lLynN5dyh9nyqYIrWt1Wm5kUCoqZE0Vs?usp=sharing) and save in folder `/data_origin`. 
1. Download modified data files from [here](https://drive.google.com/drive/folders/1HfcYFgt0H3dFzdazVL0AjpQUwNVYuNjf?usp=sharing) and save in folder `/plain/data_mod`.
1. (Optional) Download modified data files for low from Minor Allele Frequency data from [here](https://drive.google.com/drive/folders/1ZprQwkqOS4CVMe4HhAo0-831f5j6UHD6?usp=sharing) and save in folder `/plain/data_mod_lowMAF`.

* Total: `/data_origin/*.txt`, `/plain/data_mod/Total_mod/*.txt`
* EUR: `/data_origin/*_EUR.txt`, `/plain/data_mod/EUR_mod/*.txt`
* AMR: `/data_origin/*_AMR.txt`, `/plain/data_mod/AMR_mod/*.txt`
* AFR: `/data_origin/*_AFR.txt`, `/plain/data_mod/AFR_mod/*.txt`

## Build Plain Models
We generate 1-hidden layer neural network models for several datasets represented by *population*. For each dataset, one can choose different models determined by *window_size*, which denotes the number of adjacent tag SNPs for each target SNP. Experiments based on the given dataset shows that the choice *window_size = 40* provides the best accuracy of genotype imputation.
Note that the larger *window_size* implies the higher computational cost. We also provide the multi-processing option for the acceleration so that one can dynamically choose the number of processes based on his/her computer environment.

To generate the imputation models, command the following:
```bash
$ cd ./plain
$ python New_gen_model_W.py -p <population> -w <window_size> -n <number_of_processes>
```
* `population`: Total, EUR, AMR, AFR 
* `window_size`: 8, 16, 24, 32, 40, 48, 56, 64, 72 
* `number_of_processes`: 1, 2, 4, 8, 16 

The generated plain models will be saved in the `/encrypted/<population>_DNNmodels/DNNmodels_<window_size>_c` directory.

### Model for Low Minor Allele Frequency data

We have added the model for data with low minor allele frequency (MAF). This data consists of 117k number of target SNPs. To run with low MAF data, run with additional option `-l 1` with `population`=Total. Then, the generated plain models will be saved in the `/encrypted/Total_DNNmodels/DNNmodels_lowMAF_<window_size>` directory.

```bash
$ cd ./plain
$ python New_gen_model_W.py -p Total -w <window_size> -n <number_of_processes> -l 1
```


## Build ModHEaaN
For the encryption of test data, we use the [ModHEaaN](https://github.com/idashSNU/Imputation/tree/master/ModHEaaN) library, which is a light-version implementation of the approximate HE scheme [CKKS](https://eprint.iacr.org/2016/421.pdf). Contrary to the original implementation of the [CKKS](https://eprint.iacr.org/2016/421.pdf) scheme [HEAAN](https://github.com/snucrypto/HEAAN), our [ModHEaaN](https://github.com/idashSNU/Imputation/tree/master/ModHEaaN) library does not have any dependency on multi-precision libraries GMP and NTL, and only supports homomorphic addition and 1-depth constant multiplication (hence bootstrapping disabled).

To build [ModHEaaN](https://github.com/idashSNU/Imputation/tree/master/ModHEaaN), command the following:
```bash
$ ./ModHEaaN/heaan
$ cmake CMakeLists.txt
$ make all
```

## Build the HE-based Genotype Imputation executable
```bash
$ ./encrypted/impute_dnn
$ cmake CMakeLists.txt
$ make all
```
The executable file is generated as `enc_impute` in the directory `./encrypted/impute_dnn`. 

## Run `enc_impute`
The command line to execute `enc_impute` depends on the target dataset. 
### For Total Dataset 
One need to choose two input arguments `window_size` and `number_of_targetSNP`. `window_size` is the number of adjacent tag SNPs for each target SNP, and `number_of_targetSNP` represents the number of target SNPs. 
```bash
$ cd ./encrypted/impute_dnn
$ ./enc_impute <window_size> <number_of_targetSNP>
```
* window_size: 8, 16, 24, 32, 40, 48, 56, 64, 72
* number_of_targetSNP: 20, 40, 80, 117

For instance, The argument below runs the genotype imputation of 80k target SNPs, with window size 40.
```bash
$ ./enc_impute 40 80
```

When you want to deal with low MAF data, When you want to deal with low MAF data, define number_of_targetSNP = 117 and set the other parameters by the same as original Total data.

### For EUR/AMR/AFR Dataset
In the case of EUR/AMR/AFR Dataset, the command line is fixed as 
```bash
$ ./enc_impute populations 80
```
The output will be genotype scores on EUR, AMR, and AFR datasets for the fixed `window_size=40` and `number_of_targetSNP=80`. 


* NOTE: The running time linearly grows up in terms of `window_size`, but the accuracy does not. For the given datasets, the choice `window_size=40` shows the best accuracy.


## Measure the Accuracy (MicroAUC)
If you succeed to run our solution, then the genotype score results of our solution "genotype_score" will be saved in `/encrypted/impute_dnn` denoted by `score_window<window_size>_<number_of_targetSNP>k.csv`. Note that the real genotypes of test data "genotype_real" is saved in `/plain` denoted by `real_<number_of_targetSNP>k.csv`. To run  `evaluation.py` in the `./plain` directory, command
```bash
python3 evaluation.py -i <genotype_score> -t <genotype_real> -o output.png
```

e.g. After running `$ ./enc_impute 40 80` in the previous step, then
```bash
$ cd ./plain
$ python3 evaluation.py -i ../encrypted/impute_dnn/score_window40_80k.csv -t real_80k.csv -o output.png
```

e.g. After running `$ ./enc_impute populations 80` in the previous step, then
```bash
$ cd ./plain
$ python3 evaluation.py -i ../encrypted/impute_dnn/score_EUR_80k.csv -t real_EUR_80k.csv -o output_EUR.png
$ python3 evaluation.py -i ../encrypted/impute_dnn/score_AMR_80k.csv -t real_AMR_80k.csv -o output_AMR.png
$ python3 evaluation.py -i ../encrypted/impute_dnn/score_AFR_80k.csv -t real_AFR_80k.csv -o output_AFR.png
```

e.g. (For low MAF data) After running `$ ./enc_impute 40 117` in the previous step, then
```bash
$ cd ./plain
$ python3 evaluation.py -i ../encrypted/impute_dnn/score_window40_117k.csv -t real_117k.csv -o output.png
```
