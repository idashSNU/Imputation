SNU Team IDASH2019

This is a solution from the SNU team for HE-based genotype imputation. It consists of four directories: /ModHEaaN, /data_origin, /plain and /encrypted.
/ModHEaaN includes a library which implements the CKKS homomorphic encryption scheme (without RNS), and /data_origin includes all snp data. In /plain directory, there are some python files that generates model files in name 'New_gen_model_W_*.py', and these models are saved in /encrypted directory. You can find our main code for HE-based genotype imputation in /encrypted. 


The following is the way to build and run our solution:


1. Build HE-imputation executable by
$ cd
$ cd encrypted/impute_dnn
$ cmake CMakeLists.txt
$ make all

2. Run the executable with command line arguments by
$ ./enc_impute `window_size' `# target SNP in k'.

The parameter 'window_size' is the number of nearby (known) SNPs used for each target SNP, and '# target SNP in k' is literally the number of target SNPs. You can choose `window_size' among 8, 16, 24 , ... , 72, and `# target SNP in k' among 20 / 40 / 80.

e.g. `$ ./enc_impute 40 80' runs imputation of 80k target SNPS, with window size 40.

* Note : The running time of (encrypted) imputation grows with 'window_size', but the quality of imputation (accuracy) does not. For the sample data, a moderate window size (about 40) shows the best result.

* WARN: `populations' mode is only compatible with '80' (i.e., $ ./enc_impute populations 80)

3. Measure MicroAUC
If you succeed to run our enc_impute solution (Step 2), then the score output file will be saved in /encrypted/impute_dnn in name 'score_window(window_size)__(# target snp in k)k.csv". Note that the real genotypes of test data is saved in /plain named as 'real_(# target snp in k)k.csv'.

e.g. if you ran `$ ./enc_impute 40 80' in Step 2, then
$ cd
$ cd plain
$ python3 evaluation.py -i ../encrypted/impute_dnn/score_window40_80k.csv -t real_80k.csv -o output.png

e.g. if you ran `$ ./enc_impute populations 80' in Step 2, then
$ cd
$ cd plain
$ python3 evaluation.py -i ../encrypted/impute_dnn/score_EUR_80k.csv -t real_EUR_80k.csv -o output_EUR.png
$ python3 evaluation.py -i ../encrypted/impute_dnn/score_AMR_80k.csv -t real_AMR_80k.csv -o output_AMR.png
$ python3 evaluation.py -i ../encrypted/impute_dnn/score_AFR_80k.csv -t real_AFR_80k.csv -o output_AFR.png

