SNU Team IDASH2019

This is a solution from the SNU team for HE-based genotype imputation. It consists of four directories: /HEAAN_NTT, /data_origin, /plain and /encrypted.
/HEAAN_NTT includes a library which implements the CKKS homomorphic encryption scheme (without RNS), and /data_origin includes all snp data. In /plain directory, there are some model generation python files whose names start with 'New_gen_model_W_', and these models are saved in /encrypted directory. You can find our main code for HE-based genotype imputation in /encrypted. 


The following is the way to build and run our solution:


1. Build HE-imputation executable by
$ cd
$ cd encrypted/impute_dnn
$ cmake CMakeLists.txt
$ make all

2. Run the executable with command line arguments by
$ ./enc_impute `optimize_mode' `# target SNP in k'.

You can choose `optimize_mode' among fast / balanced / best / populations,
and `# target SNP in k' among 20 / 40 / 80.

e.g. `$ ./enc_impute best 80' runs the best microAUC version (but slowest), with 80k target SNPs

* WARN: `populations' mode is only compatible with '80' (i.e., $ ./enc_impute populations 80)


3. Measure MicroAUC
When you success to run our enc_impute solution (Step 2), then the score output file will be saved in /encrypted/impute_dnn as 'score_(model)__(# target snp in k)k.csv". Note that the real genotypes of test data is saved in /plain named as 'real_(# target snp in k)k.csv'.

e.g. if you ran `$ ./enc_impute best 80' in Step 2, then
$ cd
$ cd plain
$ python3 evaluation.py -i ../encrypted/impute_dnn/score_best_80k.csv -t real_80k.csv -o output.png

e.g. if you ran `$ ./enc_impute populations 80' in Step 2, then
$ cd
$ cd plain
$ python3 evaluation.py -i ../encrypted/impute_dnn/score_EUR_80k.csv -t real_EUR_80k.csv -o output_EUR.png
$ python3 evaluation.py -i ../encrypted/impute_dnn/score_AMR_80k.csv -t real_AMR_80k.csv -o output_AMR.png
$ python3 evaluation.py -i ../encrypted/impute_dnn/score_AFR_80k.csv -t real_AFR_80k.csv -o output_AFR.png

