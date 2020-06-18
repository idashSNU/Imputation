SNU Team IDASH2019

This is a plain code to build 1-hidden layer neural network models from the SNU team for HE-based genotype imputation.
It consists of two directories:  /plain, /data_mod.
In /plain directory, there are some model generation python files whose names start with 'New_gen_model_W_', and these models are saved in /encrypted directory and /data_mod includes all SNP data.

In our model generation python files provide multiprocessing module which can change the number of multiprocess in the function, DNNmodel_multiprocessing.

Our plain code consists of two part.

First, our code divides the tag SNPs adjacent to the target SNP using the function split_data.
split_data exploits the number of tag SNPs for each target SNP which is chosen by users in line 132.
You can make 1-hidden layer neural network models by changing the "window_size" according to the data. 

Second, our code generates1-hidden layer neural network models to predict target SNPs.
You can change the condition of DNN models such as the number of nodes, the rate of dropout, activation function, epochs, etc in line 163~169.
In function "DNNmodel", you change the location of 1-hidden layer neural network model W for each target SNPs in line 181.
 
Note that when using SNP data, 0,1,2 must be converted to 100, 010, 001 and used.
To split data, you make a csv files which contains the location of tag SNPs and target SNPs.



The following is the way to build and run our solution:


1. Choose the window size(8,16,24,32,40,48,56,64,72) in line 15. The default setting is 40.

2.Command $python New_gen_model_W_*.py

