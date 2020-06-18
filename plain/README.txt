SNU Team IDASH2019

Our solution is based on (1-hidden layer) neural network models.
- The python files 'New_gen_model_W_*.py' generates model csv files, which will be saved in /encrypted directory. 
- The python file 'evaluation.py' is for evaluating (encrypted) imputation results. See also README.txt in the parent folder.
- The folder /data_train contains SNP data for model training. 

Model generation files ('New_gen_model_W_*.py') mainly perform two tasks.

First, our code divides the tag SNPs adjacent to the target SNP using the function split_data.
split_data exploits the number of tag SNPs for each target SNP which is chosen by users in line 132.
You can make 1-hidden layer neural network models by changing the "window_size" according to the data. 

Second, our code generates1-hidden layer neural network models to predict target SNPs.
You can change the condition of DNN models such as the number of nodes, the rate of dropout, activation function, epochs, etc in line 163~169.
In function "DNNmodel", you change the location of 1-hidden layer neural network model W for each target SNPs in line 181.

* We also support multiprocessing by function "DNNmodel_multiprocessing", where you can change the number of multiprocess by input.
 
* To run our code with your own data, your data should be one-hot-encoded. In other words, you have to convert the SNP data consisting of 0, 1, 2 into 100, 010, 001.

* To split data, you make a csv files which contains the location of tag SNPs and target SNPs.

The following is the way to build and run our solution:


1. Choose the window size(8,16,24,32,40,48,56,64,72) in line 15. The default setting is 40.

2.Command $python New_gen_model_W_*.py

