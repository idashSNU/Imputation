SNU Team IDASH2019

(yong: 전체적으로, command line argument를 받아오도록 변경이 되면 좋을 것 같은데, 이러면 README가 확 바뀌어야될거같아서 여기 README는 많이 건드리질 못하겠음.
$python New_gen_model_W_*.py 'window_size' '# threads for multiprocess' 'path for model files' 같은 느낌으로?
찾아보니 이런 방법이 있는듯:
https://blog.naver.com/PostView.nhn?blogId=wideeyed&logNo=221400066328&from=search&redirect=Log&widgetTypeCall=true&directAccess=false)

Our solution is based on (1-hidden layer) neural network models.
- The python files 'New_gen_model_W_*.py' generates model csv files, which will be saved in /encrypted directory. 
- The python file 'evaluation.py' is for evaluating (encrypted) imputation results. See also README.txt in the parent folder.
- The folder /data_train contains SNP data for model training. 

The following is the way to generate model csv files:

1. Choose the window size(8,16,24,32,40,48,56,64,72) in line 15. The default setting is 40.

2. Command $python New_gen_model_W_*.py

Model generation files ('New_gen_model_W_*.py') mainly perform two tasks.

First, our code divides the tag SNPs adjacent to the target SNP using the function split_data.
split_data exploits the number of tag SNPs for each target SNP which is chosen by users in line 132.
You can make 1-hidden layer neural network models by changing the "window_size" according to the data. 
(yong: 가능하면 split과 training을 분리해서 command-line argument에서 선택할수 있게 제공하면 좋을 것 같은데.. )

* WARN: To split data, you make a csv files which contains the location of tag SNPs and target SNPs.
(yong: 이거 영어가 이해가 안되는데, location파일이 '있어야한다' 인건가?? ㅠ)

Second, our code generates 1-hidden layer neural network models to predict target SNPs.
You can change hyperparameters of network such as the number of nodes, the rate of dropout, activation function, epochs, etc in line 163~169.
In function "DNNmodel", you can change the location of 1-hidden layer neural network model W for each target SNPs in line 181.

* We support multiprocessing on training, where you can change the number of threads in codes.

* WARN: To run our code with your own data, it should be one-hot-encoded. In other words, you have to convert the SNP data consisting of 0, 1, 2 into 100, 010, 001.
