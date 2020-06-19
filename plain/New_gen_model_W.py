import csv
import numpy as np
import pandas as pd
import tensorflow as tf
import random
import sys, getopt
from keras.utils import np_utils
from keras import regularizers, optimizers
from sklearn.preprocessing import LabelEncoder
from sklearn.neighbors import KNeighborsClassifier
from keras import backend as K
from multiprocessing import Process, Queue

class locationList:
    def __init__(self, file_name, window_size, population):
        data_df = read_data2(file_name)
        self.population = population
        self.window_size = window_size

        self.num_list_x_start = []
        self.num_list_x_end = []
        # you can change the number of adjacent tag SNPs. This value affects accuracy and time.
        self.num_list_x_start, self.num_list_x_end = split_data(data_df, window_size)

        # print num_list_x_start
        
        self.num_list_y_start = []
        self.num_list_y_end = []

        # the number of target SNPs
        for j in range(0, 83072):
            self.num_list_y_start.append(j)
            self.num_list_y_end.append(j)

        print("file reading starts...")
        data_X_train = read_data("data_mod/"+ str(self.population)+"_mod/target_training_"+ str(self.population) +"_mod.txt")
        data_X_test = read_data("data_mod/"+ str(self.population)+"_mod/target_testing_"+ str(self.population) +"_mod.txt")
        self.data_Y_1 = [data_X_train, data_X_test]
        data_Y_train = read_data("data_mod/"+ str(self.population)+"_mod/tag_training_"+ str(self.population) +"_mod.txt")
        data_Y_test = read_data("data_mod/"+ str(self.population)+"_mod/tag_testing_"+ str(self.population) +"_mod.txt")
        self.data_X_1 = [data_Y_train, data_Y_test]
        print("file reading ends!")



# Read csv files.
def read_data(file_name):
    data = []
    with open(file_name, 'r') as txtfile:
        txtfile_reader = csv.reader(txtfile, delimiter='\t')
        for a in txtfile_reader:
            data.append(a)
    data = np.transpose(data)
    data_reshape = []
    for a in data:
        data_reshape.append(a.reshape(-1, a.shape[0]))
    data_train = data_reshape
    data_train = np.array(data_train)
    
    tmp = []
    for a in data_train:
        tmp.append(a[0])
    data_train = np.array(tmp)

    return data_train

def nogada_data(np_array, start1, end1):
    data_1 = np.zeros((len(np_array[0]), end1-start1))
    data_2 = np.zeros((len(np_array[1]), end1-start1))

    for i in range(0, len(np_array[0])):
        for j in range(0, end1-start1):
            data_1[i][j]=np_array[0][i][j+start1]
    for i in range(0, len(np_array[1])):
        for j in range(0, end1-start1):
            data_2[i][j]=np_array[1][i][j+start1]

    return [data_1, data_2]


# read the SNP location data.

def read_data2(file_name):
    data_df = pd.read_csv(file_name)
    return data_df


def max(a, b):
    if a > b:
        return a
    else:
        return b

# split tag SNPs for each target SNP

def split_data(data_df, window):

    # the name of category of location
    name_data_X = ['tag(total 16184)']   
    name_data_Y = ['target(total 83073)']
   
    X = np.array(data_df[name_data_X])
    Y = np.array(data_df[name_data_Y])

    # the number of tag SNPs
    num = 16184

    num_list_x_end = []
    num_list_x_start = []

    start_pt = 0;
    end_pt = window-1;

    index_speed = 0

    # This work find the adjacent tag SNPs for each target SNPs. Therefore, iteration number is the number of target SNPs.
    for i in range(0, 83072):
        if end_pt > num-3:
            num_list_x_start.append(start_pt)
            num_list_x_end.append(end_pt)
        else:
            min_distance = max(abs(Y[i][0]-X[start_pt][0]), abs(Y[i][0]-X[end_pt][0]))
#            print("check")
            
            while 1:
                temp = max(abs(Y[i][0]-X[start_pt+2][0]), abs(Y[i][0]-X[end_pt+2][0]))

    
                if temp > min_distance:
                    num_list_x_start.append(start_pt)
                    num_list_x_end.append(end_pt)
                    break
                else:
                    min_distance = temp
                    start_pt += 2
                    end_pt += 2
                if end_pt > num-3:
                    num_list_x_start.append(start_pt)
                    num_list_x_end.append(end_pt)
                    break
        
#        print(i, "-th iteration of sorting")
#        print(start_pt, end_pt)
  
    print("sorting ends!")
#    print(num_list_x_start)
    return num_list_x_start, num_list_x_end





def DNNmodel(numbering, ll):
    print(numbering,"-th iteration")
    real_y_num = ll.num_list_y_end[numbering]-ll.num_list_y_start[numbering]+1

    data_X = nogada_data(ll.data_X_1,ll.num_list_x_start[numbering] * 3, ll.num_list_x_end[numbering] * 3 + 3)
    data_Y = nogada_data(ll.data_Y_1,ll.num_list_y_start[numbering] * 3, ll.num_list_y_end[numbering] * 3 + 3)

    # you can change any setting of DNN models.

    model = tf.keras.models.Sequential()     
    model.add(tf.keras.layers.Dense(64, input_shape=(len(data_X[0][0]),), activation='linear', use_bias = False))    
    model.add(tf.keras.layers.Dropout(0.3))
    model.add(tf.keras.layers.Dense(3*(ll.num_list_y_end[numbering]-ll.num_list_y_start[numbering]+1), activation=tf.nn.sigmoid, use_bias = False))

    model.compile(optimizer='Nadam', loss='binary_crossentropy', metrics=['accuracy'])
    model.fit(data_X[0], data_Y[0], batch_size=32, epochs=40, validation_data=(data_X[1], data_Y[1]))  

        ####### weight matrices #######

    weight1 = model.layers[0].get_weights()
    weight2 = model.layers[2].get_weights()

    w1 = pd.DataFrame(weight1[0])           # first matrix
    w2 = pd.DataFrame(weight2[0])           # second matrix
    W = w1.dot(w2)
        # W = w1
    # The storage location of model W

    W.to_csv('../encrypted/'+str(ll.population)+'_DNNmodels/DNNmodels_' + str(ll.window_size) + '/W_New' + str(numbering)  + '.csv', header = False, index =False)
    
    del W
    del w1
    del w2
    del weight1
    del weight2
    del model
    del data_X
    del data_Y
    del real_y_num
    K.clear_session() 


# for numbering in range(1666, 82548):
# This can be changed according to your SNP data.

def DNNmodel_multiprocessing_16(ll):
    # multiprocessing number can be changed. our setting is 16.
    for i in range(0,5055):
        n1 = 1666+16*i
        n2 = n1+1
        n3 = n1+2
        n4 = n1+3
        n5 = n1+4
        n6 = n1+5
        n7 = n1+6
        n8 = n1+7
        n9 = n1+8
        n10 = n1+9
        n11 = n1+10
        n12 = n1+11
        n13 = n1+12
        n14 = n1+13
        n15 = n1+14
        n16 = n1+15       
        procs = []
        procs.append(Process(target=DNNmodel, args=(n1,ll)))
        procs.append(Process(target=DNNmodel, args=(n2,ll)))
        procs.append(Process(target=DNNmodel, args=(n3,ll)))
        procs.append(Process(target=DNNmodel, args=(n4,ll)))
        procs.append(Process(target=DNNmodel, args=(n5,ll)))
        procs.append(Process(target=DNNmodel, args=(n6,ll)))
        procs.append(Process(target=DNNmodel, args=(n7,ll)))
        procs.append(Process(target=DNNmodel, args=(n8,ll)))
        procs.append(Process(target=DNNmodel, args=(n9,ll)))
        procs.append(Process(target=DNNmodel, args=(n10,ll)))
        procs.append(Process(target=DNNmodel, args=(n11,ll)))
        procs.append(Process(target=DNNmodel, args=(n12,ll)))
        procs.append(Process(target=DNNmodel, args=(n13,ll)))
        procs.append(Process(target=DNNmodel, args=(n14,ll)))
        procs.append(Process(target=DNNmodel, args=(n15,ll)))
        procs.append(Process(target=DNNmodel, args=(n16,ll)))        
        for p in procs:
            p.start()
        for p in procs:
            p.join()

    procs = []
    n1=82546
    n2=82547
    n3=82548
    procs.append(Process(target=DNNmodel, args=(n1,ll)))
    procs.append(Process(target=DNNmodel, args=(n2,ll)))
    procs.append(Process(target=DNNmodel, args=(n3,ll)))
    for p in procs:
        p.start()
    for p in procs:
        p.join()


def DNNmodel_multiprocessing_8(ll):
    # multiprocessing number can be changed. our setting is 16.
    for i in range(0,10110):
        n1 = 1666+8*i
        n2 = n1+1
        n3 = n1+2
        n4 = n1+3
        n5 = n1+4
        n6 = n1+5
        n7 = n1+6
        n8 = n1+7
        procs = []
        procs.append(Process(target=DNNmodel, args=(n1,ll)))
        procs.append(Process(target=DNNmodel, args=(n2,ll)))
        procs.append(Process(target=DNNmodel, args=(n3,ll)))
        procs.append(Process(target=DNNmodel, args=(n4,ll)))
        procs.append(Process(target=DNNmodel, args=(n5,ll)))
        procs.append(Process(target=DNNmodel, args=(n6,ll)))
        procs.append(Process(target=DNNmodel, args=(n7,ll)))
        procs.append(Process(target=DNNmodel, args=(n8,ll)))
        for p in procs:
            p.start()
        for p in procs:
            p.join()

    procs = []
    n1=82546
    n2=82547
    n3=82548
    procs.append(Process(target=DNNmodel, args=(n1,ll)))
    procs.append(Process(target=DNNmodel, args=(n2,ll)))
    procs.append(Process(target=DNNmodel, args=(n3,ll)))
    for p in procs:
        p.start()
    for p in procs:
        p.join()

def DNNmodel_multiprocessing_4(ll):
    # multiprocessing number can be changed. our setting is 16.
    for i in range(0,20220):
        n1 = 1666+4*i
        n2 = n1+1
        n3 = n1+2
        n4 = n1+3
        procs = []
        procs.append(Process(target=DNNmodel, args=(n1,ll)))
        procs.append(Process(target=DNNmodel, args=(n2,ll)))
        procs.append(Process(target=DNNmodel, args=(n3,ll)))
        procs.append(Process(target=DNNmodel, args=(n4,ll)))
        for p in procs:
            p.start()
        for p in procs:
            p.join()

    procs = []
    n1=82546
    n2=82547
    n3=82548
    procs.append(Process(target=DNNmodel, args=(n1,ll)))
    procs.append(Process(target=DNNmodel, args=(n2,ll)))
    procs.append(Process(target=DNNmodel, args=(n3,ll)))
    for p in procs:
        p.start()
    for p in procs:
        p.join()

def DNNmodel_multiprocessing_2(ll):
    # multiprocessing number can be changed. our setting is 16.
    for i in range(0,40441):
        n1 = 1666+2*i
        n2 = n1+1
        procs = []
        procs.append(Process(target=DNNmodel, args=(n1,ll)))
        procs.append(Process(target=DNNmodel, args=(n2,ll)))
        for p in procs:
            p.start()
        for p in procs:
            p.join()

    procs = []
    n1=82548
    procs.append(Process(target=DNNmodel, args=(n1,ll)))
    for p in procs:
        p.start()
    for p in procs:
        p.join()

def DNNmodel_multiprocessing_1(ll):
    # multiprocessing number can be changed. our setting is 16.
    for i in range(1666,82549):
        n1 = i
        procs = []
        procs.append(Process(target=DNNmodel, args=(n1,ll)))
        for p in procs:
            p.start()
        for p in procs:
            p.join()


def main(argv):
    window_size = 0
    nthread = 0
    try:
        opts, args = getopt.getopt(argv,"hp:w:n:",["population=","window_size=","nthread="])
    except getopt.GetoptError:
        print('python evaluation.py -p <population> -w <window_size> -n <nthread>')
        sys.exit(2)
    for opt, arg in opts:
        if opt == '-h':
            print('python New_gen_model_W_total.py -p population -w window_size -n nthread')
            sys.exit()
        elif opt in ("-p", "--population"):
            population = arg
        elif opt in ("-w", "--window_size"):
            window_size = int(arg)
        elif opt in ("-n", "--nthread"):
            nthread = int(arg)
        
    print('Population is '+ str(population))
    print('Window_size is '+ str(window_size))
    print('The number of threads is '+ str(nthread))

    # Generate true label list for our personal experiment (Not Required for Real Evaluation)

    # data_Y_train = read_data("data_mod/target_training_mod.txt")

    # data_Y_test = read_data("data_mod/target_testing_mod.txt")

    # data_Y_2 = [data_Y_train, data_Y_test]

    file_name = "data_mod/idash_local.csv"
    ll = locationList(file_name, window_size, population)


    

    if nthread == 16 :
        DNNmodel_multiprocessing_16(ll)
    elif nthread == 8:
        DNNmodel_multiprocessing_8(ll)
    elif nthread == 4:
        DNNmodel_multiprocessing_4(ll)
    elif nthread == 2:
        DNNmodel_multiprocessing_2(ll)
    elif nthread == 1:
        DNNmodel_multiprocessing_1(ll)
    else:
        print("This number does not support multiprocessing.")
    
    

if __name__ == "__main__":
    main(sys.argv[1:])
    
    





