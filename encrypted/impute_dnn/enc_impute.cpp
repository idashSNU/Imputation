#include "../src/Evaluator.h"
#include "../src/FileReader.h"
#include "../src/PrintUtils.h"
#include "../src/TimeUtils.h"
#include "../src/MemoryUsage.h"

void runEvaluation(Evaluator& evaluator, string tag_filename, SNP_Parameters& snp_params, long num_sample__, long datatype, long data_version, bool olddata, long population, Message*** result_matrix)
{
    TimeUtils timeutils_tot;
    TimeUtils timeutils;
    cout << "[runEvaluation] start num_sample__ = " << num_sample__ << endl;
    double** mat_X = new double*[num_sample__]();
    for (int i = 0; i < num_sample__; i++) {
        mat_X[i] = new double[3 * snp_params.num_tag_snp]();
    }
    cout << "[runEvaluation] Read and Encode X from" << tag_filename << endl;
    FileReader::readAndEncodeX(tag_filename, mat_X, snp_params.num_sample, num_sample__, snp_params.num_tag_snp, olddata);
    
    cout << "[runEvaluation] Read W with datatype " << datatype << ", data_version " << data_version << ", snp_params.num_sample = " << snp_params.num_sample << endl;
    double*** W = new double**[snp_params.num_target_snp];
    #pragma omp parallel for
    for (int i = 0; i < snp_params.num_target_snp; i++) {
        W[i] = new double*[snp_params.list_nSNP_tag[i]];
        for (int k = 0; k < snp_params.list_nSNP_tag[i]; k++) {
            W[i][k] = new double[snp_params.list_nSNP_target[i]];
        }
    }
    evaluator.readW_nosave(datatype, data_version, snp_params, population, W);
    
    timeutils.start("EncryptX_imag_nosave");
    Ciphertext** EncX = new Ciphertext*[(3 * snp_params.num_tag_snp + 1) / 2];
    evaluator.encryptX_imag_nosave_newdata(mat_X, snp_params, EncX);
    timeutils.stop("EncryptX_imag_nosave");

    for (int i = 0; i < num_sample__; i++) {
        delete[] mat_X[i];
    }
    delete[] mat_X;

    timeutils.start("compute_XW_imag_nosave");
    Ciphertext*** EncXW = new Ciphertext**[snp_params.num_target_snp];
    evaluator.compute_XW_imag_nosave_newdata(EncX, W, snp_params, EncXW);
    timeutils.stop("compute_XW_imag_nosave");
    
    for (int i = 0; i < snp_params.num_target_snp; i++) {
        for (int j = 0; j < snp_params.list_nSNP_tag[i]; j++) {
            delete[] W[i][j];
        }
        delete[] W[i];
    }
    delete[] W;
    
    for (int i = 0; i < (3 * snp_params.num_tag_snp + 1) / 2; i++) {
        delete[] EncX[i];
    }
    delete[] EncX;

    evaluator.decryptXW_nosave_newdata(result_matrix, EncXW, snp_params);
}

int main(int argc, char* argv[])
{
    //default setting
    long logN = 10;
    long logQ = 27;
    long datatype = 0;

    // long logN = atoi(argv[1]);
    // long logQ = atoi(argv[2]);
    // long datatype = atoi(argv[3]);

    long data_version = 1;          // default setting
    // data_version = atoi(argv[4]);

    bool olddata = 0;               // default setting
    // olddata = atoi(argv[5]);
    bool population = 0;            // default setting
    if(!string(argv[1]).compare("populations")) population = 1;

    long num_target_snp_k = atoi(argv[2]);

    long logp = 17;
    long max_num_sample = 1024;
    
    long window;
    string mode; 

    if(!population)
	window = stoi(argv[1]);
    else
	mode = argv[1];

    // if (mode == "fast") window = 8;
    // else if (mode == "balanced") window = 24;
    // else {
   //  	window = 40;
    // 	if (mode == "populations") population = 1;
//	}

    TimeUtils timeutils;
    TimeUtils timeutils_tot;

    cout << "Encrypted Imputation using DNN model" << endl;

    timeutils.start("Scheme Generation");
    Evaluator evaluator(logN, logp, logQ);
    timeutils.stop("Scheme Generation");    


    // SNP_Parameters snp_params(datatype, data_version, window);
    SNP_Parameters snp_params(datatype, data_version, window, num_target_snp_k);
    FileReader::read_x_start(snp_params);
    snp_params.setLists_newdata();

    string encode_foldername = "../Models_Encoded";
    string data_foldername = "../Ctxts";
    string tag_filename;
    string tag_filename_AFR, tag_filename_AMR, tag_filename_EUR;
    long num_sample__, num_sample_AFR__, num_sample_AMR__, num_sample_EUR__;

    if (olddata) {
        tag_filename = "../../data_origin/sorted_tag_SNPs_" + to_string(datatype) + "k_genotypes.txt";
        snp_params.num_sample = 400;
        num_sample__ = snp_params.num_sample;
    } else if (population) {
        tag_filename_AFR = "../../data_origin/population/tag_testing_AFR.txt";
        tag_filename_AMR = "../../data_origin/population/tag_testing_AMR.txt";
        tag_filename_EUR = "../../data_origin/population/tag_testing_EUR.txt";
        num_sample_AFR__ = 272;
        num_sample_AMR__ = 135;
        num_sample_EUR__ = 210;
        // snp_params.num_target_snp = 80882;

    } else {
        if (num_target_snp_k <= 80) {
            tag_filename = "../../data_origin/tag_testing_80k.txt"; // Evaluator guys have to change here
        } else {
            tag_filename = "../../data_origin/tag_testing_117k.txt"; // Evaluator guys have to change here
        }
        cout << "check filename = " << tag_filename << endl;
        snp_params.num_sample = 1004;
        num_sample__ = snp_params.num_sample;
    }

    long Ylength;
    if (olddata) {
        Ylength = 500;
    } else { // population or newData
        if (num_target_snp_k <= 80) {
            Ylength = 83072;
        } else {
	        Ylength = 117904;
        }
    }

    Message*** dec_XW = new Message**[snp_params.num_target_snp];
    Message*** dec_XW_AFR = new Message**[snp_params.num_target_snp];
    Message*** dec_XW_AMR = new Message**[snp_params.num_target_snp];
    Message*** dec_XW_EUR = new Message**[snp_params.num_target_snp];
    for (int i = 0; i < snp_params.num_target_snp; i++) {
    	dec_XW[i] = new Message*[snp_params.list_nSNP_target[i]];
    	dec_XW_AFR[i] = new Message*[snp_params.list_nSNP_target[i]];
    	dec_XW_AMR[i] = new Message*[snp_params.list_nSNP_target[i]];
    	dec_XW_EUR[i] = new Message*[snp_params.list_nSNP_target[i]];
    }


    size_t currentmemory = getCurrentRSS() >> 20;
    size_t peakmemory = getPeakRSS() >> 20;
    cout << "------------------" << endl;
    cout << "Current Memory Usage before Genotype Imputation: " << currentmemory << "MB" << endl;
    cout << "Peak Memory Usage before Genotype Imputation: " << peakmemory << "MB" << endl;
    cout << "------------------" << endl;


    if (population) {
        snp_params.num_sample = num_sample_AFR__;
        runEvaluation(evaluator, tag_filename_AFR, snp_params, num_sample_AFR__, datatype, data_version, olddata, 1, dec_XW_AFR);
        snp_params.num_sample = num_sample_AMR__;
        runEvaluation(evaluator, tag_filename_AMR, snp_params, num_sample_AMR__, datatype, data_version, olddata, 2, dec_XW_AMR);
        snp_params.num_sample = num_sample_EUR__;
        runEvaluation(evaluator, tag_filename_EUR, snp_params, num_sample_EUR__, datatype, data_version, olddata, 3, dec_XW_EUR);
    } else {
        runEvaluation(evaluator, tag_filename, snp_params, num_sample__, datatype, data_version, olddata, 0, dec_XW);

    }

    cout << "Genotype Score Imputation Done!" << endl;



    currentmemory = getCurrentRSS() >> 20;
    peakmemory = getPeakRSS() >> 20;
    cout << "------------------" << endl;
    cout << "Current Memory Usage after Genotype Imputation: " << currentmemory << "MB" << endl;
    cout << "Peak Memory Usage after Genotype Imputation: " << peakmemory << "MB" << endl;
    cout << "------------------" << endl;

    cout << "Save genotype score matrix as csv files in order to get result via micro AUC" << endl;
    string target_filename;
    string target_filename_AFR, target_filename_AMR, target_filename_EUR;
    
    string** mat_Y = new string*[Ylength];
    string** mat_Y_AFR = new string*[Ylength];
    string** mat_Y_AMR = new string*[Ylength];
    string** mat_Y_EUR = new string*[Ylength];
    for (int i = 0; i < Ylength; i++) {
        mat_Y[i] = new string[2510];
        mat_Y_AFR[i] = new string[num_sample_AFR__ + 4];
        mat_Y_AMR[i] = new string[num_sample_AMR__ + 4];
        mat_Y_EUR[i] = new string[num_sample_EUR__ + 4];
    }
    if (olddata) {
        target_filename = "../../data_origin/sorted_target_SNP_genotypes.txt";
        cout << "read target in " << target_filename << endl;
        FileReader::readMatrix(target_filename, mat_Y, 0, '\t', 0);
    }
    else if (population) {
        target_filename_AFR = "../../data_origin/population/target_testing_AFR.txt";
        target_filename_AMR = "../../data_origin/population/target_testing_AMR.txt";
        target_filename_EUR = "../../data_origin/population/target_testing_EUR.txt";
        cout << "read target in " << target_filename_AFR << endl;
        FileReader::readMatrix(target_filename_AFR, mat_Y_AFR, 0, '\t', 0);
        cout << "read target in " << target_filename_AMR << endl;
        FileReader::readMatrix(target_filename_AMR, mat_Y_AMR, 0, '\t', 0);
        cout << "read target in " << target_filename_EUR << endl;
        FileReader::readMatrix(target_filename_EUR, mat_Y_EUR, 0, '\t', 0);
    } else {
        if (num_target_snp_k <= 80) {
            target_filename = "../../data_origin/target_testing_80k.txt";
        } else {
            target_filename = "../../data_origin/target_testing_117k.txt";
        }
        cout << "read target in " << target_filename << endl;
        FileReader::readMatrix(target_filename, mat_Y, 0, '\t', 0);
    }


    cout << "set snp_data" << endl;
    string* snp_data = new string[snp_params.num_target_snp];
    if (olddata) {
    	for (int i = 0; i < snp_params.num_target_snp; i++) {
        	snp_data[i] = mat_Y[i][3];
        }
    } else if (population) {
        for (int i = 0; i < snp_params.num_target_snp; i++) {
        	snp_data[i] = mat_Y_AMR[i + 1666][3];
        }
    } else {
    	for (int i = 0; i < snp_params.num_target_snp; i++) {
            if (num_target_snp_k <= 80) {
        	    snp_data[i] = mat_Y[i + 1666][3];
            }
            else {
                snp_data[i] = mat_Y[i][3];
            }
        }
    }
    
    if (population) {
        // evaluator.computeScore(dec_XWold, snp_data, data_version, snp_params);
        evaluator.computeScore_population(dec_XW_AFR, dec_XW_AMR, dec_XW_EUR, snp_data, snp_params);
    } else {
        evaluator.computeScore_newdata(dec_XW, snp_data, data_version, snp_params);
    }

    // currentmemory = getCurrentRSS() >> 20;
    // peakmemory = getPeakRSS() >> 20;
    // cout << "------------------" << endl;
    // cout << "Current Memory Usage after saving Genotype Score in csv files: " << currentmemory << "MB" << endl;
    // cout << "Peak Memory Usage after saving Genotype Score in csv files: " << peakmemory << "MB" << endl;
    // cout << "------------------" << endl;
  	
    cout << "Done.." << endl;

    return 0;
}
