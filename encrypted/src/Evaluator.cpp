#include "Evaluator.h"
#include "SNP_Parameters.h"
#include "TimeUtils.h"


chrono::high_resolution_clock::time_point t1, t2;

#define START() t1 = chrono::high_resolution_clock::now();
#define END() t2 = chrono::high_resolution_clock::now();
#define PRINTTIME(msg, repeat) cout << "* " << msg << " time = " << (double)chrono::duration_cast<chrono::milliseconds>(t2 - t1).count() / repeat << " ms" << endl;

Evaluator::Evaluator(long _logN, long _logp, long _logQ)
    : logN(_logN)
    
    , logp(_logp)
    , logQ(_logQ)
    , logSlots(_logN - 1)
{
    // num_packing = (1 << logSlots) / 512;
    omp_set_num_threads(16);
    params = Parameters(1ULL << logN, logQ, logp, 1);
    context = Context(params);
    secret_key = SecretKey(context);
    string key_dir_path = "./Key";
    mkdir(key_dir_path.data(), 0775);
    public_key_pack = PublicKeyPack(context, secret_key, key_dir_path);
    pk = public_key_pack.getEncKey();

    encryptor = Encryptor(context);
    decryptor = Decryptor(context);
    homEvaluator = HomEvaluator(context);
}


void Evaluator::computeScore(Message** dec_XW, string* snp_data, long data_version, SNP_Parameters snp_params)
{
    double* score = new double[snp_params.num_sample * snp_params.num_target_snp * 3];
    long* block_size = new long[snp_params.num_target_snp];
    for (int i = 0; i < snp_params.num_target_snp; i++) {
        for (int k = 0; k < snp_params.list_nSNP_target[i]; k++) {
            for (int j = 0; j < snp_params.num_sample; j++) {
                score[3 * snp_params.num_target_snp * j + k + 3 * snp_params.list_y_start[i]] = dec_XW[i][k][j].real();
            }
        }
    }

    ofstream fout;
    long length = snp_params.num_sample * snp_params.num_target_snp;
    string output_filename = "score_" + to_string(logN) + "_" + to_string(snp_params.datatype) + "k_old";
    if (data_version == 2) {
        output_filename += "_2";
    }
    output_filename += ".csv";
    fout.open(output_filename);
    fout << "Subject ID,Target SNP,0,1,2\n";
    for (int i = 0; i < snp_params.num_sample * snp_params.num_target_snp; ++i) {
        fout << i / snp_params.num_target_snp + 1 << ",";
        fout << snp_data[i % snp_params.num_target_snp] << ",";
        fout << score[3 * i] << "," << score[3 * i + 1] << "," << score[3 * i + 2] << '\n';
    }
    fout.close();

    cout << "Final scores are saved in " << output_filename << endl;
}


void Evaluator::computeScore_newdata(Message*** dec_XW, string* snp_data, long data_version, SNP_Parameters snp_params)
{
    long nslots = 1 << logSlots;

    // number of samples
    
    long num_packing = (snp_params.num_sample - 1) / nslots + 1;
    int remain = snp_params.num_sample % nslots;
    
    // cout << "[Evaluator::computeScore_newdata] num_sample = " << snp_params.num_sample << endl;
    // cout << "[Evaluator::computeScore_newdata] num_target_snp = " << snp_params.num_target_snp << endl;cout << "[Evaluator::computeScore_newdata] num_packing = " << num_packing << endl;
    // cout << "[Evaluator::computeScore_newdata] num_remain = " << remain << endl;
    
    double* score = new double[snp_params.num_sample * snp_params.num_target_snp * 3];
    long* block_size = new long[snp_params.num_target_snp];
    
    for (int i = 0; i < snp_params.num_target_snp; i++) {
        for (int k = 0; k < snp_params.list_nSNP_target[i]; k++) {
            for (int l = 0; l < num_packing - 1 ; l++) {
                for (int j = 0; j < nslots; j++) {
                    score[3 * snp_params.num_target_snp * (j + l * nslots) + k + 3 * i] = dec_XW[i][k][l][j].real();
                }
            }
            for (int j = 0; j < remain; j++) {
                score[3 * snp_params.num_target_snp * (snp_params.num_sample - remain + j) + k + 3 * i] = dec_XW[i][k][num_packing-1][j].real();
            }
        }
    }

    ofstream fout;
    long length = snp_params.num_sample * snp_params.num_target_snp;
    string output_filename = "score";

    output_filename += "_window" + to_string(snp_params.window);
  //  if (snp_params.window == 8) {
  //  	output_filename += "_fast";
  //  } else if (snp_params.window == 24) {
  //  	output_filename += "_balanced";
  //  } else if (snp_params.window == 40) {
  //  	output_filename += "_best";
  //  }
    if (snp_params.num_target_snp < 30000) {
        output_filename += "_20k";
    } else if (snp_params.num_target_snp < 50000) {
        output_filename += "_40k";
    } else if (snp_params.num_target_snp > 70000 && snp_params.num_target_snp < 90000) {
        output_filename += "_80k";
    }
    output_filename += ".csv";
    fout.open(output_filename);
    fout << "Subject ID,Target SNP,0,1,2\n";
    for (int i = 0; i < snp_params.num_sample * snp_params.num_target_snp; ++i) {
        fout << i / snp_params.num_target_snp + 1 << ",";
        fout << snp_data[i % snp_params.num_target_snp] << ",";    
        fout << score[3 * i] << "," << score[3 * i + 1] << "," << score[3 * i + 2] << '\n';
    }
    fout.close();

    // cout << "[Evaluator::computeScore_newdata] num_sample = " << snp_params.num_sample << endl;
    // cout << "[Evaluator::computeScore_newdata] num_target_snp = " << snp_params.num_target_snp << endl;
    // cout << "[Evaluator::computeScore_newdata] total length of score = " << snp_params.num_sample * snp_params.num_target_snp << endl;
    cout << "[Evaluator::computeScore_newdata] Final scores are saved in " << output_filename << endl;
}

void Evaluator::computeScore_population(Message*** dec_XW_AFR, Message*** dec_XW_AMR, Message*** dec_XW_EUR, string* snp_data, SNP_Parameters snp_params)
{
    long nslots = 1 << logSlots;

    // number of samples
    long num_sample_AFR = 272;
    long num_sample_AMR = 135;
    long num_sample_EUR = 210;
    
    long num_packing_AFR = (num_sample_AFR - 1) / nslots + 1;
    long num_packing_AMR = (num_sample_AMR - 1) / nslots + 1;
    long num_packing_EUR = (num_sample_EUR - 1) / nslots + 1;
    int remain_AFR = num_sample_AFR % nslots;
    int remain_AMR = num_sample_AMR % nslots;
    int remain_EUR = num_sample_EUR % nslots;

    double* score_AFR = new double[num_sample_AFR * snp_params.num_target_snp * 3];
    double* score_AMR = new double[num_sample_AMR * snp_params.num_target_snp * 3];
    double* score_EUR = new double[num_sample_EUR * snp_params.num_target_snp * 3];
    long* block_size = new long[snp_params.num_target_snp];

    for (int i = 0; i < snp_params.num_target_snp; i++) {
        for (int k = 0; k < snp_params.list_nSNP_target[i]; k++) {
            // **********
            // score_AFR 
            for (int l = 0; l < num_packing_AFR - 1 ; l++) {
                for (int j = 0; j < nslots; j++) {
                    score_AFR[3 * snp_params.num_target_snp * (j + l * nslots) + k + 3 * i] = dec_XW_AFR[i][k][l][j].real();
                }
            }
            for (int j = 0; j < remain_AFR; j++) {
                score_AFR[3 * snp_params.num_target_snp * (num_sample_AFR - remain_AFR + j) + k + 3 * i] = dec_XW_AFR[i][k][num_packing_AFR-1][j].real();
            }

            // **********
            // score_AMR 
            for (int l = 0; l < num_packing_AMR - 1 ; l++) {
                for (int j = 0; j < nslots; j++) {
                    score_AMR[3 * snp_params.num_target_snp * (j + l * nslots) + k + 3 * i] = dec_XW_AMR[i][k][l][j].real();
                }
            }
            for (int j = 0; j < remain_AMR; j++) {
                score_AMR[3 * snp_params.num_target_snp * (num_sample_AMR - remain_AMR + j) + k + 3 * i] = dec_XW_AMR[i][k][num_packing_AMR-1][j].real();
            }

            // **********
            // score_EUR
            for (int l = 0; l < num_packing_EUR - 1 ; l++) {
                for (int j = 0; j < nslots; j++) {
                    score_EUR[3 * snp_params.num_target_snp * (j + l * nslots) + k + 3 * i] = dec_XW_EUR[i][k][l][j].real();
                }
            }
            for (int j = 0; j < remain_EUR; j++) {
                score_EUR[3 * snp_params.num_target_snp * (num_sample_EUR - remain_EUR + j) + k + 3 * i] = dec_XW_EUR[i][k][num_packing_EUR-1][j].real();
            }
        }
    }

    ofstream fout;
    long length = snp_params.num_sample * snp_params.num_target_snp;
    // string output_filename = "score_total"
    string output_filename_AFR = "score_AFR";
    string output_filename_AMR = "score_AMR";
    string output_filename_EUR = "score_EUR";
    string suffix;
    if (snp_params.num_target_snp < 30000) {
        suffix = "_20k";
    } else if (snp_params.num_target_snp < 50000) {
        suffix = "_40k";
    } else if (snp_params.num_target_snp > 70000 && snp_params.num_target_snp < 90000) {
        suffix = "_80k";
    }
    // suffix += "_0430";
    // output_filename += suffix + ".csv";
    output_filename_AFR += suffix + ".csv";
    output_filename_AMR += suffix + ".csv";
    output_filename_EUR += suffix + ".csv";
    
    // fout.open(output_filename);
    // fout << "Subject ID,Target SNP,0,1,2\n";
    // for (int i = 0; i < num_sample_AFR * snp_params.num_target_snp; ++i) {
    //     fout << i / snp_params.num_target_snp + 1 << ",";
    //     fout << snp_data[i % snp_params.num_target_snp] << ",";    
    //     fout << score_AFR[3 * i] << "," << score_AFR[3 * i + 1] << "," << score_AFR[3 * i + 2] << '\n';
    // }
    // for (int i = 0; i < num_sample_AMR * snp_params.num_target_snp; ++i) {
    //     fout << i / snp_params.num_target_snp + 1 << ",";
    //     fout << snp_data[i % snp_params.num_target_snp] << ",";    
    //     fout << score_AMR[3 * i] << "," << score_AMR[3 * i + 1] << "," << score_AMR[3 * i + 2] << '\n';
    // }
    // for (int i = 0; i < num_sample_EUR * snp_params.num_target_snp; ++i) {
    //     fout << i / snp_params.num_target_snp + 1 << ",";
    //     fout << snp_data[i % snp_params.num_target_snp] << ",";    
    //     fout << score_EUR[3 * i] << "," << score_EUR[3 * i + 1] << "," << score_EUR[3 * i + 2] << '\n';
    // }
    // fout.close();

    fout.open(output_filename_AFR);
    fout << "Subject ID,Target SNP,0,1,2\n";
    for (int i = 0; i < num_sample_AFR * snp_params.num_target_snp; ++i) {
        fout << i / snp_params.num_target_snp + 1 << ",";
        fout << snp_data[i % snp_params.num_target_snp] << ",";    
        fout << score_AFR[3 * i] << "," << score_AFR[3 * i + 1] << "," << score_AFR[3 * i + 2] << '\n';
    }
    fout.close();

    fout.open(output_filename_AMR);
    fout << "Subject ID,Target SNP,0,1,2\n";
    for (int i = 0; i < num_sample_AMR * snp_params.num_target_snp; ++i) {
        fout << i / snp_params.num_target_snp + 1 << ",";
        fout << snp_data[i % snp_params.num_target_snp] << ",";    
        fout << score_AMR[3 * i] << "," << score_AMR[3 * i + 1] << "," << score_AMR[3 * i + 2] << '\n';
    }
    fout.close();

    fout.open(output_filename_EUR);
    fout << "Subject ID,Target SNP,0,1,2\n";
    for (int i = 0; i < num_sample_EUR * snp_params.num_target_snp; ++i) {
        fout << i / snp_params.num_target_snp + 1 << ",";
        fout << snp_data[i % snp_params.num_target_snp] << ",";    
        fout << score_EUR[3 * i] << "," << score_EUR[3 * i + 1] << "," << score_EUR[3 * i + 2] << '\n';
    }
    fout.close();

    // long num_sample_total = num_sample_AFR + num_sample_AMR + num_sample_EUR;

    // cout << "[Evaluator::computeScore_population] num_target_snp = " << snp_params.num_target_snp << endl;
    // cout << "[Evaluator::computeScore_population] num_sample_AFR = " << num_sample_AFR << endl;
    // cout << "[Evaluator::computeScore_population] AFR, length of score = " << num_sample_AFR * snp_params.num_target_snp << endl;
    // cout << "[Evaluator::computeScore_population] num_sample_AMR = " << num_sample_AMR << endl;
    // cout << "[Evaluator::computeScore_population] AMR, length of score = " << num_sample_AMR * snp_params.num_target_snp << endl;
    // cout << "[Evaluator::computeScore_population] num_sample_EUR = " << num_sample_EUR << endl;
    // cout << "[Evaluator::computeScore_population] EUR, length of score = " << num_sample_EUR * snp_params.num_target_snp << endl;
    // cout << "[Evaluator::computeScore_population] total length of score = " << num_sample_total * snp_params.num_target_snp << endl;
    cout << "[Evaluator::computeScore_population] Final scores are saved in " << output_filename_AFR 
         << ", " << output_filename_AMR << ", and " << output_filename_EUR << endl;
}


//Only using real part + encrypt columns of X in order
void Evaluator::encryptX_nosave(double** mat_X, SNP_Parameters snp_params, Ciphertext* ctxt_X)
{
    long nColumn = 3 * snp_params.num_tag_snp;
    // NTL_EXEC_RANGE(nColumn, first, last);
// #pragma omp parallel for
    for (int i = 0; i < nColumn; i++) {
        ofstream fout;
        Message msg(1 << logSlots);
        // for (int l = 0; l < (1 << logSlots) / 512; l++) {
        for (int j = 0; j < snp_params.num_sample; j++) {
            msg[j].real(mat_X[j][i]);
            // msg[j + 512 * l].imag(0);
        }
        // }
        // if(i==0){
        // cout << "[Evaluator::encryptX_nosave] msg" << endl;
        //     for(int j = 0; j < 20; j++){
        //         cout << msg[j] << ", ";
        //     }
        //     cout << endl; 
        // }

        encryptor.encrypt(msg, pk, ctxt_X[i]);
    }
    // NTL_EXEC_RANGE_END;  

    // cout << "[Evaluator::encryptX_nosave] decX" << endl;
    // Message dec_X;
    // decryptor.decrypt(ctxt_X[0], secret_key, dec_X);
    // for(int i = 0; i < 20; i++){
    //   cout << dec_X[i] << ", ";
    // }
    // cout<< endl;

    // cout << "[Evaluator::encryptX_nosave] SK" << endl;
    // for(int i = 0; i < (1 << logSlots); i++){
    //     cout << secret_key->getSxData()[i] << ", ";
    // }
    // cout << endl;

}



//Only using real part + encrypt columns of X in order
void Evaluator::encryptX_imag_nosave(double** mat_X, SNP_Parameters snp_params, Ciphertext* ctxt_X)
{
    long nColumn = (3 * snp_params.num_tag_snp + 1) / 2;
    long nslots = 1 << logSlots;
    Message* msg = new Message[nColumn];

    #pragma omp parallel for
    for(int i = 0; i < nColumn - 1; i++){
        msg[i].resize(nslots);
        for (int j = 0; j < snp_params.num_sample; j++) {
            msg[i][j].real(mat_X[j][2 * i]);
            msg[i][j].imag(-mat_X[j][2 * i + 1]);
        }
    }
    msg[nColumn-1].resize(nslots);
    if(snp_params.num_tag_snp % 2 == 0){
        for (int j = 0; j < snp_params.num_sample; j++) {
            msg[nColumn-1][j].real(mat_X[j][2 * nColumn - 2]);
            msg[nColumn-1][j].imag(-mat_X[j][2 * nColumn - 1]);
        }
    }
    else{
        for (int j = 0; j < snp_params.num_sample; j++) {
            msg[nColumn-1][j].real(mat_X[j][2 * nColumn - 2]);
        }      
    }


    #pragma omp parallel for
        for (int i = 0; i < nColumn; i++) {
        	encryptor.encrypt(msg[i], pk, ctxt_X[i]);
        }
    



        // cout << "[Evaluator::encryptX_nosave] decX" << endl;
        // Message dec_X;
        // decryptor.decrypt(ctxt_X[0], secret_key, dec_X);
        // for(int i = 0; i < 20; i++){
        //   cout << dec_X[i] << ", ";
        // }
        // cout<< endl;

        // cout << "[Evaluator::encryptX_nosave] SK" << endl;
        // for(int i = 0; i < (1 << logSlots); i++){
        //     cout << secret_key->getSxData()[i] << ", ";
        // }
        // cout << endl;
}

// TODO: Correct Check is NOT Done
void Evaluator::encryptX_imag_nosave_newdata(double** mat_X, SNP_Parameters snp_params, Ciphertext** ctxt_X)
{
    // cout << "tag_snp = " << snp_params.num_tag_snp << endl;
    long nColumn = (3 * snp_params.num_tag_snp + 1) / 2;
    // cout << "nColumn = " << nColumn << endl;
    long nslots = 1 << logSlots;
    long num_packing = (snp_params.num_sample - 1) / nslots + 1;
    // cout << "num_packing = " << num_packing << endl; 
    int remain = snp_params.num_sample % nslots;
    Message** msg = new Message*[nColumn];

#pragma omp parallel for
    for(int i = 0; i < nColumn - 1; i++){
        msg[i] = new Message[num_packing];
        for (int k = 0; k < num_packing - 1; k++) {
            msg[i][k].resize(nslots);
            for (int j = 0; j < nslots; j++) {
                msg[i][k][j].real(mat_X[nslots * k + j][2 * i]);
                msg[i][k][j].imag(-mat_X[nslots * k + j][2 * i + 1]);
            }
        }
        msg[i][num_packing - 1].resize(nslots);
        for (int j = 0; j < remain; j++) {
            msg[i][num_packing - 1][j].real(mat_X[snp_params.num_sample - remain + j][2 * i]);
            msg[i][num_packing - 1][j].imag(-mat_X[snp_params.num_sample - remain + j][2 * i + 1]);
        }
    }
    msg[nColumn-1] = new Message[num_packing];
    if(snp_params.num_tag_snp % 2 == 0){
        for (int k = 0; k < num_packing - 1; k++) {
            msg[nColumn-1][k].resize(nslots);
            for (int j = 0; j < nslots; j++) {
                msg[nColumn-1][k][j].real(mat_X[nslots * k + j][2 * nColumn - 2]);
                msg[nColumn-1][k][j].imag(-mat_X[nslots * k + j][2 * nColumn - 1]);
            }
        }
        msg[nColumn-1][num_packing-1].resize(nslots);
        for (int j = 0; j < remain; j++) {
            msg[nColumn-1][num_packing - 1][j].real(mat_X[snp_params.num_sample - remain + j][2 * nColumn - 2]);
            msg[nColumn-1][num_packing - 1][j].imag(-mat_X[snp_params.num_sample - remain + j][2 * nColumn - 1]);
        }
    } else {
        for (int k = 0; k < num_packing - 1; k++) {
            msg[nColumn-1][k].resize(nslots);
            for (int j = 0; j < nslots; j++) {
                msg[nColumn-1][k][j].real(mat_X[nslots * k + j][2 * nColumn - 2]);
            }
        }
        msg[nColumn-1][num_packing-1].resize(nslots);
        for (int j = 0; j < remain; j++) {
            msg[nColumn-1][num_packing - 1][j].real(mat_X[snp_params.num_sample - remain + j][2 * nColumn - 2]);
        }    
    }
    #pragma omp parallel for
        for (int i = 0; i < nColumn; i++) {
            ctxt_X[i] = new Ciphertext[num_packing];
            for (int k = 0; k < num_packing; k++) {
                encryptor.encrypt(msg[i][k], pk, ctxt_X[i][k]);
            }
        }
}


//Only using real part + encrypt columns of X in order
void Evaluator::compute_XW_nosave(Ciphertext* EncX, double*** W, SNP_Parameters snp_params, Ciphertext** EncXW)
{
#pragma omp parallel for
    for (int i = 0; i < snp_params.num_target_snp; i++) {
        EncXW[i] = new Ciphertext[snp_params.list_nSNP_target[i]];
    }
#pragma omp parallel for
    for (int i = 0; i < snp_params.num_target_snp; i++) {
        long nSNP_tag = snp_params.list_nSNP_tag[i];
        long nSNP_target = snp_params.list_nSNP_target[i];

        // TimeUtils timeutils;
        for (int k = 0; k < nSNP_tag; k++) {
            Ciphertext EncX_ik = EncX[3 * snp_params.list_x_start[i] + k];
            // timeutils.stop("fileopen");
            // timeutils.start("constmult");
            for (int j = 0; j < nSNP_target; j++) {
                if (k == 0) {
                    homEvaluator.constmultWithoutRescale(EncX_ik, W[i][k][j], 6, EncXW[i][j]);
                } else {
                    Ciphertext tmp;
                    // timeutils.start("constmult");
                    homEvaluator.constmultWithoutRescale(EncX_ik, W[i][k][j], 6, tmp);
                    // timeutils.stop("constmult");
                    // timeutils.start("add");
                    homEvaluator.add(EncXW[i][j], tmp, EncXW[i][j]);
                    // timeutils.stop("add");
                }
            }
            // timeutils.stop("constmult");
        }
   }
}



//Using both real and imaginary part + encrypt columns of X in order
//NOTE: Elements of list_x_start should be even numbers
//NOTE: Elements of list_x_end should be odd numbers
void Evaluator::compute_XW_imag_nosave(Ciphertext* EncX, double*** W, SNP_Parameters snp_params, Ciphertext** EncXW)
{
#pragma omp parallel for
    for (int i = 0; i < snp_params.num_target_snp; i++) {
        EncXW[i] = new Ciphertext[snp_params.list_nSNP_target[i]];
    }
#pragma omp parallel for
    for (int i = 0; i < snp_params.num_target_snp; i++) {
        long nSNP_tag = snp_params.list_nSNP_tag[i];
        long nSNP_target = snp_params.list_nSNP_target[i];

        for (int k = 0; k < nSNP_tag / 2; k++) {

            Ciphertext EncX_ik = EncX[3 * snp_params.list_x_start[i] / 2 + k];

            for (int j = 0; j < nSNP_target; j++) {

                if (k == 0) {
                    Ciphertext tmp;
                    homEvaluator.constmultWithoutRescale(EncX_ik, W[i][2 * k][j], 5, EncXW[i][j]);
                    homEvaluator.monomialmultWithoutRescale(EncX_ik, 1 << (logN - 1), W[i][2 * k + 1][j], 5, tmp);
                    homEvaluator.add(EncXW[i][j], tmp, EncXW[i][j]);
                } 
                else {
                    Ciphertext tmp1, tmp2;
                    homEvaluator.constmultWithoutRescale(EncX_ik, W[i][2 * k][j], 5, tmp1);
                    homEvaluator.monomialmultWithoutRescale(EncX_ik, 1 << (logN - 1), W[i][2 * k + 1][j], 5, tmp2);
                    homEvaluator.add(EncXW[i][j], tmp1, EncXW[i][j]);
                    homEvaluator.add(EncXW[i][j], tmp2, EncXW[i][j]);
                }
            }
        }
   }
}

void Evaluator::compute_XW_imag_nosave_newdata(Ciphertext** EncX, double*** W, SNP_Parameters snp_params, Ciphertext*** EncXW)
{
	long nslots = 1 << logSlots;
    long num_packing = (snp_params.num_sample - 1) / nslots + 1;
#pragma omp parallel for
    for (int i = 0; i < snp_params.num_target_snp; i++) {
        EncXW[i] = new Ciphertext*[snp_params.list_nSNP_target[i]];
        for (int j = 0; j < snp_params.list_nSNP_target[i]; j++) {
        	EncXW[i][j] = new Ciphertext[num_packing];
        }
    }
#pragma omp parallel for
    for (int i = 0; i < snp_params.num_target_snp; i++) {
        long nSNP_tag = snp_params.list_nSNP_tag[i];
        long nSNP_target = snp_params.list_nSNP_target[i];
        for (int l = 0; l < num_packing; l++) {
	        for (int k = 0; k < nSNP_tag / 2; k++) {

	            Ciphertext EncX_ik = EncX[3 * snp_params.list_x_start[i] / 2 + k][l];

	            for (int j = 0; j < nSNP_target; j++) {

	                if (k == 0) {
	                    Ciphertext tmp;
	                    homEvaluator.constmultWithoutRescale(EncX_ik, W[i][2 * k][j], 5, EncXW[i][j][l]);
	                    homEvaluator.monomialmultWithoutRescale(EncX_ik, 1 << (logN - 1), W[i][2 * k + 1][j], 5, tmp);
	                    homEvaluator.add(EncXW[i][j][l], tmp, EncXW[i][j][l]);
	                } 
	                else {
	                    Ciphertext tmp1, tmp2;
	                    homEvaluator.constmultWithoutRescale(EncX_ik, W[i][2 * k][j], 5, tmp1);
	                    homEvaluator.monomialmultWithoutRescale(EncX_ik, 1 << (logN - 1), W[i][2 * k + 1][j], 5, tmp2);
	                    homEvaluator.add(EncXW[i][j][l], tmp1, EncXW[i][j][l]);
	                    homEvaluator.add(EncXW[i][j][l], tmp2, EncXW[i][j][l]);
	                }
	            }
	        }
	    }
	}
}


void Evaluator::decryptXW_nosave(Message** dec_XW, Ciphertext** EncXW, SNP_Parameters snp_params)
{
    // NTL_EXEC_RANGE(snp_params.num_target_snp, first, last);
#pragma omp parallel for
    for (int i = 0; i < snp_params.num_target_snp; i++) {
	for (int j = 0; j < snp_params.list_nSNP_target[i]; j++) {
            decryptor.decrypt(EncXW[i][j], secret_key, dec_XW[i][j]);
        }
    }
    // NTL_EXEC_RANGE_END;
}

void Evaluator::decryptXW_nosave_newdata(Message*** dec_XW, Ciphertext*** EncXW, SNP_Parameters snp_params)
{
	long nslots = 1 << logSlots;
    long degree = context.getDegree();
    long num_packing = (snp_params.num_sample - 1) / nslots + 1;
    
// #pragma omp parallel for
//     for (int i = 0; i < snp_params.num_target_snp; i++) {
//     	for (int j = 0; j < snp_params.list_nSNP_target[i]; j++) {
//         	dec_XW[i][j] = new Message[num_packing];
//             for (int l = 0; l < num_packing; l++) {
//             dec_XW[i][j][l].resize(nslots);
// 		    decryptor.decrypt(EncXW[i][j][l], secret_key, dec_XW[i][j][l]);
//     	}
//         }
//     }


    //START();
    int nthreads = 16;

    PrimePoly* pp1 = new PrimePoly[nthreads];
    PrimePoly* pp2 = new PrimePoly[nthreads];
    LargePoly* lp = new LargePoly[nthreads];
    Plaintext* ptxt = new Plaintext[nthreads];
    secret_key.save("sk.txt");
    SecretKey* sk = new SecretKey[nthreads];

    for (int l = 0; l < nthreads; l++) {
        pp1[l].resize(degree);
        pp2[l].resize(degree);
        lp[l].resize(degree);
        ptxt[l].mx__.resize(degree);
        sk[l].load("sk.txt");
    }
    
    //END();
    //PRINTTIME("alloc finished", 1);
    //START();
#pragma omp parallel for
    for (int i = 0; i < snp_params.num_target_snp; i++) {
    	for (int j = 0; j < snp_params.list_nSNP_target[i]; j++) {
            dec_XW[i][j] = new Message[num_packing];
            for (int l = 0; l < num_packing; l++) {
                dec_XW[i][j][l].resize(nslots);
            }
        }
    }

    //END();
    //PRINTTIME("dec_XW finished", 1);
	cout << "Memory allocation finished..." << endl;
    // START();
    // cout << "Dgree = " << degree << endl;
    // PrimePoly*** pp1 = new PrimePoly**[snp_params.num_target_snp];
    // PrimePoly*** pp2 = new PrimePoly**[snp_params.num_target_snp];
    // LargePoly*** lp = new LargePoly**[snp_params.num_target_snp];
    // Plaintext*** ptxt = new Plaintext**[snp_params.num_target_snp];
    
    // for (int i = 0; i < snp_params.num_target_snp; i++) {
    //     pp1[i] = new PrimePoly*[snp_params.list_nSNP_target[i]];
    //     pp2[i] = new PrimePoly*[snp_params.list_nSNP_target[i]];
    //     lp[i] = new LargePoly*[snp_params.list_nSNP_target[i]];
    //     ptxt[i] = new Plaintext*[snp_params.list_nSNP_target[i]];
    // 	for (int j = 0; j < snp_params.list_nSNP_target[i]; j++) {
    //         pp1[i][j] = new PrimePoly[num_packing];
    //         pp2[i][j] = new PrimePoly[num_packing];
    //         lp[i][j] = new LargePoly[num_packing];
    //         ptxt[i][j] = new Plaintext[num_packing];
    //     	dec_XW[i][j] = new Message[num_packing];
    //         for (int l = 0; l < num_packing; l++) {
    //             pp1[i][j][l].resize(degree);
    //             pp2[i][j][l].resize(degree);
    //             lp[i][j][l].resize(degree);
    //             dec_XW[i][j][l].resize(nslots);
    //             ptxt[i][j][l].mx__.resize(degree);
    //         }
    //     }
    // }

    // END();
    // PRINTTIME("alloc finished", 1);



    START();
#pragma omp parallel for
    for (int i = 0; i < snp_params.num_target_snp; i++) {
    	for (int j = 0; j < snp_params.list_nSNP_target[i]; j++) {
            for (int l = 0; l < num_packing; l++) {
                long cthread = omp_get_thread_num();
                decryptor.decrypt_alloc(EncXW[i][j][l], sk[cthread], pp1[cthread], pp2[cthread], lp[cthread], ptxt[cthread], dec_XW[i][j][l]);
    	    }
        }
    }
    
    END();
    PRINTTIME("Decypt XW", 1);
}


void Evaluator::readW_nosave(long datatype, long data_version, SNP_Parameters snp_params, long population, double*** W)
{
    // NTL_EXEC_RANGE(snp_params.num_target_snp, first, last);
    string population_names[3] = {"AFR", "AMR", "EUR"};
#pragma omp parallel for
    for (int i = 0; i < snp_params.num_target_snp; i++) {
        string input_filename;
        if (datatype != 0) {
            input_filename = "../DNNmodels_New" + to_string(datatype) + "k/W_New" + to_string(i) + ".csv";
        } else if (population != 0) {
            input_filename = "../" + population_names[population - 1] + "_DNNmodels/DNNmodels_" + to_string(snp_params.window) + "_c/W_New" + to_string(i + 1666) + ".csv";
        } else {
            input_filename = "../Total_DNNmodels/DNNmodels_" + to_string(snp_params.window) + "_c/W_New" + to_string(i + 1666) + ".csv";
            // input_filename = "../../../idash2019_final/enc/Total_DNNmodels/DNNmodels_48_epoc50_0.2/W_New" + to_string(i + 1666) + ".csv";
            // input_filename = "../../../idash2019_final/enc/EXP_DNNmodels/DNNmodels_" + to_string(snp_params.window) + "_epoc40/W_New" + to_string(i + 1666) + ".csv";
        }
        if (i == 0) cout << "readW from " << input_filename << endl;
	    FileReader::readMatrix(input_filename, W[i], 0, ',', 0);
    }
    // NTL_EXEC_RANGE_END;
}
