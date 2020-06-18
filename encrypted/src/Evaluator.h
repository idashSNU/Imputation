#ifndef EVALUATOR_H_
#define EVALUATOR_H_

#include "FileReader.h"
#include "HEAAN.h"
#include "SNP_Parameters.h"
#include "TimeUtils.h"
#include <fstream>
#include <omp.h>

#include <iostream>
#include <sys/stat.h>
#include <sys/types.h>
#include <vector>
#include <ctime>
#include <chrono>
#include <omp.h>

using namespace std;
using namespace heaan;

class Evaluator {
public:
    long logp;
    long logQ;
    long logSlots;
    long logN;

    long num_packing;

    Parameters params;
    Context context;

    // secret_key is for debugging
    // Should be removed before submission!
    SecretKey secret_key;

    PublicKeyPack public_key_pack;
    PublicKey pk;
    Encryptor encryptor;
    Decryptor decryptor;
    HomEvaluator homEvaluator;

    Evaluator(long logN, long logp, long logQ);


    void encryptX_nosave(double** mat_X, SNP_Parameters snp_params, Ciphertext* ctxt_X);
    void compute_XW_nosave(Ciphertext* EncX, double*** W, SNP_Parameters snp_params, Ciphertext** EncXW);
    void encryptX_imag_nosave(double** mat_X, SNP_Parameters snp_params, Ciphertext* ctxt_X);
    void compute_XW_imag_nosave(Ciphertext* EncX, double*** W, SNP_Parameters snp_params, Ciphertext** EncXW);
    void encryptX_imag_nosave_newdata(double** mat_X, SNP_Parameters snp_params, Ciphertext** ctxt_X);
    void compute_XW_imag_nosave_newdata(Ciphertext** EncX, double*** W, SNP_Parameters snp_params, Ciphertext*** EncXW);
    void decryptXW_nosave(Message** dec_XW, Ciphertext** EncXW, SNP_Parameters snp_params);
    void decryptXW_nosave_newdata(Message*** dec_XW, Ciphertext*** EncXW, SNP_Parameters snp_params);
    void readW_nosave(long datatype, long data_version, SNP_Parameters snp_params, long population, double*** W);
    void computeScore(Message** dec_XW, string* snp_data, long data_version, SNP_Parameters snp_params);
    void computeScore_newdata(Message*** dec_XW, string* snp_data, long data_version, SNP_Parameters snp_params);
    void computeScore_population(Message*** dec_XW, Message*** dec_XW2, Message*** dec_XW3, string* snp_data, SNP_Parameters snp_params);
};

#endif
