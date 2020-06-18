#ifndef FILEREADER_H_
#define FILEREADER_H_

// #include <direct.h>
#include <sys/stat.h>
#include <complex>
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>
#include "SNP_Parameters.h"

using namespace std;

class FileReader {
public:
    static void readMatrix(const string input_filename, double** mat, long = 0, char = ',', long = 0);
    static void readMatrix(const string input_filename, string** mat, long = 0, char = ',', long = 0);

    static void readMatrix(const string input_filename, double** mat, long& col_num, long = 0, char = ',', long = 0);

    static void readMatrixWithData(const string input_filename, double** mat, vector<double>& snp_data, long& col_num, long column_start, char delimiter, long remove_first_col, long snp_col);

    static void readAndEncodeX(const string input_filename, double** matX, long& col_num, long data_num, long snp_num, bool test);

    static void transpose(double** mat_trans, double** mat, double row_trans, double col_trans);

    static void encode_snp(double** mat_encode, double** mat, double row, double col);

    static vector<string> split(string str, char delimiter);

    static void read_x_start(SNP_Parameters& snp_params);
};

#endif