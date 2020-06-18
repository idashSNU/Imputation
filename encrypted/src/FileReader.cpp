#include "FileReader.h"

void FileReader::readMatrix(const string input_filename, double** mat, long column_start, char delimiter, long remove_first_col)
{
    ifstream openFile(input_filename.data());
    vector<vector<double>> file_data;
    long count_row = 0;
    long count_col = 0;
    if (openFile.is_open()) {
        string line, temp;
        while (getline(openFile, line)) {
            vector<double> xline;
            vector<string> lineBreak = split(line, delimiter);
            count_col = lineBreak.size() - remove_first_col - column_start;
            for (long i = remove_first_col + column_start; i < lineBreak.size(); ++i) {
                if (lineBreak[i].compare("NaN") == 0) {
                    xline.push_back(-1.);
                } else {
                    xline.push_back(stod(lineBreak[i]));
                }
            }
            file_data.push_back(xline);
            count_row++;
        }

        for (int i = 0; i < count_row; i++) {
            for (int j = 0; j < count_col; j++) {
                mat[i][j] = file_data[i][j];
            }
        }
    } else {
        cout << "Error: cannot read file" << endl;
    }
}

void FileReader::readMatrix(const string input_filename, string** mat, long column_start, char delimiter, long remove_first_col)
{
    ifstream openFile(input_filename.data());
    vector<vector<string>> file_data;
    long count_row = 0;
    long count_col = 0;
    if (openFile.is_open()) {
        string line, temp;
        while (getline(openFile, line)) {
            vector<string> xline;
            vector<string> lineBreak = split(line, delimiter);
            count_col = lineBreak.size() - remove_first_col - column_start;
            for (long i = remove_first_col + column_start; i < lineBreak.size(); ++i) {
                xline.push_back(lineBreak[i]);
            }
            file_data.push_back(xline);
            count_row++;
        }

        for (int i = 0; i < count_row; i++) {
            for (int j = 0; j < count_col; j++) {
                mat[i][j] = file_data[i][j];
            }
        }
    } else {
        cerr << "Error: cannot read file" << endl;
    }
}

void FileReader::readMatrixWithData(const string input_filename, double** mat, vector<double>& snp_data, long& col_num, long column_start, char delimiter, long remove_first_col, long snp_col)
{
    ifstream openFile(input_filename.data());
    vector<vector<double>> file_data;
    long count_row = 0;
    long count_col = 0;
    if (openFile.is_open()) {
        string line, temp;
        while (getline(openFile, line)) {
            vector<double> xline;
            vector<string> lineBreak = split(line, delimiter);
            count_col = lineBreak.size() - remove_first_col - column_start;
            for (long i = remove_first_col + column_start; i < lineBreak.size(); ++i) {
                if (lineBreak[i].compare("NaN") == 0) {
                    xline.push_back(-1.);
                } else {
                    xline.push_back(stod(lineBreak[i]));
                }
            }
            snp_data.push_back(stod(lineBreak[snp_col]));
            file_data.push_back(xline);
            count_row++;
        }

        col_num = count_col;

        for (int i = 0; i < count_row; i++) {
            for (int j = 0; j < count_col; j++) {
                mat[i][j] = file_data[i][j];
            }
        }
    } else {
        cerr << "Error: cannot read file" << endl;
    }
}

void FileReader::readMatrix(const string input_filename, double** mat, long& col_num, long column_start, char delimiter, long remove_first_col)
{
    ifstream openFile(input_filename.data());
    vector<vector<double>> file_data;
    long count_row = 0;
    long count_col = 0;
    if (openFile.is_open()) {
        string line, temp;
        while (getline(openFile, line)) {
            vector<double> xline;
            vector<string> lineBreak = split(line, delimiter);
            count_col = lineBreak.size() - remove_first_col - column_start;
            for (long i = remove_first_col + column_start; i < lineBreak.size(); ++i) {
                if (lineBreak[i].compare("NaN") == 0) {
                    xline.push_back(-1.);
                } else {
                    xline.push_back(stod(lineBreak[i]));
                }
            }
            file_data.push_back(xline);
            count_row++;
        }
        
        col_num = count_col;
        for (int i = 0; i < count_row; i++) {
            for (int j = 0; j < count_col; j++) {
                mat[i][j] = file_data[i][j];
            }
        }
    } else {
        cerr << "Error: cannot read file" << endl;
    }
}

void FileReader::readAndEncodeX(const string input_filename, double** matX, long& col_num, long data_num, long snp_num, bool test)
{
    double** mat = new double*[snp_num]();
    for (int i = 0; i < snp_num; i++) {
        mat[i] = new double[data_num]();
    }
    
    if (test) {
        // FileReader::readMatrixWithData(input_filename, mat, snp_data_vec, col_num, 2504 - data_num, '\t', 4, 2);
        FileReader::readMatrix(input_filename, mat, col_num, 2504 - data_num, '\t', 4);
    } else {
        // FileReader::readMatrixWithData(input_filename, mat, snp_data_vec, col_num, 0, '\t', 4, 2);
        FileReader::readMatrix(input_filename, mat, col_num, 0, '\t', 4);
    }

    double** mat_trans = new double*[data_num]();
    for (int i = 0; i < data_num; i++) {
        mat_trans[i] = new double[snp_num]();
    }

    FileReader::transpose(mat_trans, mat, data_num, snp_num);

    FileReader::encode_snp(matX, mat_trans, data_num, snp_num);
}

void FileReader::encode_snp(double** mat_encode, double** mat, double row_origin, double col_origin)
{
    for (int i = 0; i < row_origin; i++) {
        for (int j = 0; j < col_origin; j++) {
            double* encode = new double[3];
            if (mat[i][j] > -0.01 && mat[i][j] < 0.01) {
                encode[0] = 1.;
                encode[1] = 0;
                encode[2] = 0;
            } else if (mat[i][j] > 0.99 && mat[i][j] < 1.01) {
                encode[0] = 0;
                encode[1] = 1.;
                encode[2] = 0;
            } else if (mat[i][j] > 1.99 && mat[i][j] < 2.01) {
                encode[0] = 0;
                encode[1] = 0;
                encode[2] = 1.;
            } else if (mat[i][j] < 0) {
                // case for NaN
                encode[0] = 0;
                encode[1] = 0;
                encode[2] = 0;
            } else {
                cout << " !! ERROR in FileReader::encode_snp !!" << endl;
            }
            for (int k = 0; k < 3; k++) {
                mat_encode[i][3 * j + k] = encode[k];
            }
        }
    }
}

void FileReader::transpose(double** mat_trans, double** mat, double row_trans, double col_trans)
{
    for (int i = 0; i < row_trans; i++) {
        for (int j = 0; j < col_trans; j++) {
            mat_trans[i][j] = mat[j][i];
        }
    }
}

vector<string> FileReader::split(string str, char delimiter)
{
    vector<string> internal;
    stringstream ss(str); // Turn the string into a stream.
    string tok;

    while (getline(ss, tok, delimiter)) {
        internal.push_back(tok);
    }

    return internal;
}

void FileReader::read_x_start(SNP_Parameters& snp_params)
{
    string filename = "../window_location/window" + to_string(snp_params.window) + ".txt";
    cout << "[SNP_Parameters::read_x_start] read window " << snp_params.window << " from " << filename << endl;
    ifstream openFile(filename.data());
    long count_row = 0;
    long count_col = 0;
    if (openFile.is_open()) {
        string line, temp;
        vector<long> xline;
        while (getline(openFile, line)) {
            vector<string> lineBreak = split(line, ',');
            count_col = lineBreak.size();
            for (long i = 0; i < lineBreak.size(); ++i) {
                xline.push_back(stoi(lineBreak[i]));
            }
            count_row++;
        }
        cout << "[SNP_Parameters::read_x_start] check difference" << endl;
        for (int j = 0; j < count_col; j++) {
            // long origin = snp_params.list_x_start_full[i];
            // if (origin != xline[i]) {
            //     cout << "diff in i = " << i << ", " << origin << ", " << xline[i] << endl;
            // }
            snp_params.list_x_start_full[j] = xline[j];
        }
        cout << "number of column = " << count_col << endl;
    } else {
        cout << "[SNP_Parameters::read_x_start] Error: cannot read file for window " << snp_params.window << endl;
    }
}
