#include "SNP_Parameters.h"
#include "iostream"


// SNP_Parameters::SNP_Parameters(long _datatype, long _data_version)
// {
//     window = 64; // modify here
//     SNP_Parameters(_datatype, _data_version, window, );
//     std::cout << "2snp_params.num_target_snp = " << num_target_snp << std::endl;
// }

SNP_Parameters::SNP_Parameters(long _datatype, long _data_version, long _window, long _num_target_snp_k)
    : datatype(_datatype)
    , data_version(_data_version)
{
    if (datatype == 10) {
        num_tag_snp = 1045;
        num_target_snp = 500;
        list_x_start = list_x10k_start_3;
        list_x_end = new long[500];
        for(int i = 0; i < 500; i++){
            list_x_end[i] = list_x_start[i] + 25;
        }
        list_y_start = new long[500];
        for(int i = 0; i < 500; i++){
            list_y_start[i] = i;
        }
        list_y_end = list_y_start; 
    } 

    else if (datatype == 1) {
        num_tag_snp = 9746;
        num_target_snp = 500;
        list_x_start = list_x1k_start_1;
        list_x_end = new long[500];
        for(int i = 0; i < 500; i++){
            list_x_end[i] = list_x_start[i] + 69;
        }
        list_y_start = new long[500];
        for(int i = 0; i < 500; i++){
            list_y_start[i] = i;
        }
        list_y_end = list_y_start; 
    }

    // else if (datatype == 0) {
    num_tag_snp = 16184;
    num_target_snp_k = _num_target_snp_k;
    if (num_target_snp_k == 20) {
        num_target_snp = 20000;
    } else if (num_target_snp_k == 40) {
        num_target_snp = 40000;
    } else if (num_target_snp_k == 80) {
        num_target_snp = 80882;
    } else if (num_target_snp_k == 117) {
        num_target_snp = 117904;
    } else {
    std::cerr << "ERROR :: No Matching num_target_snp_k" << std::endl;
    }
    window = _window ;  
    std::cout << "num_target_snp = " << num_target_snp << ", window = " << window << std::endl;
    // list_x_start = new long[num_target_snp];
    // list_x_end = new long[num_target_snp];
    // for(int i = 0; i < num_target_snp; i++){
    //     list_x_start[i] = list_x_start_full[i + 1666];
    //     list_x_end[i] = list_x_start[i] + window - 1; 
    // }
    // list_y_start = new long[num_target_snp];
    // for(int i = 0; i < num_target_snp; i++){
    //     list_y_start[i] = i + 1666;
    // }
    // list_y_end = list_y_start; 
    // }
    // list_nSNP_target = new long[num_target_snp];
    // list_nSNP_tag = new long[num_target_snp];

    // for (int i = 0; i < num_target_snp; i++) {
    //     list_nSNP_target[i] = (list_y_end[i] - list_y_start[i] + 1) * 3;
    //     list_nSNP_tag[i] = (list_x_end[i] - list_x_start[i] + 1) * 3;
    // }
    std::cout << "snp_params.num_target_snp = " << num_target_snp << std::endl;
    
}

void SNP_Parameters::setLists_newdata()
{
    list_x_start = new long[num_target_snp];
    list_x_end = new long[num_target_snp];
    for(int i = 0; i < num_target_snp; i++){
        if (num_target_snp_k <= 80) {
            list_x_start[i] = list_x_start_full[i + 1666];
        } else {
            list_x_start[i] = list_x_start_full[i];
        }
        list_x_end[i] = list_x_start[i] + window - 1; 
    }
    list_y_start = new long[num_target_snp];
    for(int i = 0; i < num_target_snp; i++){
        if (num_target_snp_k <= 80) {
            list_y_start[i] = i + 1666;
        } else {
            list_y_start[i] = i;
        }
    }
    
    list_y_end = list_y_start; 

    list_nSNP_target = new long[num_target_snp];
    list_nSNP_tag = new long[num_target_snp];
    for (int i = 0; i < num_target_snp; i++) {
        list_nSNP_target[i] = (list_y_end[i] - list_y_start[i] + 1) * 3;
        list_nSNP_tag[i] = (list_x_end[i] - list_x_start[i] + 1) * 3;
    }
}