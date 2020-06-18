#ifndef PRINTUTILS_H_
#define PRINTUTILS_H_

#include "iostream"
#include <complex>
#include <math.h>
#include <string>
#include "HEAAN.h"

using namespace std;
using namespace heaan;

class PrintUtils {
public:
    static void averageDifference(double* a1, std::complex<double>* a2, long n);
    static void averageDifference(std::complex<double>* a1, std::complex<double>* a2, long n);
    static void averageDifference(std::vector<complex<double>> a1, std::vector<complex<double>> a2, long n);

    static void printSingleArray(std::string str, double* array, long n);
    static void printSingleArray(std::string str, complex<double>* array, long n);
    static void printSingleArraySmall(std::string str, double* array, long n);
    static void printSingleArraySmall(std::string str, complex<double>* array, long n);
    static void printSingleMatrix(std::string str, double** matrix, long row, long col);
    static void printSingleMatrix(std::string str, complex<double>** matrix, long row, long col);
    static void printArrays(double* a1, std::complex<double>* a2, long n);
    static void printArrays(std::complex<double>* a1, std::complex<double>* a2, long n);
    static void printFewArrays(double* a1, std::complex<double>* a2, long n);
    static void printFewArrays(std::complex<double>* a1, std::complex<double>* a2, long n);
    static void printArraysWithDataNum(double* a1, double* a2, long n, long logDataNum, long colNum);
    static void printArraysWithDataNum(double* a1, std::complex<double>* a2, long n, long logDataNum, long colNum);
    static void printArraysWithDataNum(std::complex<double>* a1, std::complex<double>* a2, long n, long logDataNum);

    static void nprint(std::string str, bool isPrint);

    // ************************
    static void printMessage(Message msg);
};

;

#endif // !PRINTUTILS_H_