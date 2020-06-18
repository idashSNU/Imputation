#include "PrintUtils.h"

void PrintUtils::averageDifference(double* a1, std::complex<double>* a2, long n)
{
    double avg = 0.;
    for (long i = 0; i < n; i++) {
        avg += abs(a1[i] - a2[i].real());
    }
    avg /= n;
    std::cout << "log2(avg of error) = " << log2(avg) << std::endl;
}

void PrintUtils::averageDifference(std::complex<double>* a1, std::complex<double>* a2, long n)
{
    double avg = 0.;
    for (long i = 0; i < n; i++) {
        avg += abs(a1[i].real() - a2[i].real());
    }
    avg /= n;
    std::cout << "log2(avg of error) = " << log2(avg) << std::endl;
}

void PrintUtils::averageDifference(std::vector<complex<double>> a1, std::vector<complex<double>> a2, long n) {
    double avg_real = 0.;
    double avg_imag = 0.;
    for (long i = 0; i < n; i++) {
        avg_real += abs(a1[i].real() - a2[i].real());
        avg_imag += abs(a1[i].imag() - a2[i].imag());
    }
    avg_real /= n;
    avg_imag /= n;
    std::cout << "log2(avg of error) = " << log2(avg_real) << ", " << log2(avg_imag) << std::endl;
}

void PrintUtils::printSingleArray(std::string str, double* array, long n)
{
    for (int i = 0; i < n; i++) {
        std::cout << str << "[" << i << "] : " << array[i] << std::endl;
    }
}

void PrintUtils::printSingleArray(std::string str, complex<double>* array, long n)
{
    for (int i = 0; i < n; i++) {
        std::cout << str << "[" << i << "] = " << array[i] << std::endl;
    }
}

void PrintUtils::printSingleArraySmall(std::string str, double* array, long n)
{
    if (n <= 1024) {
        PrintUtils::printSingleArray(str, array, n);
    } else {
        for (int i = 0; i < n; i++) {
            if (i % 1000 == 0) {
                std::cout << str << "[" << i << "] = " << array[i] << std::endl;
            }
        }
    }
}

void PrintUtils::printSingleArraySmall(std::string str, complex<double>* array, long n)
{
    if (n <= 1024) {
        PrintUtils::printSingleArray(str, array, n);
    } else {
        for (int i = 0; i < n; i++) {
            if (i % 1000 == 0) {
                std::cout << str << "[" << i << "] = " << array[i] << std::endl;
            }
        }
    }
}

void PrintUtils::printSingleMatrix(std::string str, double** matrix, long row, long col)
{
    for (int i = 0; i < row; i++) {
        std::cout << str << "[" << i << "] = [";
        for (int j = 0; j < col; j++) {
            std::cout << matrix[i][j] << ", ";
        }
        std::cout << "]" << std::endl;
    }
}

void PrintUtils::printSingleMatrix(std::string str, complex<double>** matrix, long row, long col)
{
    for (int i = 0; i < row; i++) {
        std::cout << str << "[" << i << "] = [";
        for (int j = 0; j < col; j++) {
            std::cout << matrix[i][j] << ", ";
        }
        std::cout << "]" << std::endl;
    }
}

void PrintUtils::printArrays(double* a1, std::complex<double>* a2, long n)
{
    for (int i = 0; i < n; i++) {
        std::cout << i << " : " << a1[i] << " // " << a2[i] << std::endl;
    }
}

void PrintUtils::printFewArrays(double* a1, std::complex<double>* a2, long n)
{
    for (int i = 0; i < n; i++) {
        if (n < 1000 || i % 1000 == 0) {
            std::cout << i << " : " << a1[i] << " // " << a2[i] << std::endl;
        }
    }
}

void PrintUtils::printFewArrays(std::complex<double>* a1, std::complex<double>* a2, long n)
{
    for (int i = 0; i < n; i++) {
        if (n < 1000 || i % 1000 == 0) {
            std::cout << i << " : " << a1[i] << " // " << a2[i] << std::endl;
        }
    }
}

void PrintUtils::printArraysWithDataNum(double* a1, double* a2, long n, long logDataNum, long colNum)
{
    long dataNum = 1 << logDataNum;
    for (int i = 0; i < n; i++) {
        if (i % dataNum == 0) {
            std::cout << "-----------------------------" << endl;
        }
        if (i % dataNum == colNum) {
            std::cout << "<<";
        }
        std::cout << i << " : " << a1[i] << " // " << a2[i];
        if (i % dataNum == colNum) {
            std::cout << ">>";
        }
        std::cout << std::endl;
    }
}

void PrintUtils::printArraysWithDataNum(double* a1, std::complex<double>* a2, long n, long logDataNum, long colNum)
{
    long dataNum = 1 << logDataNum;
    for (int i = 0; i < n; i++) {
        if (i % dataNum == 0) {
            std::cout << "-----------------------------" << endl;
        }
        if (i % dataNum == colNum) {
            std::cout << "<<";
        }
        std::cout << i << " : " << a1[i] << " // " << a2[i].real();
        if (i % dataNum == colNum) {
            std::cout << ">>";
        }
        std::cout << std::endl;
    }
}

void PrintUtils::printArrays(std::complex<double>* a1, std::complex<double>* a2, long n)
{
    for (int i = 0; i < n; i++) {
        std::cout << i << " : " << a1[i] << " // " << a2[i] << std::endl;
    }
}

void PrintUtils::printArraysWithDataNum(std::complex<double>* a1, std::complex<double>* a2, long n, long logDataNum)
{
    long dataNum = 1 << logDataNum;
    for (int i = 0; i < n; i++) {
        if (i % dataNum == 0) {
            std::cout << "<<";
        }
        std::cout << i << " : " << a1[i].real() << " // " << a2[i].real();
        if (i % dataNum == 0) {
            std::cout << ">>";
        }
        std::cout << std::endl;
    }
}

void PrintUtils::nprint(std::string str, bool isPrint)
{
    if (isPrint)
        std::cout << str << std::endl;
}

// *********** *********** ***********
// *********** *********** ***********
// *********** *********** ***********

void PrintUtils::printMessage(Message msg)
{
    if (msg.size() > 5) {
        cout << "[";
        for (long i = 0; i < 4; i++) {
            cout << "(" << msg[i].real() << "," << msg[i].imag() << "), ";
        }
        cout << "... ,(" << msg[msg.size() - 1].real() << "," << msg[msg.size() - 1].imag() << ")]" << endl;
    } else {
        cout << "[";
        for (long i = 0; i < msg.size() - 1; i++) {
            cout << "(" << msg[i].real() << "," << msg[i].imag() << "), ";
        }
        cout << "(" << msg[msg.size() - 1].real() << "," << msg[msg.size() - 1].imag() << ")]" << endl;
    }
}